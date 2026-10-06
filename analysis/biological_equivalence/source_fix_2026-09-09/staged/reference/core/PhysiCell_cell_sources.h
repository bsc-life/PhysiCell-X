#ifndef __PhysiCell_cell_sources_h__
#define __PhysiCell_cell_sources_h__

#include "PhysiCell_cell.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <stdexcept>

namespace PhysiCell {
namespace cell_source_detail {

// Optional per-rank accounting of source updates, measured before the next
// diffusion step. Interval sums avoid writing a row for every diffusion tick.
class Budget_Log {
    std::ofstream output;
    unsigned long every = 1000, steps = 0;
    double start = 0., end = 0.;
    std::vector<double> net, requested, applied, max_error;
public:
    explicit Budget_Log(int rank) {
        const char* prefix = std::getenv("PHYSICELL_SUBSTRATE_BUDGET_PREFIX");
        if (!prefix || !*prefix) return;
        const char* interval = std::getenv("PHYSICELL_SUBSTRATE_BUDGET_EVERY");
        if (interval) every = std::max(1UL, std::strtoul(interval, NULL, 10));
        const std::string path = std::string(prefix) + ".rank" + std::to_string(rank) + ".csv";
        output.open(path);
        if (!output) throw std::runtime_error("Cannot open substrate budget: " + path);
        output << "interval_start_min,interval_end_min,steps,substrate_index,"
                  "requested_net_export,requested_total,applied_total,max_abs_step_error\n";
        output << std::setprecision(17);
    }
    bool enabled() const { return output.is_open(); }
    void flush() {
        if (!steps) return;
        for (size_t d = 0; d < requested.size(); ++d)
            output << start << ',' << end << ',' << steps << ',' << d << ','
                   << net[d] << ',' << requested[d] << ',' << applied[d] << ',' << max_error[d] << '\n';
        output.flush();
        std::fill(net.begin(), net.end(), 0.);
        std::fill(requested.begin(), requested.end(), 0.);
        std::fill(applied.begin(), applied.end(), 0.);
        std::fill(max_error.begin(), max_error.end(), 0.);
        steps = 0;
    }
    void record(double time, double dt, const std::vector<double>& n,
                const std::vector<double>& r, const std::vector<double>& a) {
        if (!steps) start = time;
        end = time + dt;
        const size_t densities = std::max(net.size(), r.size());
        net.resize(densities, 0.); requested.resize(densities, 0.);
        applied.resize(densities, 0.); max_error.resize(densities, 0.);
        for (size_t d = 0; d < r.size(); ++d) {
            net[d] += n[d]; requested[d] += r[d]; applied[d] += a[d];
            max_error[d] = std::max(max_error[d], std::abs(a[d] - r[d]));
        }
        if (++steps == every) flush();
    }
    ~Budget_Log() { flush(); }
};

struct Group {
    size_t begin, end;
    std::vector<double> net, requested, applied;
};

} // namespace cell_source_detail

// One writer per substrate voxel. Stable cell-ID ordering also removes
// dependence on OpenMP scheduling for sequential implicit secretion/uptake.
inline void advance_cell_sources(const std::vector<Cell*>& all, double time,
                                 double dt, int rank = 0) {
    static cell_source_detail::Budget_Log budget(rank);
    const bool audit = budget.enabled();
    std::vector<Cell*> cells;
    cells.reserve(all.size());
    for (Cell* cell : all)
        if (!cell->is_out_of_domain && cell->is_active && cell->get_current_voxel_index() >= 0)
            cells.push_back(cell);
    std::sort(cells.begin(), cells.end(), [](Cell* a, Cell* b) {
        if (a->get_microenvironment() != b->get_microenvironment())
            return std::less<BioFVM::Microenvironment*>()(a->get_microenvironment(), b->get_microenvironment());
        if (a->get_current_voxel_index() != b->get_current_voxel_index())
            return a->get_current_voxel_index() < b->get_current_voxel_index();
        return a->ID < b->ID;
    });
    std::vector<cell_source_detail::Group> groups;
    for (size_t begin = 0; begin < cells.size();) {
        size_t end = begin + 1;
        while (end < cells.size() && cells[end]->get_microenvironment() == cells[begin]->get_microenvironment()
               && cells[end]->get_current_voxel_index() == cells[begin]->get_current_voxel_index()) ++end;
        cell_source_detail::Group group;
        group.begin = begin; group.end = end;
        groups.push_back(group);
        begin = end;
    }

    #pragma omp parallel for schedule(static)
    for (size_t g = 0; g < groups.size(); ++g) {
        auto& group = groups[g];
        Cell* first = cells[group.begin];
        auto* environment = first->get_microenvironment();
        const int voxel = first->get_current_voxel_index();
        const int densities = environment->number_of_densities();
        const double voxel_volume = environment->voxels(voxel).volume;
        std::vector<double> before, predicted;
        if (audit) {
            before.resize(densities);
            for (int d = 0; d < densities; ++d) before[d] = environment->density_vector(voxel)[d];
            predicted = before;
            group.net.assign(densities, 0.);
            group.requested.resize(densities);
            group.applied.resize(densities);
        }
        for (size_t i = group.begin; i < group.end; ++i) {
            Cell* cell = cells[i];
            if (audit) {
                const auto& secretion = cell->phenotype.secretion;
                // Secretion::advance synchronizes the volume on first association.
                const double volume = cell->secretion_rates != &secretion.secretion_rates
                    ? cell->phenotype.volume.total : static_cast<BioFVM::Basic_Agent*>(cell)->get_total_volume();
                const double factor = dt * volume / voxel_volume;
                for (int d = 0; d < densities; ++d) {
                    const double net = dt * secretion.net_export_rates[d];
                    group.net[d] += net;
                    predicted[d] = (predicted[d] + factor * secretion.secretion_rates[d]
                                    * secretion.saturation_densities[d])
                        / (1. + factor * (secretion.secretion_rates[d] + secretion.uptake_rates[d]))
                        + net / voxel_volume;
                }
            }
            cell->phenotype.secretion.advance(cell, cell->phenotype, dt);
        }
        if (audit) {
            for (int d = 0; d < densities; ++d) {
                group.requested[d] = (predicted[d] - before[d]) * voxel_volume;
                group.applied[d] = (environment->density_vector(voxel)[d] - before[d]) * voxel_volume;
            }
        }
    }
    if (audit) {
        // Deterministic reduction; no atomics or MPI reductions in the hot loop.
        std::vector<double> net, requested, applied;
        for (const auto& group : groups) {
            net.resize(std::max(net.size(), group.net.size()), 0.);
            requested.resize(net.size(), 0.); applied.resize(net.size(), 0.);
            for (size_t d = 0; d < group.net.size(); ++d) {
                net[d] += group.net[d]; requested[d] += group.requested[d]; applied[d] += group.applied[d];
            }
        }
        budget.record(time, dt, net, requested, applied);
    }
}
} // namespace PhysiCell
#endif
