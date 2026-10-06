#include "core/PhysiCell_cell_sources.h"
#include "core/PhysiCell_utilities.h"
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace BioFVM;
using namespace PhysiCell;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

int main() {
    omp_set_dynamic(0);
    omp_set_num_threads(8);
    SeedRandom(50);
    Microenvironment& m = microenvironment;
    m.add_density("test_product", "dimensionless", 0., 0.);
    m.resize_space(-4., 4., -4., 4., -2., 2., 2., 2., 2.);
#ifdef SOURCE_TEST_MPI
    m.mesh.local_bounding_box = m.mesh.bounding_box;
#endif
    default_microenvironment_options.track_internalized_substrates_in_each_agent = true;
    std::vector<Cell*> cells;
    for (int id = 0; id < 64; ++id) {
#ifdef SOURCE_TEST_MPI
        Cell* cell = new Cell(id);
#else
        Cell* cell = new Cell();
        cell->ID = id;
#endif
        cell->register_microenvironment(&m);
        cell->position = {-3. + 2.*(id % 4), -1., 1.};
        cell->update_voxel_index();
        cell->set_total_volume(1.);
        cell->phenotype.secretion.sync_to_microenvironment(&m);
        auto& s = cell->phenotype.secretion;
        // One purely additive field plus one field with implicit secretion/uptake.
        s.secretion_rates[0] = 0.; s.uptake_rates[0] = 0.;
        s.net_export_rates[0] = 8. + id;
        s.secretion_rates[1] = .1 * (1 + id % 3);
        s.uptake_rates[1] = .2;
        s.saturation_densities[1] = 3.;
        s.net_export_rates[1] = .5;
        cells.push_back(cell);
    }
    auto reset = [&]() {
        for (int n = 0; n < m.number_of_voxels(); ++n)
            for (int d = 0; d < m.number_of_densities(); ++d) m.density_vector(n)[d] = 1.;
        for (Cell* cell : cells) cell->internalized_substrates->assign(2, 0.);
    };
    auto capture = [&]() {
        std::vector<double> values;
        for (int n = 0; n < m.number_of_voxels(); ++n)
            for (int d = 0; d < m.number_of_densities(); ++d) values.push_back(m.density_vector(n)[d]);
        for (Cell* cell : cells)
            values.insert(values.end(), cell->internalized_substrates->begin(), cell->internalized_substrates->end());
        return values;
    };
    // Associate phenotype rates with the Basic_Agent before deliberately poisoning
    // its timestep cache, as happens after an intracellular update.
    reset();
    for (Cell* cell : cells) cell->phenotype.secretion.advance(cell, cell->phenotype, .001);
    reset();
    for (Cell* cell : cells) cell->set_internal_uptake_constants(.01);
    for (Cell* cell : cells) cell->phenotype.secretion.advance(cell, cell->phenotype, .001);
    const auto expected = capture();
    double actual = 0., requested = 0., internal = 0.;
    for (int n = 0; n < m.number_of_voxels(); ++n)
        actual += (m.density_vector(n)[0] - 1.) * m.voxels(n).volume;
    for (Cell* cell : cells) {
        requested += .001 * cell->phenotype.secretion.net_export_rates[0];
        internal += (*cell->internalized_substrates)[0];
    }
    require(std::abs(actual - requested) < 1e-12, "Wrong secretion timestep or exported amount");
    require(std::abs(actual + internal) < 1e-12, "Internal/external substrate accounting mismatch");
    for (int threads : {1, 2, 8}) {
        omp_set_num_threads(threads);
        auto shuffled = cells;
        std::reverse(shuffled.begin(), shuffled.end());
        for (int repeat = 0; repeat < 20; ++repeat) {
            reset();
            for (Cell* cell : cells) cell->set_internal_uptake_constants(.01);
            advance_cell_sources(shuffled, repeat * .001, .001);
            require(capture() == expected, "Source update depends on threads or input cell ordering");
        }
    }
    // Changed rates without a volume change must also invalidate any old increment.
    Cell* cell = cells[0];
    cell->set_internal_uptake_constants(.01);
    cell->phenotype.secretion.net_export_rates[0] = 80.;
    reset();
    cell->phenotype.secretion.advance(cell, cell->phenotype, .002);
    const int voxel = cell->get_current_voxel_index();
    require(std::abs((m.density_vector(voxel)[0] - 1.) * m.voxels(voxel).volume - .16) < 1e-12,
            "Source cache ignored new rates or timestep");
    std::cout << "PASS: 64 cells / 4 shared voxels, 1/2/8 threads, reversed input order, "
                 "cached timestep, changed rates, implicit uptake and internalized accounting\n";
}
