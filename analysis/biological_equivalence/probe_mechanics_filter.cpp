// Isolate the production neighbor-voxel filter with the dFBA cell radius.
#include "../../core/PhysiCell_cell.h"
#include <cmath>
#include <iostream>

int main()
{
    PhysiCell::Cell cell(900000);
    cell.position = {-1., 0., 0.};
    cell.phenotype.geometry.radius = 0.7777484771851859;
    cell.phenotype.mechanics.relative_maximum_adhesion_distance = 1.;
    DistPhy::mpi::mpi_Environment world{};
    DistPhy::mpi::mpi_Cartesian topo{};
    const std::vector<double> local_center = {-2.5, 2.5, -1.5};
    const std::vector<double> remote_center = {2.5, 2.5, -1.5};
    const bool empty_cache = PhysiCell::is_neighbor_voxel(
        &cell, local_center, remote_center, 0., world, topo);
    const bool populated_cache = PhysiCell::is_neighbor_voxel(
        &cell, local_center, remote_center, cell.phenotype.geometry.radius, world, topo);
    // A second equal-radius cell at (+0.4, 0, 0) physically overlaps this one.
    const double overlap = 1. - 1.4/(2*cell.phenotype.geometry.radius);
    std::cout << "radius=" << cell.phenotype.geometry.radius
              << " pair_overlap=" << overlap
              << " accepted_with_zero_cache=" << empty_cache
              << " accepted_with_correct_cache=" << populated_cache << '\n';
    return !empty_cache && populated_cache && overlap > 0 ? 0 : 1;
}
