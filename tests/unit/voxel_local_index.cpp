// Exercise the production MPI-local indexing function without starting MPI.
#include "../../BioFVM/BioFVM_mesh.h"
#include <iostream>

int main()
{
    unsigned checks = 0, failures = 0;
    auto check = [&](BioFVM::Cartesian_Mesh& mesh,
                     DistPhy::mpi::mpi_Environment& world,
                     DistPhy::mpi::mpi_Cartesian& topology,
                     std::vector<double> position, int expected) {
        const int actual = mesh.nearest_voxel_local_index(position, world, topology);
        ++checks;
        if (actual != expected) {
            if (failures < 12)
                std::cerr << "rank=" << world.rank << " position=("
                          << position[0] << ',' << position[1] << ',' << position[2]
                          << ") expected=" << expected << " actual=" << actual << '\n';
            ++failures;
        }
    };

    // Actual run geometry, uneven x partitioning, and single-voxel local axes.
    for (int ranks : {1, 2, 3}) {
        for (bool singleton : {false, true}) {
            for (int rank = 0; rank < ranks; ++rank) {
                DistPhy::mpi::mpi_Environment world{};
                world.rank = rank;
                world.size = ranks;
                DistPhy::mpi::mpi_Cartesian topology{};
                topology.mpi_dims[0] = topology.mpi_dims[2] = 1;
                topology.mpi_dims[1] = ranks;
                topology.mpi_coords[1] = rank;
                BioFVM::Cartesian_Mesh mesh;
                if (singleton)
                    mesh.resize(0., 2.*ranks, 0., 2., 0., 2.,
                                2., 2., 2., world, topology);
                else
                    mesh.resize(-250., 250., -250., 250., -4., 4.,
                                2., 2., 2., world, topology);

                // Every stored voxel centre must map back to its storage index.
                for (std::size_t i = 0; i < mesh.voxels.size(); ++i)
                    check(mesh, world, topology, mesh.voxels[i].center, i);

                const int sizes[] = {static_cast<int>(mesh.x_coordinates.size()),
                                     static_cast<int>(mesh.y_coordinates.size()),
                                     static_cast<int>(mesh.z_coordinates.size())};
                const int strides[] = {sizes[1]*sizes[2], sizes[2], 1};
                const double spacing[] = {mesh.dx, mesh.dy, mesh.dz};
                const auto first = mesh.voxels.front().center;
                for (int axis = 0; axis < 3; ++axis) {
                    auto position = first;
                    const double lo = mesh.local_bounding_box[axis];
                    const double hi = mesh.local_bounding_box[axis+3];
                    // Small and large undershoots must clamp to the first voxel.
                    for (double coordinate : {lo-1e-8, lo-10*spacing[axis], lo}) {
                        position[axis] = coordinate;
                        check(mesh, world, topology, position, 0);
                    }
                    for (double coordinate : {hi-1e-8, hi, hi+10*spacing[axis]}) {
                        position[axis] = coordinate;
                        check(mesh, world, topology, position,
                              (sizes[axis]-1)*strides[axis]);
                    }
                    if (sizes[axis] > 1) {
                        position[axis] = lo+spacing[axis]-1e-8;
                        check(mesh, world, topology, position, 0);
                        position[axis] = lo+spacing[axis];
                        check(mesh, world, topology, position, strides[axis]);
                    }
                }
                auto corner = first;
                for (int axis = 0; axis < 3; ++axis)
                    corner[axis] = mesh.local_bounding_box[axis]-1e-8;
                check(mesh, world, topology, corner, 0);
                for (int axis = 0; axis < 3; ++axis)
                    corner[axis] = mesh.local_bounding_box[axis+3]+1e-8;
                check(mesh, world, topology, corner, mesh.voxels.size()-1);
            }
        }
    }
    std::cout << checks << " checks, " << failures << " failures\n";
    return failures ? 1 : 0;
}
