// Exercise production Cartesian_Mesh construction without launching MPI.
// The constructor only needs the rank/topology values; it communicates nothing.
#include "../../BioFVM/BioFVM_mesh.h"
#include <algorithm>
#include <iostream>

int main(int argc, char** argv)
{
    const double z_extent = argc > 1 ? std::stod(argv[1]) : 8.0;
    BioFVM::Cartesian_Mesh meshes[2];
    for (int rank = 0; rank < 2; ++rank)
    {
        DistPhy::mpi::mpi_Environment world{};
        world.rank = rank;
        world.size = 2;
        DistPhy::mpi::mpi_Cartesian topo{};
        topo.mpi_dims[0] = 1;
        topo.mpi_dims[1] = 2;
        topo.mpi_dims[2] = 1;
        topo.mpi_coords[1] = rank;
        meshes[rank].resize(-250., 250., -250., 250., -z_extent/2, z_extent/2,
                            5., 5., 5., world, topo);
    }
    int absent = 0, total = 0;
    for (int rank = 0; rank < 2; ++rank)
    {
        auto& mine = meshes[rank];
        auto& other = meshes[1-rank];
        int nx = mine.x_coordinates.size();
        int ny = mine.y_coordinates.size();
        int nz = mine.z_coordinates.size();
        const auto& lists = rank == 0 ? mine.moore_connected_voxel_global_indices_right
                                      : mine.moore_connected_voxel_global_indices_left;
        for (int j = 0; j < ny; ++j)
        for (int k = 0; k < nz; ++k)
        {
            const int local = mine.voxel_index(rank == 0 ? nx-1 : 0, j, k);
            const auto& list = lists[j*nz+k];
            for (int jj = std::max(0, j-1); jj <= std::min(ny-1, j+1); ++jj)
            for (int kk = std::max(0, k-1); kk <= std::min(nz-1, k+1); ++kk)
            {
                const int expected = other.voxels[other.voxel_index(rank == 0 ? 0 : other.x_coordinates.size()-1, jj, kk)].global_mesh_index;
                ++total;
                if (std::find(list.begin(), list.end(), expected) == list.end()) ++absent;
            }
            if (j == ny/2 && k == 0)
            {
                const int expected = other.voxels[other.voxel_index(rank == 0 ? 0 : other.x_coordinates.size()-1, j, k)].global_mesh_index;
                std::cout << "rank=" << rank << " actual_z_layers=" << nz
                          << " voxel=" << mine.voxels[local].global_mesh_index
                          << " expected_face_neighbor=" << expected << " candidates=";
                for (int id : list) std::cout << id << ',';
                std::cout << '\n';
            }
        }
    }
    std::cout << "z_extent=" << z_extent << " missing_directed_neighbors=" << absent
              << " expected_directed_neighbors=" << total << '\n';
    return absent == 0 ? 0 : 1;
}
