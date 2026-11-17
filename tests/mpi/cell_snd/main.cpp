#include <mpi.h>
#include <iostream>
#include <vector>
#include <omp.h>
#include <filesystem>
#include <cmath>

#include "../../../core/PhysiCell.h"
#include "../../../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM;
using namespace PhysiCell;

using namespace DistPhy::mpi; 

void create_cells(std::vector<Cell>& cells, int n_cells){
    cells.clear();
    cells.resize(n_cells);
    for(int i = 0; i < n_cells; ++i){
        // Use a realistic cell-definition-backed fixture. initialize_random()
        // deliberately constructs arbitrary, internally inconsistent models and
        // types, which cannot be restored through the normal MPI receive path.
        cells[i].ID = i + 1000;
        cells[i].position = {0.25 * i, -0.5 * i, 0.75 * i};
        cells[i].velocity = {0.01 * i, 0.02 * i, -0.03 * i};
        cells[i].displacement = {0.04 * i, -0.05 * i, 0.06 * i};
        cells[i].state.simple_pressure = 0.125 * i;
        cells[i].phenotype.volume.total += i;
        cells[i].custom_data["serialization_marker"] = 1.5 + i;
    }
}

void snd_cells_parallel(std::vector<Cell>& cells){
    int n_threads = omp_get_max_threads();
    std::vector<MPI_Request> requests(2 * n_threads);
    long int total_bytes = 0;
    std::vector<std::vector<char>> buffer(n_threads);
    int base = cells.size()/n_threads;
    int remainder = cells.size()%n_threads;

    auto start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel reduction(+:total_bytes)
    {
        int tid = omp_get_thread_num();
        int start_cell = tid * base + std::min(tid, remainder);
        int end_cell = start_cell + base + (tid < remainder ? 1 : 0);
        int position = 0;
        int length = 0;
        std::vector<char> bf;
        for(int i = start_cell; i < end_cell; ++i){
            cells[i].pack(bf, length, position);
        }
        total_bytes += position;
        buffer[tid] = std::move(bf);
        int size = buffer[tid].size();
        MPI_Send(&size, 1, MPI_INT, 1, (tid*2), MPI_COMM_WORLD);
        //send data
        MPI_Send(buffer[tid].data(), size, MPI_CHAR, 1, (tid*2) + 1, MPI_COMM_WORLD);
    }
   
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "[Rank 0] Packed " << cells.size() << " cells  ("<< total_bytes/1024 <<"KB) in " << duration.count() << " us\n"; 
}

bool cell_matches_fixture(Cell& cell, int index) {
    const auto close = [](double first, double second) {
        return std::fabs(first - second) <= 1e-12;
    };
    return cell.ID == index + 1000 &&
           cell.type == cell_defaults.type &&
           cell.type_name == cell_defaults.name &&
           close(cell.position[0], 0.25 * index) &&
           close(cell.position[1], -0.5 * index) &&
           close(cell.position[2], 0.75 * index) &&
           close(cell.velocity[0], 0.01 * index) &&
           close(cell.velocity[1], 0.02 * index) &&
           close(cell.velocity[2], -0.03 * index) &&
           close(cell.displacement[0], 0.04 * index) &&
           close(cell.displacement[1], -0.05 * index) &&
           close(cell.displacement[2], 0.06 * index) &&
           close(cell.state.simple_pressure, 0.125 * index) &&
           close(cell.custom_data["serialization_marker"], 1.5 + index) &&
           close(cell.phenotype.volume.total,
                 cell_defaults.phenotype.volume.total + index);
}

int rcv_cells_parallel_v2(std::vector<Cell>& cells, int n_cells) {
    int n_threads = omp_get_max_threads();
    cells.resize(n_cells);

    int base = n_cells / n_threads;
    int remainder = n_cells % n_threads;

    auto start = std::chrono::high_resolution_clock::now();

    // 2. Unpack buffers in parallel
    int validation_failed = 0;
    #pragma omp parallel reduction(max:validation_failed)
    {
        int tid = omp_get_thread_num();
        int size;
        MPI_Recv(&size, 1, MPI_INT, 0, (tid * 2), MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::vector<char> bf(size);
        MPI_Recv(bf.data(), size, MPI_CHAR, 0, (tid * 2) + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int chunk_start = tid * base + std::min(tid, remainder);
        int chunk_end = chunk_start + base + (tid < remainder ? 1 : 0);

       
        int position = 0;
        int length = bf.size();

        for (int i = chunk_start; i < chunk_end; ++i) {
            Cell cell;
            cell.unpack(bf, length, position);
            cells[i] = std::move(cell);
#ifdef VALIDATE
            if (!cell_matches_fixture(cells[i], i)) {
                validation_failed = 1;
            }
#endif
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "[Rank 1] Unpacked " << cells.size() << " cells in " << duration.count() << " us\n";
    return validation_failed;
}

//This one is slower!
void rcv_cells_parallel(std::vector<Cell>& cells, int n_cells) {
    int n_threads = omp_get_max_threads();
    cells.resize(n_cells);

    std::vector<std::vector<char>> buffer(n_threads);
    std::vector<int> sizes(n_threads);

    int base = n_cells / n_threads;
    int remainder = n_cells % n_threads;

    auto start = std::chrono::high_resolution_clock::now();

    // 1. Receive all buffers serially per-thread (for safety)
    
    for (int tid = 0; tid < n_threads; ++tid) {
        int size;
        MPI_Recv(&size, 1, MPI_INT, 0, (tid * 2), MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        sizes[tid] = size;

        buffer[tid].resize(size);
        MPI_Recv(buffer[tid].data(), size, MPI_CHAR, 0, (tid * 2) + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // 2. Unpack buffers in parallel
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int chunk_start = tid * base + std::min(tid, remainder);
        int chunk_end = chunk_start + base + (tid < remainder ? 1 : 0);

        std::vector<char>& bf = buffer[tid];
        int position = 0;
        int length = bf.size();

        for (int i = chunk_start; i < chunk_end; ++i) {
            Cell cell;
            cell.unpack(bf, length, position);
            cells[i] = std::move(cell);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "[Rank 1] Unpacked " << cells.size() << " cells in " << duration.count() << " us\n";
}

//Serial snd
void pack_cells( std::vector<Cell>& cells){
    std::vector<char> buffer;
    int position = 0;
    int length = 0;
    for (int i = 0; i < cells.size(); ++i) {
        cells[i].pack(buffer, length, position);
        //std::cout << "  Packed cell " << i << " | pos: " << position << " | length: " << length << std::endl;
    }

    int snd = buffer.size();
    std::cout << "[Rank 0] Buffer size is " << buffer.size() << "\n";
        
    // Send buffer length first
    MPI_Send(&snd, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    // Then send the actual buffer
    MPI_Send(buffer.data(), buffer.size(), MPI_CHAR, 1, 1, MPI_COMM_WORLD);
}

//Serial rcv
void unpack_cells(std::vector<Cell>& cells) {
    // Receive buffer length
    std::vector<char> buffer;
    int recv_len = 0;
    MPI_Recv(&recv_len, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
    buffer.resize(recv_len);
    //std::cout << "[Rank "<< rank << "] Buffer size is " << buffer.size() << "\n";

    // Receive buffer data
    MPI_Recv(buffer.data(), recv_len, MPI_CHAR, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    int position = 0;
    int length = buffer.size();
    int cell_num = 0;
    while (position < buffer.size()) {
        Cell cell;
        cell.unpack(buffer, length, position);
        cells.push_back(cell);
        //std::cout << "  Unpacked cell " << cell_num << " | pos: " << position << " | length: " << length << std::endl;
        ++cell_num;
    }
}

namespace fs = std::filesystem;

bool create_directory(const std::string& folder_name) {
    if (!fs::exists(folder_name)) {
        return fs::create_directory(folder_name);
    }
    return true;  // already exists
}

std::vector<fs::path> get_files_in_directory(const std::string& folder_name) {
    std::vector<fs::path> files;
    for (const auto& entry : fs::directory_iterator(folder_name)) {
        if (fs::is_regular_file(entry.path())) {
            files.push_back(entry.path());
        }
    }
    return files;
}

bool files_are_equal(const fs::path& file1, const fs::path& file2) {
    std::string command = "diff -q \"" + file1.string() + "\" \"" + file2.string() + "\" > /dev/null";
    int result = system(command.c_str());
    return result == 0; // 0 means files are identical
}


int main(int argc, char** argv) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided < MPI_THREAD_MULTIPLE) {
        
        std::cerr << "Error: MPI does not support MPI_THREAD_MULTIPLE\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    initialize_default_cell_definition();
    cell_defaults.custom_data.add_variable("serialization_marker", "dimensionless", 1.5);
    build_cell_definitions_maps();

    int n_cells;
    if (argc != 2) {
        if (rank == 0) std::cout << "Error: expercted 1 argument: number of cells" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 2);
    } else {
        n_cells = stoi(argv[1]);
    }

    if (rank == 0) std::cout << "SND/RCV unit test will use " << n_cells << " cells" << std::endl;
    
    int buffer_len = 0;
    int position = 0;

    int validation_failed = 0;

    if (rank == 0) {
        
        vector<Cell> cells;
        auto start = std::chrono::high_resolution_clock::now();
        create_cells(cells, n_cells);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "[Rank 0] Created " << cells.size() << " cells in " << duration.count() << " us" << std::endl;

        start = std::chrono::high_resolution_clock::now();
        //pack_cells(cells);
        snd_cells_parallel(cells);
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "[Rank 0] Packed " << cells.size() << " cells in " << duration.count() << " us \n";
        
        
#ifdef VALIDATE
        MPI_Recv(&validation_failed, 1, MPI_INT, 1, 1000, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
#endif
    }
    else if (rank == 1) {
       
        vector<Cell> cells;
        auto start = std::chrono::high_resolution_clock::now();
        //unpack_cells(cells);
        validation_failed = rcv_cells_parallel_v2(cells, n_cells);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "[Rank 1] Unpacked " << cells.size() << " cells in " << duration.count() << " us" << std::endl;
        std::cout << "Test completed!" << std::endl;

#ifdef VALIDATE
        MPI_Send(&validation_failed, 1, MPI_INT, 0, 1000, MPI_COMM_WORLD);
#endif
    }

    MPI_Barrier(MPI_COMM_WORLD);

#ifdef VALIDATE
    MPI_Bcast(&validation_failed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0 && validation_failed == 0) {
        std::cout << "PASS: every serialized cell preserved its identity, kinematics, phenotype, and custom data"
                  << std::endl;
    }
#endif
    MPI_Finalize();
    return validation_failed;
}
