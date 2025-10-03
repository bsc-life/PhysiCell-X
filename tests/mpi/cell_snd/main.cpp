#include <mpi.h>
#include <iostream>
#include <vector>
#include <omp.h>
#include <filesystem>

#include "../../../core/PhysiCell.h"
#include "../../../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM;
using namespace PhysiCell;

using namespace DistPhy::mpi; 

void create_cells(std::vector<Cell>& cells, int n_cells){
    cells.clear();
    cells.resize(n_cells);
    for(int i = 0; i < n_cells; ++i){
        cells[i].initialize_random();
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

void rcv_cells_parallel_v2(std::vector<Cell>& cells, int n_cells) {
    int n_threads = omp_get_max_threads();
    cells.resize(n_cells);

    int base = n_cells / n_threads;
    int remainder = n_cells % n_threads;

    auto start = std::chrono::high_resolution_clock::now();

    // 2. Unpack buffers in parallel
    #pragma omp parallel
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
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "[Rank 1] Unpacked " << cells.size() << " cells in " << duration.count() << " us\n";
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
            create_directory("./tmp_cells/");
            for (int i = 0; i < n_cells; i++){
                std::string path = "./tmp_cells/snd_" + to_string(i) + ".cell";
                std::ofstream outFile(path);
                cells[i].print_parameters(outFile); 
            }
        #endif
    }
    else if (rank == 1) {
       
        vector<Cell> cells;
        auto start = std::chrono::high_resolution_clock::now();
        //unpack_cells(cells);
        rcv_cells_parallel_v2(cells, n_cells);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "[Rank 1] Unpacked " << cells.size() << " cells in " << duration.count() << " us" << std::endl;
        std::cout << "Test completed!" << std::endl;

        #ifdef VALIDATE
            create_directory("./tmp_cells/");
            for (int i = 0; i < n_cells; i++){
                std::string path = "./tmp_cells/rcv_" + to_string(i) + ".cell";
                std::ofstream outFile(path);
                cells[i].print_parameters(outFile); 
            }
        #endif
    }

    MPI_Barrier(MPI_COMM_WORLD);

    #ifdef VALIDATE
    if (rank == 0) {
        int cells_different = 0;
        for (int i = 0; i < n_cells; ++i){
            std::string cell_snd = "./tmp_cells/snd_" + to_string(i) + ".cell";
            std::string cell_rcv = "./tmp_cells/rcv_" + to_string(i) + ".cell";

            bool equal = files_are_equal(cell_snd, cell_rcv);
            if (equal){
                std::remove(cell_snd.c_str());
                std::remove(cell_rcv.c_str());
            } else {
                ++cells_different;
            }
        }
        if (cells_different == 0) std::cout << "Validation succesful: no cells have been corrupted" << std::endl;
        else std::cout << "Validation error! Cells that have been corrupted: "<< cells_different << std::endl;
    }
    #endif
    MPI_Finalize();
    return 0;
}
