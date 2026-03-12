

#ifndef __DistPhy_Utils_h__
#define __DistPhy_Utils_h__

#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "./DistPhy_Environment.h"
#include "./DistPhy_Cartesian.h"

using namespace std; 

namespace DistPhy
{
    namespace mpi
    {
        class mpi_Environment;  //Forward declaration
        class mpi_Cartesian;    //Forward declaration
        
        class mpi_CellPositions
        {
        public:
            
            std::vector<std::vector<double>> cell_coords_all_procs; 
            std::vector<std::vector<int>>    cell_IDs_all_procs; 
            // std::vector<std::vector<int>> cell_types_all_procs;
            std::vector<int>                 no_of_IDs_all_procs; 
            
        public:
            void positions_to_rank_list(std::vector<std::vector<double>> &generated_list, 
                                        int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int dx, int dy, int dz, 
                                        mpi_Environment &world, mpi_Cartesian &cart_topo, int strt_cell_ID ); 
            
        };
        
        /*---------------------------------------------------------------------------------------------------------------*/
        /* Class contains data-structures needed to store no. of cell IDs, cell IDs and cell coordinates at each process */
        /* Each process will create an object of this class, even the root process as well                               */
        /*---------------------------------------------------------------------------------------------------------------*/
        class mpi_MyCells 
        {
        public:
            std::vector<double> my_cell_coords;
            std::vector<int>    my_cell_IDs;
            // std::vector<int> my_cell_types
            int                 my_no_of_cell_IDs; 

        };
        
        void distribute_cell_positions(mpi_CellPositions &cp, mpi_MyCells &mc, mpi_Environment &world, mpi_Cartesian &cart_topo); 
        int request_IDs_from_root(int no_of_requested_IDs, int *range_ID, int strt_cell_ID, mpi_Environment &world, mpi_Cartesian &cart_topo); 
    }
}

#endif
