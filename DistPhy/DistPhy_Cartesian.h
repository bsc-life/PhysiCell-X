#ifndef __DistPhy_Cartesian_h__
#define __DistPhy_Cartesian_h__

#include<mpi.h>
#include "DistPhy_Environment.h"                    //When compiled separately, only forward declaration was creating problems, so possibly requires full class.
                                                    //See if we can possibly remove this
#include<vector>
#include<iostream>


namespace DistPhy
{
    namespace mpi
    {
        class mpi_Environment;                      //Forward declaration
        
        class mpi_Cartesian
        {
            public:
                MPI_Comm mpi_cart_comm;             //The new Cartesian communicator built from old communicator "init_comm"
                int mpi_dims[3];                    //Dimensions of the new communicator
                int mpi_coords[3];                  //MPI Cartesian coordinates of each MPI process
                int X_LEFT, X_RIGHT; 								//Cartesian left and right (as just 1-D x-decomposition) 
            public:
                void Build_Cartesian_Topology(mpi_Environment &world);          //Pass by reference to prevent creation of new object
                void Find_Cartesian_Coordinates(mpi_Environment &world);        //If new object created then possibly need to write copy constructor !
            		void Find_Left_Right_Neighbours(mpi_Environment &world);				//Finds left/right neighbour processes 
            		void Send_and_Receive_Cells
            				 (
            					int cells_to_left,  int size_left,  std::vector<char> & snd_buf_left, 
       								int cells_to_right, int size_right, std::vector<char> & snd_buf_right,
       								int & cells_from_left,  std::vector<char> &rcv_buf_left,
       								int & cells_from_right, std::vector<char> &rcv_buf_right,
       								std::vector<int> & snd_pos_left, std::vector<int> & snd_pos_right,
       								std::vector<int> & rcv_pos_left, std::vector<int> & rcv_pos_right,
       								mpi_Environment &world
       							); 
        };
    }
}

#endif
