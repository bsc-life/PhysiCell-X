#include<mpi.h>

#include "DistPhy_Cartesian.h"     //See if we can possibly remove this             

namespace DistPhy
{
    namespace mpi
    {
        void mpi_Cartesian::Build_Cartesian_Topology(mpi_Environment &world)
        {
            int is_periodic[3]; 
            int mpi_reorder;
            int dimensions = 3; 
            
            mpi_dims[0] = 1;                    //mpi_dims[] is a member data of mpi_Environment
            mpi_dims[1] = world.size;           //To create X-decomposition, we need to specify MPI processes in the Y-direction
            mpi_dims[2] = 1;                    //PhysiCell/BioFVM data direction is left to right, hence Y-direction MPI produce cuts in X-direction of PhysiCell/BioFVM data i.e. X-decomposition
            
            is_periodic[0] = 0;                 //No direction is periodic
            is_periodic[1] = 0;
            is_periodic[2] = 0;
            
            mpi_reorder = 0;                    //No need to reorder MPI processes in new communicator
            MPI_Cart_create(world.init_comm, dimensions, mpi_dims, is_periodic, mpi_reorder, &mpi_cart_comm); 
        }
        
        void mpi_Cartesian::Find_Cartesian_Coordinates(mpi_Environment &world)
        {
            int dimensions = 3; 
            
            MPI_Cart_coords(mpi_cart_comm, world.rank, dimensions, mpi_coords);         //The class public data mpi_coords[] is filled with coordinates in 3D space
        }
        
        void mpi_Cartesian::Find_Left_Right_Neighbours(mpi_Environment &world)
        {
        	int dimension 	 = 1; //MPI Y direction is X direction in PhysiCell
        	int displacement = 1; //Need to find immediate neighbours only
        	
        	MPI_Cart_shift(mpi_cart_comm, dimension, displacement, &X_LEFT, &X_RIGHT);
        }
        
       void mpi_Cartesian::Send_and_Receive_Cells
       														(
       														 int cells_to_left,  int size_left,  std::vector<char> & snd_buf_left, 
       														 int cells_to_right, int size_right, std::vector<char> & snd_buf_right,
       														 int & cells_from_left,  std::vector<char> & rcv_buf_left,
       														 int & cells_from_right, std::vector<char> & rcv_buf_right,
       														 mpi_Environment &world
       														)
				{
						/*==========================================================================================*/
						/* cells_to_left    - Number of cells going to left neighbour (Pass by value)								*/
						/* size_left  	    - Total size in bytes of cells going to left neighbour (Pass by value)	*/
						/* snd_buf_left     - Actual buffer which contains packed data	(Pass by Ref)								*/
						/* cells_to_right   - Number of cells going to right neighbour (Pass by value)							*/
						/* size_right  	    - Total size in bytes of going to right neighbour (Pass by value)				*/
						/* snd_buf_right 	  - Actual buffer which contains packed data	(Pass by Ref)								*/
						/* cells_from_left  - Cells coming from left neighbour (Pass by Ref)												*/
						/* rcv_buf_left 	  - Buffer where data from left neighbour is received (Pass by Ref)				*/
						/* cells_from_right - Cells coming from right neighbour (Pass by Ref)												*/
						/* rcv_buf_right 	  - Buffer where data from right neighbour is received (Pass by Ref)			*/
						/*==========================================================================================*/
								
				  MPI_Request snd_req[2], rcv_req[2];
								
					int snd_cells_and_size_to_left[2]; 
					int rcv_cells_and_size_from_right[2];
					
					int snd_cells_and_size_to_right[2];
					int rcv_cells_and_size_from_left[2];
					
					/*-------------------------------------------------------------------*/			
					/* Will send both no_of_cells and size together in one communication */
				  /*-------------------------------------------------------------------*/			
 
								
					snd_cells_and_size_to_left[0]  = cells_to_left;
					snd_cells_and_size_to_left[1]  = size_left;
					 
					
					snd_cells_and_size_to_right[0] = cells_to_right;
					snd_cells_and_size_to_right[1] = size_right;
					
					if(cells_to_left > 0)
						std::cout<<"RANK="<<world.rank<<" IS SENDING "<<cells_to_left<<" CELLS TO LEFT WITH TOTAL SIZE "<<size_left<<" BYTES"<<std::endl;
					
					if(cells_to_right > 0)
						std::cout<<"RANK="<<world.rank<<" IS SENDING "<<cells_to_right<<" CELLS TO RIGHT WITH TOTAL SIZE "<<size_right<<" BYTES"<<std::endl;
					
					
					/*------------------------------------------------------------------------------------------*/
					/* Will need to make rcv_cells_and size_from_right[0]=0 because if this is not initialized 	*/
					/* then cells_from_right = rcv_cells_and size_from_right[0] gets a garbage value on Rank 1	*/
					/* Similarly cells_from_left = rcv_cells_and_size_from_left[0] can get a garbage value on 	*/
					/* Rank 0 (I am assuming only 2 process scenario but this holds for any 'n' rank scenario)	*/
					/*------------------------------------------------------------------------------------------*/
					
					rcv_cells_and_size_from_right[0] = 0;
					rcv_cells_and_size_from_left[0]	 = 0; 
					
					
					/*------------------------------------------------------------------------------------------*/
					/* This communication is necessary for ALL processes even if no cell crosses the sub-domain */							
					/*------------------------------------------------------------------------------------------*/
					
					MPI_Irecv(rcv_cells_and_size_from_right, 2, MPI_INT, X_RIGHT, 1111, MPI_COMM_WORLD, &rcv_req[0]); 
					MPI_Irecv(rcv_cells_and_size_from_left,  2, MPI_INT, X_LEFT,  2222, MPI_COMM_WORLD, &rcv_req[1]);
								 
					MPI_Isend(snd_cells_and_size_to_left,  2, MPI_INT, X_LEFT,  1111, MPI_COMM_WORLD, &snd_req[0]);
					MPI_Isend(snd_cells_and_size_to_right, 2, MPI_INT, X_RIGHT, 2222, MPI_COMM_WORLD, &snd_req[1]);		

					MPI_Waitall(2, snd_req, MPI_STATUSES_IGNORE); 
					MPI_Waitall(2, rcv_req, MPI_STATUSES_IGNORE);
					
					/*================================================================================*/
					/* These 2 statements were missing i.e. the cells_from_left/right data members 		*/
					/* were not being set and hence the unpack() function always WAS seeing the value */
					/* of these 2 data members as 0. 																									*/
					/*================================================================================*/
					
					cells_from_right = rcv_cells_and_size_from_right[0];
					cells_from_left  = rcv_cells_and_size_from_left[0];
					
					/*------------------------------------------------------------------------------*/
					/* First handle data flowing to left i.e. sending TO left, receiving FROM right */
					/* process_0 <---- process_1 <-----process_2 <----......<----process_n 					*/
					/*------------------------------------------------------------------------------*/
					
					if(cells_from_right > 0)
						rcv_buf_right.resize(rcv_cells_and_size_from_right[1]);
					
					/*------------------------------------------------------------------------------*/
					/* Next handle data flowing to right i.e. sending TO right, receiving FROM left */
					/* process_0---> process_1--->process_2--->...--->process_n 										*/
					/*------------------------------------------------------------------------------*/
			
					if(cells_from_left > 0)
						rcv_buf_left.resize(rcv_cells_and_size_from_left[1]); 
						
					/*-------------------------------------------------------------------------*/
					/* Now the actual data i.e. packed buffers will be sent 									 */
					/* process_0 <---- process_1 <-----process_2 <----......<----process_n		 */
					/* Send data towards left IFF snd_cells_and_size_to_left[0] > 0 					 */
					/* Receive data coming from right IFF rcv_cells_and_size_from_right[0] > 0 */
					/*-------------------------------------------------------------------------*/
					
					if(cells_from_right > 0)
						MPI_Irecv(&rcv_buf_right[0], rcv_cells_and_size_from_right[1], MPI_PACKED, X_RIGHT, 3333, MPI_COMM_WORLD, &rcv_req[0]);
					
					if(cells_to_left > 0)
						MPI_Isend(&snd_buf_left[0], size_left, MPI_PACKED, X_LEFT, 3333, MPI_COMM_WORLD, &snd_req[0]);
					
					/*----------------------------------------------------------------*/
					/* Now the actual data i.e. packed buffers will be sent 					*/
					/* process_0---> process_1--->process_2--->...--->process_n				*/
					/* Send data towards right IFF snd_cells_and_size_to_right[0] > 0 */ 
					/* Receive data from left IFF rcv_cells_and_size_from_left[0] > 0 */
					/*----------------------------------------------------------------*/
					
					if(cells_from_left > 0)
						MPI_Irecv(&rcv_buf_left[0], rcv_cells_and_size_from_left[1], MPI_PACKED, X_LEFT, 4444, MPI_COMM_WORLD, &rcv_req[1]);
					
					if(cells_to_right > 0)
						MPI_Isend(&snd_buf_right[0], size_right, MPI_PACKED, X_RIGHT, 4444, MPI_COMM_WORLD, &snd_req[1]);
						
				  /*---------------------------------------------------------------------*/
				  /* We need to wait ONLY IF some communication has taken place, else NO */ 
					/*---------------------------------------------------------------------*/
				  
				  if(cells_to_left  > 0)
				  	MPI_Wait(&snd_req[0], MPI_STATUS_IGNORE);
				  	
				  if(cells_to_right > 0)
				  	MPI_Wait(&snd_req[1], MPI_STATUS_IGNORE);
				  	
				  if(cells_from_right > 0)
				  	MPI_Wait(&rcv_req[0], MPI_STATUS_IGNORE);
				  	
				  if(cells_from_left > 0)
				  	MPI_Wait(&rcv_req[1], MPI_STATUS_IGNORE);				  
								
				}
    }
}
