#include "DistPhy_Utils.h"

using namespace std;
using namespace DistPhy::mpi; 

namespace DistPhy
{
    namespace mpi
    {
        void mpi_CellPositions::positions_to_rank_list(std::vector<std::vector<double>> &generated_list, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int dx, int dy, int dz, mpi_Environment &world, mpi_Cartesian &cart_topo, int strt_cell_ID )
        {
            int dims[3];
            int proc_x_coord, proc_y_coord, proc_z_coord; 
            int proc_rank;
            double temp_point[3]; 
            
            cell_coords_all_procs.resize(world.size);               //List of coordinates of each process, each element maintains a list 3 * coordinates for that rank
            cell_IDs_all_procs.resize(world.size);                  //IDs corresponding to each position generated for each rank. 
            no_of_IDs_all_procs.resize(world.size,0);               //No. of positions for each rank, 0 if no positions for a particular rank. 
            
            
            dims[0] = cart_topo.mpi_dims[0]; 
            dims[1] = cart_topo.mpi_dims[1];
            dims[2] = cart_topo.mpi_dims[2];
            
            int global_x_voxels = (xmax-xmin)/dx; 
            int global_y_voxels = (ymax-ymin)/dy;
            int global_z_voxels = (zmax-zmin)/dz;
            
            int local_x_voxels = global_x_voxels/dims[1]; 
            int local_y_voxels = global_y_voxels/dims[0];
            int local_z_voxels = global_z_voxels/dims[2]; 
            
            for(size_t i=0; i < generated_list.size(); i++)
            {
                /*-----------------------------------------------------------*/
                /* temp_points[] are only used to calculate the process rank */
                /* they should NOT be put into cell_coords_all_procs[rank]   */
                /*-----------------------------------------------------------*/
                temp_point[0] = generated_list[i][0] - xmin; 
                temp_point[1] = generated_list[i][1] - ymin; 
                temp_point[2] = generated_list[i][2] - zmin; 
                
                proc_y_coord = temp_point[0]/(local_x_voxels * dx);                                     //Can use floor() here but no need
                proc_x_coord = (dims[0] - 1) - floor(temp_point[1]/(local_y_voxels * dy));              //Automatic promotion to double was creating problems
                proc_z_coord = temp_point[2]/(local_z_voxels * dz);                                     //Can use floor() here but no need 
                proc_rank = proc_x_coord * dims[1] * dims[2] + proc_y_coord * dims[2] + proc_z_coord;
                
                //std::cout<<"CELL ID="<<strt_cell_ID<<" ("<<generated_list[i][0]<<","<<generated_list[i][1]<<","<<generated_list[i][2]<<")"<<" Rank="<<proc_rank<<std::endl;
                
                cell_coords_all_procs[proc_rank].push_back(generated_list[i][0]);
                cell_coords_all_procs[proc_rank].push_back(generated_list[i][1]);
                cell_coords_all_procs[proc_rank].push_back(generated_list[i][2]);
                
                cell_IDs_all_procs[proc_rank].push_back(strt_cell_ID);
                strt_cell_ID++; 
                
                no_of_IDs_all_procs[proc_rank]++; 
                
                
            }
        }
        
        void distribute_cell_positions(mpi_CellPositions &cp, mpi_MyCells &mc, mpi_Environment &world, mpi_Cartesian &cart_topo)
        {
           MPI_Request snd_req[world.size]; 
           MPI_Request rcv_req; 
           
           /*------------------------------------------------------------------------*/
           /* First Scatter the number of cell IDs to each process (even if it is 0) */
           /*------------------------------------------------------------------------*/
           
           int send_count = 1; 
           int recv_count = 1;
           int root = 0; 
           
           MPI_Scatter(&cp.no_of_IDs_all_procs[0], send_count, MPI_INT, &mc.my_no_of_cell_IDs, recv_count, MPI_INT, root, cart_topo.mpi_cart_comm);
           std::cout<<"MPI Rank = "<<world.rank<<" No of cells = "<<mc.my_no_of_cell_IDs<<std::endl; 
           
           /*--------------------------------------------------------------------------------------------------------------*/
           /* Now resize the vectors inside mpi_MyCells on each process according to the no. of IDs that they will receive */
           /* (at the root process as well as IDs, coords will be copied to these arrays for homogeneity                   */
           /*--------------------------------------------------------------------------------------------------------------*/
           
           mc.my_cell_IDs.resize(mc.my_no_of_cell_IDs);
           mc.my_cell_coords.resize(3 * mc.my_no_of_cell_IDs);
           
           /*-----------------------------------------------------------------------------------------------------*/
           /* Now send the IDs from the root process (rank 0) to other processes (rank > 0) and then copy the IDs */
           /* on the root process from cp.cell_IDs_all_procs[0] to mc.my_cell_IDs                                 */
           /*-----------------------------------------------------------------------------------------------------*/
           
           if(world.rank == 0)
           {
               for(int i=1; i<=world.size-1; i++)
               {
                   if(cp.no_of_IDs_all_procs[i] > 0)
                   {
                       MPI_Isend(&cp.cell_IDs_all_procs[i][0], cp.no_of_IDs_all_procs[i], MPI_INT, i, 1111, cart_topo.mpi_cart_comm, &snd_req[i]); 
                       MPI_Wait(&snd_req[i], MPI_STATUS_IGNORE);
                   }
               }
               /*--------------------------------------------------------------------------------------*/
               /* At the root copy the IDs into the correct array for homogeneity with other processes */
               /*--------------------------------------------------------------------------------------*/
               
               for(int i=0; i<=cp.no_of_IDs_all_procs[0]-1; i++)
                   mc.my_cell_IDs[i] = cp.cell_IDs_all_procs[0][i]; 
           }
           else
           {
                if(mc.my_no_of_cell_IDs > 0)
                {
                    MPI_Irecv(&mc.my_cell_IDs[0], mc.my_no_of_cell_IDs, MPI_INT, 0, 1111, cart_topo.mpi_cart_comm, &rcv_req);
                    MPI_Wait(&rcv_req, MPI_STATUS_IGNORE);
                }
           }
           
           /*---------------------------------------------------------------------------------------------------------------------*/
           /* Now send the coordinates from the root process (rank 0) to other processes (rank > 1) and then copy the coordinates */
           /* on the root process to the appropriate array.                                                                       */
           /*---------------------------------------------------------------------------------------------------------------------*/ 
           
           if(world.rank == 0)
           {
               for(int i=1; i<=world.size-1; i++)
               {
                   if(cp.no_of_IDs_all_procs[i] > 0)
                   {
                       MPI_Isend(&cp.cell_coords_all_procs[i][0], 3 * cp.no_of_IDs_all_procs[i], MPI_DOUBLE, i, 2222, cart_topo.mpi_cart_comm, &snd_req[i]); 
                       MPI_Wait(&snd_req[i], MPI_STATUS_IGNORE);
                   }
               }
               /*--------------------------------------------------------------------------------------*/
               /* At the root copy the coords into the correct array for homogeneity with other processes */
               /*--------------------------------------------------------------------------------------*/
               
               for(int i=0; i<= 3 * cp.no_of_IDs_all_procs[0]-1; i++)
                   mc.my_cell_coords[i] = cp.cell_coords_all_procs[0][i]; 
           }
           else
           {
              if(mc.my_no_of_cell_IDs > 0)
              {
                  MPI_Irecv(&mc.my_cell_coords[0], 3 * mc.my_no_of_cell_IDs, MPI_DOUBLE, 0, 2222, cart_topo.mpi_cart_comm, &rcv_req);
                  MPI_Wait(&rcv_req, MPI_STATUS_IGNORE);
              }
           }
		// the above was commented out in development branch, keeping it in case its needed
           /*
            if(world.rank == 0) {
                for(int i=1; i<=world.size-1; i++) {
                    if(cp.no_of_IDs_all_procs[i] > 0) {
                        MPI_Isend(&cp.cell_types_all_procs[i][0], cp.no_of_IDs_all_procs[i], MPI_INT, i, 3333, cart_topo.mpi_cart_comm, &snd_req[i]); 
                        MPI_Wait(&snd_req[i], MPI_STATUS_IGNORE);
                    }
                }

               
               for(int i=0; i<=cp.no_of_IDs_all_procs[0]-1; i++)
                   mc.my_cell_types[i] = cp.cell_types_all_procs[0][i]; 
           }
           else {
                if(mc.my_no_of_cell_IDs > 0) {
                    MPI_Irecv(&mc.my_cell_types[0], mc.my_no_of_cell_IDs, MPI_INT, 0, 3333, cart_topo.mpi_cart_comm, &rcv_req);
                    MPI_Wait(&rcv_req, MPI_STATUS_IGNORE);
                }
           }
           */
           
           
        }
        
        int request_IDs_from_root(int requested_IDs, int *range_ID, int strt_cell_ID, mpi_Environment &world, mpi_Cartesian &cart_topo)
        {
        	/*--------------------------------------------------------------------------------------*/
          /* Declare gathering array at all processes but allocate only at root									  */
        	/*--------------------------------------------------------------------------------------*/
        	
        	int *gather_array_at_root, *scatter_array_at_root;
        	
        	if(world.rank == 0)
        	{
        		gather_array_at_root  = new int[world.size];      //A single value collected from each process
        		scatter_array_at_root = new int[2 * world.size]; 	//2 values sent to each process
        	}
        	
        	/*--------------------------------------------------------------------------------------*/
          /* Collect the number of requested IDs from each process on root array 									*/
          /*--------------------------------------------------------------------------------------*/
        	
        	MPI_Gather(&requested_IDs, 1, MPI_INT, gather_array_at_root, 1, MPI_INT, 0, cart_topo.mpi_cart_comm);
        	 
        	if(world.rank == 0)
        	{
        		for(int i=0; i<= 2 * world.size-1; i=i+2)
        		{
        			scatter_array_at_root[i] 		= strt_cell_ID;
        			scatter_array_at_root[i+1] 	= strt_cell_ID + gather_array_at_root[i/2] - 1;
        			strt_cell_ID 								= scatter_array_at_root[i+1] + 1;  
        		}
        	}
        	
        	MPI_Scatter(scatter_array_at_root, 2, MPI_INT, range_ID, 2, MPI_INT, 0, cart_topo.mpi_cart_comm);
        	return(strt_cell_ID);  
        }
    }
}
