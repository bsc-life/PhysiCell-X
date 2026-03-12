void Cell_Container::pack_cell_interact_info(mpi_Environment &world, mpi_Cartesian &cart_topo)
{

	/*--------------------------------------------------------------*/
	/* 4 variables below are data members of class Cell_Container 	*/
	/* These are used in sending packed cells also - use carefully	*/
	/*--------------------------------------------------------------*/
	position_left = 0;									//Must be initialized to 0
	position_right = 0;									//Must be initialized to 0
	snd_buf_left.resize(0); 											//When we enter function again, this is reset
	snd_buf_right.resize(0);											//Reset the others also

	int z_dim = underlying_mesh.z_coordinates.size();
	int y_dim = underlying_mesh.y_coordinates.size();
	int x_dim = underlying_mesh.x_coordinates.size();

	int len_snd_buf_left  	 = 0;
	int len_snd_buf_right 	 = 0;

	// Reserve enough space up-front to avoid repeated reallocations while packing.
	if(world.rank > 0)
	{
		size_t estimated = 0;
		for(int k=0; k<z_dim; k++)
		{
			for(int j=0; j<y_dim; j++)
			{
				uint local_vxl_index = underlying_mesh.voxel_index(0, j, k);
				int no_of_cells_in_vxl = agent_grid[local_vxl_index].size();
				estimated += 2 * sizeof(int) + 4 * sizeof(double);
				estimated += static_cast<size_t>(no_of_cells_in_vxl) * (2 * sizeof(int) + 8 * sizeof(double));
			}
		}
		snd_buf_left.reserve(estimated);
	}
	if(world.rank < world.size-1)
	{
		size_t estimated = 0;
		for(int k=0; k<z_dim; k++)
		{
			for(int j=0; j<y_dim; j++)
			{
				uint local_vxl_index = underlying_mesh.voxel_index(x_dim-1, j, k);
				int no_of_cells_in_vxl = agent_grid[local_vxl_index].size();
				estimated += 2 * sizeof(int) + 4 * sizeof(double);
				estimated += static_cast<size_t>(no_of_cells_in_vxl) * (2 * sizeof(int) + 8 * sizeof(double));
			}
		}
		snd_buf_right.reserve(estimated);
	}

	// per_cell_bytes = 2 ints (ID, type) + 1 bool (dead) + 3 doubles (position) + 1 int (number_of_nuclei) + 
	//                  8 doubles (volume vars) + n_substrates doubles (internalized_substrates) +
	//                  2 doubles (target_solid) + n_substrates doubles (fraction_transferred_when_ingested)
	//                = 3*sizeof(int) + 1*sizeof(bool) + (3 + 8 + 2)*sizeof(double) + 2*n_substrates*sizeof(double)
	//                = 3*sizeof(int) + 1*sizeof(bool) + 13*sizeof(double) + 2*n_substrates*sizeof(double)
	const int per_cell_bytes = 3 * sizeof(int) + 1 * sizeof(bool) + 13 * sizeof(double) + 2 * underlying_mesh.n_substrates * sizeof(double);

	// First pass: compute total buffer sizes to avoid repeated reallocations while packing.
	if(world.rank > 0) 
	{
		for(int j=0; j< y_dim; j++) 
		{
			for(int k=0; k<z_dim; k++)
			{
				int local_vxl_inex = underlying_mesh.voxel_index(0, j, k);
				int no_of_cells_in_vxl 	= agent_grid[local_vxl_inex].size();
				len_snd_buf_left += 2 * sizeof(int);
				len_snd_buf_left += per_cell_bytes * no_of_cells_in_vxl;
			}
		}
	}

	if(world.rank < world.size-1)
	{
		for(int j=0; j<y_dim; j++)
		{
			for(int k=0; k<z_dim; k++)
			{
				int local_vxl_inex 	= underlying_mesh.voxel_index(underlying_mesh.x_coordinates.size()-1, j, k);
				int no_of_cells_in_vxl 	= agent_grid[local_vxl_inex].size();
				len_snd_buf_right += 2 * sizeof(int);
				len_snd_buf_right += per_cell_bytes * no_of_cells_in_vxl;
			}
		}
	}

	snd_buf_left.resize(len_snd_buf_left);
	snd_buf_right.resize(len_snd_buf_right);
	position_left = 0;
	position_right = 0;

	/*---------------------*/
	/* Data Packing first  */
	/*---------------------*/

	/* Data Packing of the left boundary voxel cells of all processes except rank = 0 */

	if(world.rank > 0) 
	{
		for(int j=0; j< y_dim; j++) 
		{
			for(int k=0; k<z_dim; k++)
			{
				int local_vxl_inex = underlying_mesh.voxel_index(0, j, k);
				int global_mesh_index = underlying_mesh.voxels[local_vxl_inex].global_mesh_index;
				int no_of_cells_in_vxl 	= agent_grid[local_vxl_inex].size();

				MPI_Pack(&global_mesh_index,  1, MPI_INT,  snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
				MPI_Pack(&no_of_cells_in_vxl, 1, MPI_INT,  snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);

				for(int vec_len=0; vec_len < no_of_cells_in_vxl; vec_len++)
				{
					Cell *pCell = agent_grid[local_vxl_inex][vec_len];

					MPI_Pack(&pCell->ID, 1, MPI_INT, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->type, 1, MPI_INT, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.death.dead, 1, MPI_C_BOOL, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					if (pCell->position.size() != 3) std::cout << "Ghost cell rank " << world.rank << " ID " << pCell->ID << std::endl;
 					MPI_Pack(&pCell->position[0] , 3, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->state.number_of_nuclei, 1, MPI_INT, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.total, 1, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.cytoplasmic_fluid,  1, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.nuclear_fluid,  1, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.cytoplasmic_solid,  1, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.nuclear_solid,  1, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&(*pCell->internalized_substrates)[0],underlying_mesh.n_substrates, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.target_solid_cytoplasmic,  1, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.target_solid_nuclear,  1, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.molecular.fraction_transferred_when_ingested[0],  underlying_mesh.n_substrates, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
				}
			}
		}
	}

/* Data Packing of the right boundary voxel cells of all processes except rank = world.size-1 */

if(world.rank < world.size-1)
	{
		for(int j=0; j<y_dim; j++)
		{
			for(int k=0; k<z_dim; k++)
			{
				int local_vxl_inex 	= underlying_mesh.voxel_index(underlying_mesh.x_coordinates.size()-1, j, k);
				int global_mesh_index = underlying_mesh.voxels[local_vxl_inex].global_mesh_index;
				int no_of_cells_in_vxl 	= agent_grid[local_vxl_inex].size();

				MPI_Pack(&global_mesh_index,  1, MPI_INT, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
				MPI_Pack(&no_of_cells_in_vxl, 1, MPI_INT, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);

				for(int vec_len=0; vec_len < no_of_cells_in_vxl; vec_len++)
				{
					Cell *pCell = agent_grid[local_vxl_inex][vec_len];

					MPI_Pack(&pCell->ID, 1, MPI_INT, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->type, 1, MPI_INT, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.death.dead, 1, MPI_C_BOOL, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					if (pCell->position.size() != 3) std::cout << "Ghost cell rank " << world.rank << " ID " << pCell->ID << std::endl;
					MPI_Pack(&pCell->position[0] , 3, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->state.number_of_nuclei, 1, MPI_INT, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.total, 1, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.cytoplasmic_fluid,  1, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.nuclear_fluid,  1, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.cytoplasmic_solid,  1, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.nuclear_solid,  1, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&(*pCell->internalized_substrates)[0],underlying_mesh.n_substrates, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.target_solid_cytoplasmic,  1, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.target_solid_nuclear,  1, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.molecular.fraction_transferred_when_ingested[0], underlying_mesh.n_substrates, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD );
				}
			}
		}
	}

	/*---------------------------------------------------*/
	/* Sending and Receiving of buffers across processes */
	/* Send to MPI_PROC_NULL as well, don't use 'if'		 */
	/*---------------------------------------------------*/

	int size_of_data_recvd_from_right_process = 0;
	int size_of_data_recvd_from_left_process  = 0;
	MPI_Request snd_req[2], rcv_req[2];

	/* Send to left, Receive from right: MPI_PROC_NULL<----R0<-----R1<----R2<----R3 */

	MPI_Irecv(&size_of_data_recvd_from_right_process, 1, MPI_INT, cart_topo.X_RIGHT, 1111, cart_topo.mpi_cart_comm, &rcv_req[0]);
	MPI_Isend(&position_left,  1, MPI_INT, cart_topo.X_LEFT,  1111, cart_topo.mpi_cart_comm, &snd_req[0]);

	/* Send to right, Receive from left: R0----->R1---->R2---->R3--->MPI_PROC_NULL */

	MPI_Irecv(&size_of_data_recvd_from_left_process, 1, MPI_INT, cart_topo.X_LEFT, 2222, cart_topo.mpi_cart_comm, &rcv_req[1]);
	MPI_Isend(&position_right,  1, MPI_INT, cart_topo.X_RIGHT, 2222, cart_topo.mpi_cart_comm, &snd_req[1]);

	MPI_Waitall(2, snd_req, MPI_STATUSES_IGNORE);
	MPI_Waitall(2, rcv_req, MPI_STATUSES_IGNORE);

	/* Resize the actual buffers that will contain the data */

	if(world.rank < world.size-1)
		rcv_buf_right.resize(size_of_data_recvd_from_right_process);

	if(world.rank > 0)
		rcv_buf_left.resize(size_of_data_recvd_from_left_process);

	/* Now send the actual data in snd_buf_left and snd_buf_right */

	char* rcv_right_ptr = rcv_buf_right.empty() ? nullptr : rcv_buf_right.data();
	char* rcv_left_ptr  = rcv_buf_left.empty()  ? nullptr : rcv_buf_left.data();
	char* snd_left_ptr  = snd_buf_left.empty()  ? nullptr : snd_buf_left.data();
	char* snd_right_ptr = snd_buf_right.empty() ? nullptr : snd_buf_right.data();

	MPI_Irecv(rcv_right_ptr, size_of_data_recvd_from_right_process , MPI_PACKED, cart_topo.X_RIGHT, 3333, cart_topo.mpi_cart_comm, &rcv_req[0]);
	MPI_Isend(snd_left_ptr, len_snd_buf_left, MPI_PACKED, cart_topo.X_LEFT,  3333, cart_topo.mpi_cart_comm, &snd_req[0]);

	MPI_Irecv(rcv_left_ptr, size_of_data_recvd_from_left_process, MPI_PACKED, cart_topo.X_LEFT, 4444, cart_topo.mpi_cart_comm, &rcv_req[1]);
	MPI_Isend(snd_right_ptr, len_snd_buf_right, MPI_PACKED, cart_topo.X_RIGHT, 4444, cart_topo.mpi_cart_comm, &snd_req[1]);

	MPI_Waitall(2, snd_req, MPI_STATUSES_IGNORE);
	MPI_Waitall(2, rcv_req, MPI_STATUSES_IGNORE);

	/*--------------------------------------------------------------------------*/
	/* Now tackle unpacking of cells in a special class called Moore_Voxel_Info */
	/*--------------------------------------------------------------------------*/

	/* Declare list of Moore_Voxel_Info - one element each for every boundary voxel i.e. y_dim * z_dim voxels */

	/* std::vector<Moore_Voxel_Info> mbfr (moore boundary from right [process]), similarly mbfl */
	/* is declared in PhysiCell_cell_container.h as a public data member 												*/

	ivfr.resize(y_dim * z_dim);
	ivfl.resize(y_dim * z_dim);

	/* Additionally, we have um_mbfl = unordered map moore boundary from left ---> 								 */
	/* maps the global mesh index of a voxel 																											 */
	/* to a specific Moore_Voxel_Info object. The idea is to take the global mesh index from 			 */
	/* moore_connected_voxel_global_indices_left/right and use this as a key to hash into the Moore*/
	/* Voxel_Info object associated with this global_mesh_index in O(1) time											 */
	/* No need to resize() um_mbfl and um_mbfr 																										 */


	position_right = 0;		//Was used in Packing but re-used here, its a data member (remember)
	position_left = 0; 		//Was used in Packing but r-used here, its a data member (remember)

	int size_right = size_of_data_recvd_from_right_process;	//For convenience
	int size_left =  size_of_data_recvd_from_left_process; 	//For convenience

	um_ivfr.clear();
	um_ivfl.clear();


	/* first unpack data coming FROM right process in 'rcv_buf_right' into 'mbfr' 					 */
	/* Outer loop runs till y_dim * z_dim voxels and each voxel can have any number of cells */
	/* IMPORTANT: rcv_buf_right will be unpacked ONLY on ranks < world.size-1 							 */
	/* REMEMBER: size_right = 0 on rank=world.size-1, hence error 													 */

	if(world.rank < world.size-1)
	{
		for(int vxl_ctr=0; vxl_ctr<y_dim*z_dim; vxl_ctr++)
		{
			MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].global_mesh_index, 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].no_of_cells_in_vxl, 1, MPI_INT, MPI_COMM_WORLD);

			/* Now resize vector moore_cells using no_of_cells_in_vxl unpacked above */

			int no_of_cells_in_vxl = ivfr[vxl_ctr].no_of_cells_in_vxl;
			ivfr[vxl_ctr].cells.resize(no_of_cells_in_vxl);
			//std::cout << "Rank " << world.rank << " rcv voxel " <<  ivfr[vxl_ctr].global_mesh_index << " with " << ivfr[vxl_ctr].no_of_cells_in_vxl << " cells" << std::endl;
			for(int cell_ctr=0; cell_ctr<no_of_cells_in_vxl; cell_ctr++)
			{
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].ID, 1, MPI_INT, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].type, 1, MPI_INT, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].dead, 1, MPI_C_BOOL, MPI_COMM_WORLD);
				ivfr[vxl_ctr].cells[cell_ctr].position.resize(3);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].position[0], 3, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].number_of_nuclei, 1, MPI_INT, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].phenotype_volume, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].cytoplasmic_fluid, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].nuclear_fluid, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].cytoplasmic_solid, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].nuclear_solid, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				ivfr[vxl_ctr].cells[cell_ctr].internalized_substrates.resize(underlying_mesh.n_substrates);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].internalized_substrates[0], underlying_mesh.n_substrates, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].target_solid_cytoplasmic, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].target_solid_nuclear, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				ivfr[vxl_ctr].cells[cell_ctr].fraction_transferred_when_ingested.resize(underlying_mesh.n_substrates);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].fraction_transferred_when_ingested[0], underlying_mesh.n_substrates, MPI_DOUBLE, MPI_COMM_WORLD);
				
			}

			  um_ivfr[ivfr[vxl_ctr].global_mesh_index]=ivfr[vxl_ctr];

		}
	}

	/* Now unpack data coming FROM left process in 'rcv_buf_left' into 'mbfl' 					 		 */
	/* Outer loop runs till y_dim * z_dim voxels and each voxel can have any number of cells */
	/* IMPORTANT: rcv_buf_left will be unpacked ONLY on ranks > 0 							 						 */
	/* REMEMBER: size_left = 0 on rank=0, hence error 													 						 */

	if(world.rank > 0)
	{
		for(int vxl_ctr=0; vxl_ctr<y_dim*z_dim; vxl_ctr++)
		{
			MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].global_mesh_index, 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].no_of_cells_in_vxl, 1, MPI_INT, MPI_COMM_WORLD);

			/* Now resize vector moore_cells using no_of_cells_in_vxl unpacked above */

			int no_of_cells_in_vxl = ivfl[vxl_ctr].no_of_cells_in_vxl;
			ivfl[vxl_ctr].cells.resize(no_of_cells_in_vxl);
			//std::cout << "Rank " << world.rank << " rcv voxel " <<  ivfl[vxl_ctr].global_mesh_index << " with " << ivfl[vxl_ctr].no_of_cells_in_vxl << " cells" << std::endl;
			for(int cell_ctr=0; cell_ctr<no_of_cells_in_vxl; cell_ctr++)
			{
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].ID, 1, MPI_INT, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].type, 1, MPI_INT, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].dead, 1, MPI_C_BOOL, MPI_COMM_WORLD);
				ivfl[vxl_ctr].cells[cell_ctr].position.resize(3);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].position[0], 3, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].number_of_nuclei, 1, MPI_INT, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].phenotype_volume, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].cytoplasmic_fluid, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].nuclear_fluid, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].cytoplasmic_solid, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].nuclear_solid, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				ivfl[vxl_ctr].cells[cell_ctr].internalized_substrates.resize(underlying_mesh.n_substrates);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].internalized_substrates[0], underlying_mesh.n_substrates, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].target_solid_cytoplasmic, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].target_solid_nuclear, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				ivfl[vxl_ctr].cells[cell_ctr].fraction_transferred_when_ingested.resize(underlying_mesh.n_substrates);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].fraction_transferred_when_ingested[0], underlying_mesh.n_substrates, MPI_DOUBLE, MPI_COMM_WORLD);
			}

			um_ivfl[ivfl[vxl_ctr].global_mesh_index]=ivfl[vxl_ctr];
		}
	}
	//std::cout << "Rank " << world.rank << " has received cell information" << std::endl;
	/* Now one step is left: Need to map the global mesh index to a specific Moor_Voxel_Info object */
	/* i.e. declare map_glbl_indx_to_Moore_Voxel_Info_object[global_mesh_index]=mbfl[vxl_ctr] 			*/
	/* Similarly for mbfr as well 																																	*/
	/* Declare: std::unordered_map<int, Moore_Voxel_Info> in PhysiCell_cell_container.h 						*/
}

void Cell_Container::pack_cell_interact_info(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	position_left = 0;
	position_right = 0;
	snd_buf_left.resize(0);
	snd_buf_right.resize(0);

	std::vector<vector<char>> snd_buf_left_per_voxel(y_dim * z_dim);
	std::vector<vector<char>> snd_buf_right_per_voxel(y_dim * z_dim);

	int z_dim = underlying_mesh.z_coordinates.size();
	int y_dim = underlying_mesh.y_coordinates.size();
	int x_dim = underlying_mesh.x_coordinates.size();

	int len_snd_buf_left  	 = 0;
	int len_snd_buf_right 	 = 0;

	auto start = std::chrono::high_resolution_clock::now();

	const size_t header_voxel = 2 * sizeof(int); // global_mesh_index
	const size_t per_cell_bytes = 3 * sizeof(int) + 1 * sizeof(bool) + 10 * sizeof(double) + 2 * underlying_mesh.n_substrates * sizeof(double);
	// Reserve enough space up-front to avoid repeated reallocations while packing.
	if (world.rank > 0) { //Send information to left neighbor
		size_t estimated = 0;
		size_t cells = 0;
		#pragma omp parallel for collapse(2) reduction(+:cells)
		for(int k=0; k<z_dim; k++)
			for(int j=0; j<y_dim; j++)
			{
				uint local_vxl_index = underlying_mesh.voxel_index(0, j, k);
				cells += agent_grid[local_vxl_index].size();
			}
		estimated += static_cast<size_t>(y_dim * z_dim) * header_voxel;
		estimated += cells * per_cell_bytes;
		snd_buf_left.resize(estimated);
		len_snd_buf_left = estimated; 
	}

	if(world.rank < world.size-1){ //Send information to right neighbor
		size_t estimated = 0;
		size_t cells = 0;
		#pragma omp parallel for collapse(2) reduction(+:cells)
		for(int k=0; k<z_dim; k++)
			for(int j=0; j<y_dim; j++)
			{
				uint local_vxl_index = underlying_mesh.voxel_index(x_dim-1, j, k);
				cells += agent_grid[local_vxl_index].size();
			}
		estimated += static_cast<size_t>(y_dim * z_dim) * header_voxel;
		estimated += cells * per_cell_bytes;
		snd_buf_right.resize(estimated);
		len_snd_buf_right = estimated;
	}

	if(world.rank > 0) 
	{
		
		for(int j=0; j< y_dim; j++) 
		{
			for(int k=0; k<z_dim; k++)
			{
				int thread_position = ()
				int local_vxl_inex = underlying_mesh.voxel_index(0, j, k);
				int global_mesh_index = underlying_mesh.voxels[local_vxl_inex].global_mesh_index;
				int no_of_cells_in_vxl 	= agent_grid[local_vxl_inex].size();

				MPI_Pack(&global_mesh_index,  1, MPI_INT,  snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
				MPI_Pack(&no_of_cells_in_vxl, 1, MPI_INT,  snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);

				for(int vec_len=0; vec_len < no_of_cells_in_vxl; vec_len++)
				{
					Cell *pCell = agent_grid[local_vxl_inex][vec_len];

					MPI_Pack(&pCell->ID, 1, MPI_INT, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->type, 1, MPI_INT, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.death.dead, 1, MPI_C_BOOL, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					if (pCell->position.size() != 3) std::cout << "Ghost cell rank " << world.rank << " ID " << pCell->ID << std::endl;
 					MPI_Pack(&pCell->position[0] , 3, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->state.number_of_nuclei, 1, MPI_INT, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.total, 1, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.cytoplasmic_fluid,  1, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.nuclear_fluid,  1, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.cytoplasmic_solid,  1, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.nuclear_solid,  1, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&(*pCell->internalized_substrates)[0],underlying_mesh.n_substrates, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.target_solid_cytoplasmic,  1, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.target_solid_nuclear,  1, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.molecular.fraction_transferred_when_ingested[0],  underlying_mesh.n_substrates, MPI_DOUBLE, snd_buf_left.data(), len_snd_buf_left, &position_left, MPI_COMM_WORLD);
				}
			}
		}
	}

	if(world.rank < world.size-1)
	{
		for(int j=0; j<y_dim; j++)
		{
			for(int k=0; k<z_dim; k++)
			{
				int local_vxl_inex 	= underlying_mesh.voxel_index(underlying_mesh.x_coordinates.size()-1, j, k);
				int global_mesh_index = underlying_mesh.voxels[local_vxl_inex].global_mesh_index;
				int no_of_cells_in_vxl 	= agent_grid[local_vxl_inex].size();

				MPI_Pack(&global_mesh_index,  1, MPI_INT, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
				MPI_Pack(&no_of_cells_in_vxl, 1, MPI_INT, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);

				for(int vec_len=0; vec_len < no_of_cells_in_vxl; vec_len++)
				{
					Cell *pCell = agent_grid[local_vxl_inex][vec_len];

					MPI_Pack(&pCell->ID, 1, MPI_INT, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->type, 1, MPI_INT, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.death.dead, 1, MPI_C_BOOL, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					if (pCell->position.size() != 3) std::cout << "Ghost cell rank " << world.rank << " ID " << pCell->ID << std::endl;
					MPI_Pack(&pCell->position[0] , 3, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->state.number_of_nuclei, 1, MPI_INT, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.total, 1, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.cytoplasmic_fluid,  1, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.nuclear_fluid,  1, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.cytoplasmic_solid,  1, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.nuclear_solid,  1, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&(*pCell->internalized_substrates)[0],underlying_mesh.n_substrates, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.target_solid_cytoplasmic,  1, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.volume.target_solid_nuclear,  1, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.molecular.fraction_transferred_when_ingested[0], underlying_mesh.n_substrates, MPI_DOUBLE, snd_buf_right.data(), len_snd_buf_right, &position_right, MPI_COMM_WORLD );
				}
			}
		}
	}



	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	std::cout << "[Rank " << world.rank << "]:  pack cell interaction info: " << duration.count() << " microseconds" << std::endl;
	/*---------------------------------------------------*/
	/* Sending and Receiving of buffers across processes */
	/* Send to MPI_PROC_NULL as well, don't use 'if'		 */
	/*---------------------------------------------------*/

	int size_of_data_recvd_from_right_process = 0;
	int size_of_data_recvd_from_left_process  = 0;
	MPI_Request snd_req[2], rcv_req[2];

	start = std::chrono::high_resolution_clock::now();

	/* Send to left, Receive from right: MPI_PROC_NULL<----R0<-----R1<----R2<----R3 */

	MPI_Irecv(&size_of_data_recvd_from_right_process, 1, MPI_INT, cart_topo.X_RIGHT, 1111, cart_topo.mpi_cart_comm, &rcv_req[0]);
	MPI_Isend(&position_left,  1, MPI_INT, cart_topo.X_LEFT,  1111, cart_topo.mpi_cart_comm, &snd_req[0]);

	/* Send to right, Receive from left: R0----->R1---->R2---->R3--->MPI_PROC_NULL */

	MPI_Irecv(&size_of_data_recvd_from_left_process, 1, MPI_INT, cart_topo.X_LEFT, 2222, cart_topo.mpi_cart_comm, &rcv_req[1]);
	MPI_Isend(&position_right,  1, MPI_INT, cart_topo.X_RIGHT, 2222, cart_topo.mpi_cart_comm, &snd_req[1]);

	MPI_Waitall(2, snd_req, MPI_STATUSES_IGNORE);
	MPI_Waitall(2, rcv_req, MPI_STATUSES_IGNORE);

	/* Resize the actual buffers that will contain the data */

	if(world.rank < world.size-1)
		rcv_buf_right.resize(size_of_data_recvd_from_right_process);

	if(world.rank > 0)
		rcv_buf_left.resize(size_of_data_recvd_from_left_process);

	/* Now send the actual data in snd_buf_left and snd_buf_right */

	char* rcv_right_ptr = rcv_buf_right.empty() ? nullptr : rcv_buf_right.data();
	char* rcv_left_ptr  = rcv_buf_left.empty()  ? nullptr : rcv_buf_left.data();
	char* snd_left_ptr  = snd_buf_left.empty()  ? nullptr : snd_buf_left.data();
	char* snd_right_ptr = snd_buf_right.empty() ? nullptr : snd_buf_right.data();

	MPI_Irecv(rcv_right_ptr, size_of_data_recvd_from_right_process , MPI_PACKED, cart_topo.X_RIGHT, 3333, cart_topo.mpi_cart_comm, &rcv_req[0]);
	MPI_Isend(snd_left_ptr, len_snd_buf_left, MPI_PACKED, cart_topo.X_LEFT,  3333, cart_topo.mpi_cart_comm, &snd_req[0]);

	MPI_Irecv(rcv_left_ptr, size_of_data_recvd_from_left_process, MPI_PACKED, cart_topo.X_LEFT, 4444, cart_topo.mpi_cart_comm, &rcv_req[1]);
	MPI_Isend(snd_right_ptr, len_snd_buf_right, MPI_PACKED, cart_topo.X_RIGHT, 4444, cart_topo.mpi_cart_comm, &snd_req[1]);

	MPI_Waitall(2, snd_req, MPI_STATUSES_IGNORE);
	MPI_Waitall(2, rcv_req, MPI_STATUSES_IGNORE);

	end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	std::cout << "[Rank " << world.rank << "]:  transfer cell interaction info: " << duration.count() << " microseconds" << std::endl;

	/*--------------------------------------------------------------------------*/
	/* Now tackle unpacking of cells in a special class called Moore_Voxel_Info */
	/*--------------------------------------------------------------------------*/

	/* Declare list of Moore_Voxel_Info - one element each for every boundary voxel i.e. y_dim * z_dim voxels */

	/* std::vector<Moore_Voxel_Info> mbfr (moore boundary from right [process]), similarly mbfl */
	/* is declared in PhysiCell_cell_container.h as a public data member 												*/
	start = std::chrono::high_resolution_clock::now();

	ivfr.resize(y_dim * z_dim);
	ivfl.resize(y_dim * z_dim);

	/* Additionally, we have um_mbfl = unordered map moore boundary from left ---> 								 */
	/* maps the global mesh index of a voxel 																											 */
	/* to a specific Moore_Voxel_Info object. The idea is to take the global mesh index from 			 */
	/* moore_connected_voxel_global_indices_left/right and use this as a key to hash into the Moore*/
	/* Voxel_Info object associated with this global_mesh_index in O(1) time											 */
	/* No need to resize() um_mbfl and um_mbfr 																										 */


	position_right = 0;		//Was used in Packing but re-used here, its a data member (remember)
	position_left = 0; 		//Was used in Packing but r-used here, its a data member (remember)

	int size_right = size_of_data_recvd_from_right_process;	//For convenience
	int size_left =  size_of_data_recvd_from_left_process; 	//For convenience

	um_ivfr.clear();
	um_ivfl.clear();


	/* first unpack data coming FROM right process in 'rcv_buf_right' into 'mbfr' 					 */
	/* Outer loop runs till y_dim * z_dim voxels and each voxel can have any number of cells */
	/* IMPORTANT: rcv_buf_right will be unpacked ONLY on ranks < world.size-1 							 */
	/* REMEMBER: size_right = 0 on rank=world.size-1, hence error 													 */

	if(world.rank < world.size-1)
	{
		for(int vxl_ctr=0; vxl_ctr<y_dim*z_dim; vxl_ctr++)
		{
			MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].global_mesh_index, 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].no_of_cells_in_vxl, 1, MPI_INT, MPI_COMM_WORLD);

			/* Now resize vector moore_cells using no_of_cells_in_vxl unpacked above */

			int no_of_cells_in_vxl = ivfr[vxl_ctr].no_of_cells_in_vxl;
			ivfr[vxl_ctr].cells.resize(no_of_cells_in_vxl);
			//std::cout << "Rank " << world.rank << " rcv voxel " <<  ivfr[vxl_ctr].global_mesh_index << " with " << ivfr[vxl_ctr].no_of_cells_in_vxl << " cells" << std::endl;
			for(int cell_ctr=0; cell_ctr<no_of_cells_in_vxl; cell_ctr++)
			{
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].ID, 1, MPI_INT, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].type, 1, MPI_INT, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].dead, 1, MPI_C_BOOL, MPI_COMM_WORLD);
				ivfr[vxl_ctr].cells[cell_ctr].position.resize(3);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].position[0], 3, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].number_of_nuclei, 1, MPI_INT, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].phenotype_volume, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].cytoplasmic_fluid, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].nuclear_fluid, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].cytoplasmic_solid, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].nuclear_solid, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				ivfr[vxl_ctr].cells[cell_ctr].internalized_substrates.resize(underlying_mesh.n_substrates);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].internalized_substrates[0], underlying_mesh.n_substrates, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].target_solid_cytoplasmic, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].target_solid_nuclear, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				ivfr[vxl_ctr].cells[cell_ctr].fraction_transferred_when_ingested.resize(underlying_mesh.n_substrates);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &ivfr[vxl_ctr].cells[cell_ctr].fraction_transferred_when_ingested[0], underlying_mesh.n_substrates, MPI_DOUBLE, MPI_COMM_WORLD);
				
			}

			  um_ivfr[ivfr[vxl_ctr].global_mesh_index]=ivfr[vxl_ctr];

		}
	}

	/* Now unpack data coming FROM left process in 'rcv_buf_left' into 'mbfl' 					 		 */
	/* Outer loop runs till y_dim * z_dim voxels and each voxel can have any number of cells */
	/* IMPORTANT: rcv_buf_left will be unpacked ONLY on ranks > 0 							 						 */
	/* REMEMBER: size_left = 0 on rank=0, hence error 													 						 */

	if(world.rank > 0)
	{
		for(int vxl_ctr=0; vxl_ctr<y_dim*z_dim; vxl_ctr++)
		{
			MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].global_mesh_index, 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].no_of_cells_in_vxl, 1, MPI_INT, MPI_COMM_WORLD);

			/* Now resize vector moore_cells using no_of_cells_in_vxl unpacked above */

			int no_of_cells_in_vxl = ivfl[vxl_ctr].no_of_cells_in_vxl;
			ivfl[vxl_ctr].cells.resize(no_of_cells_in_vxl);
			//std::cout << "Rank " << world.rank << " rcv voxel " <<  ivfl[vxl_ctr].global_mesh_index << " with " << ivfl[vxl_ctr].no_of_cells_in_vxl << " cells" << std::endl;
			for(int cell_ctr=0; cell_ctr<no_of_cells_in_vxl; cell_ctr++)
			{
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].ID, 1, MPI_INT, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].type, 1, MPI_INT, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].dead, 1, MPI_C_BOOL, MPI_COMM_WORLD);
				ivfl[vxl_ctr].cells[cell_ctr].position.resize(3);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].position[0], 3, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].number_of_nuclei, 1, MPI_INT, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].phenotype_volume, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].cytoplasmic_fluid, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].nuclear_fluid, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].cytoplasmic_solid, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].nuclear_solid, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				ivfl[vxl_ctr].cells[cell_ctr].internalized_substrates.resize(underlying_mesh.n_substrates);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].internalized_substrates[0], underlying_mesh.n_substrates, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].target_solid_cytoplasmic, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].target_solid_nuclear, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				ivfl[vxl_ctr].cells[cell_ctr].fraction_transferred_when_ingested.resize(underlying_mesh.n_substrates);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &ivfl[vxl_ctr].cells[cell_ctr].fraction_transferred_when_ingested[0], underlying_mesh.n_substrates, MPI_DOUBLE, MPI_COMM_WORLD);
			}

			um_ivfl[ivfl[vxl_ctr].global_mesh_index]=ivfl[vxl_ctr];
		}
	}

	end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	std::cout << "[Rank " << world.rank << "]:  unpack cell interaction info: " << duration.count() << " microseconds" << std::endl;
	//std::cout << "Rank " << world.rank << " has received cell information" << std::endl;
	/* Now one step is left: Need to map the global mesh index to a specific Moor_Voxel_Info object */
	/* i.e. declare map_glbl_indx_to_Moore_Voxel_Info_object[global_mesh_index]=mbfl[vxl_ctr] 			*/
	/* Similarly for mbfr as well 																																	*/
	/* Declare: std::unordered_map<int, Moore_Voxel_Info> in PhysiCell_cell_container.h 						*/
	
}