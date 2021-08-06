/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "../BioFVM/BioFVM_agent_container.h"
#include "PhysiCell_constants.h"
#include "../BioFVM/BioFVM_vector.h"
#include "PhysiCell_cell.h"


using namespace BioFVM;


namespace PhysiCell{

std::vector<Cell*> *all_cells;

Cell_Container::Cell_Container()
{
	all_cells = (std::vector<Cell*> *) &all_basic_agents;
	boundary_condition_for_pushed_out_agents= PhysiCell_constants::default_boundary_condition_for_pushed_out_agents;
	std::vector<Cell*> cells_ready_to_divide;
	std::vector<Cell*> cells_ready_to_die;

	return;
}

void Cell_Container::initialize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double voxel_size)
{
	initialize(x_start, x_end, y_start, y_end, z_start, z_end , voxel_size, voxel_size, voxel_size);

	return;
}

/*-------------------------------------------------*/
/* Parallel version of initialize() function above */
/*-------------------------------------------------*/

void Cell_Container::initialize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double voxel_size, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	initialize(x_start, x_end, y_start, y_end, z_start, z_end , voxel_size, voxel_size, voxel_size, world, cart_topo);

	return;
}


void Cell_Container::initialize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx, double dy, double dz)
{
	all_cells = (std::vector<Cell*> *) &all_basic_agents;
	boundary_condition_for_pushed_out_agents= PhysiCell_constants::default_boundary_condition_for_pushed_out_agents;
	std::vector<Cell*> cells_ready_to_divide;
	std::vector<Cell*> cells_ready_to_die;

	underlying_mesh.resize(x_start, x_end, y_start, y_end, z_start, z_end , dx, dy, dz);
	agent_grid.resize(underlying_mesh.voxels.size());
	max_cell_interactive_distance_in_voxel.resize(underlying_mesh.voxels.size(), 0.0);
	agents_in_outer_voxels.resize(6);

	return;
}

/*-------------------------------------------------*/
/* Parallel version of initialize() function above */
/*-------------------------------------------------*/

void Cell_Container::initialize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx, double dy, double dz, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	all_cells = (std::vector<Cell*> *) &all_basic_agents;
	boundary_condition_for_pushed_out_agents= PhysiCell_constants::default_boundary_condition_for_pushed_out_agents;
	std::vector<Cell*> cells_ready_to_divide;
	std::vector<Cell*> cells_ready_to_die;

	underlying_mesh.resize(x_start, x_end, y_start, y_end, z_start, z_end , dx, dy, dz, world, cart_topo);
	agent_grid.resize(underlying_mesh.voxels.size());
	max_cell_interactive_distance_in_voxel.resize(underlying_mesh.voxels.size(), 0.0);
	agents_in_outer_voxels.resize(6);

	return;
}

void Cell_Container::update_all_cells(double t)
{
	// update_all_cells(t, dt_settings.cell_cycle_dt_default, dt_settings.mechanics_dt_default);

	update_all_cells(t, phenotype_dt, mechanics_dt , diffusion_dt );

	return;
}

/*----------------------------------------*/
/* Parallel version of the function above	*/
/*----------------------------------------*/

void Cell_Container::update_all_cells(double t, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	// update_all_cells(t, dt_settings.cell_cycle_dt_default, dt_settings.mechanics_dt_default);

	update_all_cells(t, phenotype_dt, mechanics_dt , diffusion_dt, world, cart_topo );

	return;
}

// deprecate me JULY 2017
void Cell_Container::update_all_cells(double t, double dt)
{
	std::cout << "WARNING : " << __FUNCTION__ << " in " << __FILE__ << " is deprecated."
		<< "\tIt returns without execution." << std::endl;
	return;
	update_all_cells(t, dt,dt);

	return;
}

// deprecate me JULY 2017
void Cell_Container::update_all_cells(double t, double phenotype_dt_ , double mechanics_dt_ )
{
	std::cout << "WARNING : " << __FUNCTION__ << " in " << __FILE__ << " is deprecated."
		<< "\tIt returns without execution." << std::endl;
	return;

	//if it is the time for running cell cycle, do it!
	double time_since_last_cycle= t- last_cell_cycle_time;

	static double phenotype_tolerance = 0.001 * phenotype_dt_;
	static double mechanics_tolerance = 0.001 * mechanics_dt_;

	if( fabs(time_since_last_cycle- phenotype_dt_ ) < phenotype_tolerance || !initialzed)
	{
		// Reset the max_radius in each voxel. It will be filled in set_total_volume
		// It might be better if we calculate it before mechanics each time
		std::fill(max_cell_interactive_distance_in_voxel.begin(), max_cell_interactive_distance_in_voxel.end(), 0.0);

		if(!initialzed)
		{
			time_since_last_cycle = phenotype_dt_;
		}

		// old functions prior to 1.2.1
		/*
		#pragma omp parallel for
		for( int i=0; i < (*all_cells).size(); i++ )
		{
			if((*all_cells)[i]->is_out_of_domain)
				continue;
			(*all_cells)[i]->update_cell_and_death_parameters((*all_cells)[i],cell_cycle_dt);
			(*all_cells)[i]->advance_cell_current_phase((*all_cells)[i],time_since_last_cycle);
			(*all_cells)[i]->update_volume((*all_cells)[i], time_since_last_cycle );
		}
		*/

		// new as of 1.2.1 -- bundles cell phenotype parameter update, volume update, geometry update,
		// checking for death, and advancing the cell cycle. Not motility, though. (that's in mechanics)
		#pragma omp parallel for
		for( int i=0; i < (*all_cells).size(); i++ )
		{
			if( (*all_cells)[i]->is_out_of_domain == false )
			{
				(*all_cells)[i]->advance_bundled_phenotype_functions( time_since_last_cycle );
			}
		}

		// process divides / removes
		for( int i=0; i < cells_ready_to_divide.size(); i++ )
		{
			cells_ready_to_divide[i]->divide();
		}
		for( int i=0; i < cells_ready_to_die.size(); i++ )
		{
			cells_ready_to_die[i]->die();
		}
		num_divisions_in_current_step+=  cells_ready_to_divide.size();
		num_deaths_in_current_step+=  cells_ready_to_die.size();

		cells_ready_to_die.clear();
		cells_ready_to_divide.clear();
		last_cell_cycle_time= t;
	}

	double time_since_last_mechanics= t- last_mechanics_time;

	// if( time_since_last_mechanics>= mechanics_dt || !initialzed)
	if( fabs(time_since_last_mechanics - mechanics_dt_)< mechanics_tolerance || !initialzed)
	{
		if(!initialzed)
		{
			time_since_last_mechanics = mechanics_dt_;
		}
		// Compute velocities
		#pragma omp parallel for
		for( int i=0; i < (*all_cells).size(); i++ )
		{

			if(!(*all_cells)[i]->is_out_of_domain && (*all_cells)[i]->is_movable)
			{
				// update_velocity already includes the motility update
				//(*all_cells)[i]->phenotype.motility.update_motility_vector( (*all_cells)[i] ,(*all_cells)[i]->phenotype , time_since_last_mechanics );
				(*all_cells)[i]->functions.update_velocity( (*all_cells)[i], (*all_cells)[i]->phenotype, time_since_last_mechanics);
			}

			if( (*all_cells)[i]->functions.custom_cell_rule )
			{
				(*all_cells)[i]->functions.custom_cell_rule((*all_cells)[i], (*all_cells)[i]->phenotype, time_since_last_mechanics);
			}
		}
		// Calculate new positions
		#pragma omp parallel for
		for( int i=0; i < (*all_cells).size(); i++ )
		{
			if(!(*all_cells)[i]->is_out_of_domain && (*all_cells)[i]->is_movable)
			{
				(*all_cells)[i]->update_position(time_since_last_mechanics);
			}
		}

		// When somebody reviews this code, let's add proper braces for clarity!!!

		// Update cell indices in the container
		for( int i=0; i < (*all_cells).size(); i++ )
			if(!(*all_cells)[i]->is_out_of_domain && (*all_cells)[i]->is_movable)
				(*all_cells)[i]->update_voxel_in_container();
		last_mechanics_time=t;
	}
	initialzed=true;
	return;
}

void Cell_Container::update_all_cells(double t, double phenotype_dt_ , double mechanics_dt_ , double diffusion_dt_ )
{
	// secretions and uptakes. Syncing with BioFVM is automated.


	#pragma omp parallel for
	for( int i=0; i < (*all_cells).size(); i++ )
	{
		(*all_cells)[i]->phenotype.secretion.advance( (*all_cells)[i], (*all_cells)[i]->phenotype , diffusion_dt_ );
	}

	//if it is the time for running cell cycle, do it!
	double time_since_last_cycle= t- last_cell_cycle_time;

	static double phenotype_dt_tolerance = 0.001 * phenotype_dt_;
	static double mechanics_dt_tolerance = 0.001 * mechanics_dt_;

	if( fabs(time_since_last_cycle-phenotype_dt_ ) < phenotype_dt_tolerance || !initialzed)
	{
		// Reset the max_radius in each voxel. It will be filled in set_total_volume
		// It might be better if we calculate it before mechanics each time
		std::fill(max_cell_interactive_distance_in_voxel.begin(), max_cell_interactive_distance_in_voxel.end(), 0.0);

		if(!initialzed)
		{
			time_since_last_cycle = phenotype_dt_;
		}

		// new as of 1.2.1 -- bundles cell phenotype parameter update, volume update, geometry update,
		// checking for death, and advancing the cell cycle. Not motility, though. (that's in mechanics)
		#pragma omp parallel for
		for( int i=0; i < (*all_cells).size(); i++ )
		{
			if( (*all_cells)[i]->is_out_of_domain == false )
			{
				(*all_cells)[i]->advance_bundled_phenotype_functions( time_since_last_cycle );
			}
		}


		// process divides / removes

		for( int i=0; i < cells_ready_to_divide.size(); i++ )
		{
			cells_ready_to_divide[i]->divide();
		}
		for( int i=0; i < cells_ready_to_die.size(); i++ )
		{
			cells_ready_to_die[i]->die();
		}
		num_divisions_in_current_step+=  cells_ready_to_divide.size();
		num_deaths_in_current_step+=  cells_ready_to_die.size();

		cells_ready_to_die.clear();
		cells_ready_to_divide.clear();
		last_cell_cycle_time= t;
	}

	double time_since_last_mechanics= t- last_mechanics_time;

	// if( time_since_last_mechanics>= mechanics_dt || !initialzed)
	if( fabs(time_since_last_mechanics - mechanics_dt_) < mechanics_dt_tolerance || !initialzed)
	{
		if(!initialzed)
		{
			time_since_last_mechanics = mechanics_dt_;
		}

		// new February 2018
		// if we need gradients, compute them
		if( default_microenvironment_options.calculate_gradients )
		{ microenvironment.compute_all_gradient_vectors();  }
		// end of new in Feb 2018

		// Compute velocities
		#pragma omp parallel for
		for( int i=0; i < (*all_cells).size(); i++ )
		{

			if(!(*all_cells)[i]->is_out_of_domain && (*all_cells)[i]->is_movable && (*all_cells)[i]->functions.update_velocity )
			{
				// update_velocity already includes the motility update
				//(*all_cells)[i]->phenotype.motility.update_motility_vector( (*all_cells)[i] ,(*all_cells)[i]->phenotype , time_since_last_mechanics );
				(*all_cells)[i]->functions.update_velocity( (*all_cells)[i], (*all_cells)[i]->phenotype, time_since_last_mechanics);
			}

			if( (*all_cells)[i]->functions.custom_cell_rule )
			{
				(*all_cells)[i]->functions.custom_cell_rule((*all_cells)[i], (*all_cells)[i]->phenotype, time_since_last_mechanics);
			}
		}
		// Calculate new positions
		#pragma omp parallel for
		for( int i=0; i < (*all_cells).size(); i++ )
		{
			if(!(*all_cells)[i]->is_out_of_domain && (*all_cells)[i]->is_movable)
			{
				(*all_cells)[i]->update_position(time_since_last_mechanics);
			}
		}

		// When somebody reviews this code, let's add proper braces for clarity!!!

		// Update cell indices in the container
		for( int i=0; i < (*all_cells).size(); i++ )
			if(!(*all_cells)[i]->is_out_of_domain && (*all_cells)[i]->is_movable)
				(*all_cells)[i]->update_voxel_in_container();
		last_mechanics_time=t;
	}

	initialzed=true;
	return;
}

/*----------------------------------------*/
/* Parallel version of the function above	*/
/*----------------------------------------*/

void Cell_Container::update_all_cells(double t, double phenotype_dt_ , double mechanics_dt_ , double diffusion_dt_ , mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	// secretions and uptakes. Syncing with BioFVM is automated.

	/*===============================================================================*/
	/* This is Section I to be parallelized, but I am parallelizing Section II first */
	/*===============================================================================*/

	/* Nothing to be parallelized in phenotype.secretion.advance(...) function below */
	#pragma omp parallel for
	for( int i=0; i < (*all_cells).size(); i++ )
	{
		(*all_cells)[i]->phenotype.secretion.advance( (*all_cells)[i], (*all_cells)[i]->phenotype , diffusion_dt_ );
	}

	//if it is the time for running cell cycle, do it!
	double time_since_last_cycle= t- last_cell_cycle_time;

	static double phenotype_dt_tolerance = 0.001 * phenotype_dt_;
	static double mechanics_dt_tolerance = 0.001 * mechanics_dt_;

	if( fabs(time_since_last_cycle-phenotype_dt_ ) < phenotype_dt_tolerance || !initialzed)
	{
		// Reset the max_radius in each voxel. It will be filled in set_total_volume
		// It might be better if we calculate it before mechanics each time
		std::fill(max_cell_interactive_distance_in_voxel.begin(), max_cell_interactive_distance_in_voxel.end(), 0.0);

		if(!initialzed)
		{
			time_since_last_cycle = phenotype_dt_;
		}

		// new as of 1.2.1 -- bundles cell phenotype parameter update, volume update, geometry update,
		// checking for death, and advancing the cell cycle. Not motility, though. (that's in mechanics)

		/* Nothing to be parallelized in advance_bundled_phenotype_functions(...) */
		#pragma omp parallel for
		for( int i=0; i < (*all_cells).size(); i++ )
		{
			if( (*all_cells)[i]->is_out_of_domain == false )
			{
				(*all_cells)[i]->advance_bundled_phenotype_functions( time_since_last_cycle );
			}
		}
		/*=================================================================================================*/
		/* Section II - division(), death(), compute_gradients(), update_veclocity() and update_position() */
		/*=================================================================================================*/

		// process divides / removes

		/*-------------------------------------------------------------------------------------------*/
		/* Each process must request cell_ready_to_divide.size() number of IDs from the root process */
		/*-------------------------------------------------------------------------------------------*/

		int no_of_requested_IDs = cells_ready_to_divide.size(); 		//no. of IDs requested by each process
		int range_ID[2]; 																						//range_ID[0] contains start ID, range_ID[1] contains end ID


		int strt_cell_ID = Basic_Agent::get_max_ID_in_parallel();		//Get the starting ID (the last agent ID plus one)

		int last_ID = request_IDs_from_root(no_of_requested_IDs, range_ID, strt_cell_ID, world, cart_topo);

		//Maybe I can just ask rank 0 to set the maxi
		Basic_Agent::set_max_ID_in_parallel(last_ID);								//This is the ID with which we start next

		int p_ID = range_ID[0];

		for( int i=0; i < cells_ready_to_divide.size(); i++ )
		{
				//cells_ready_to_divide[i]->divide();											//GS commented this out

				cells_ready_to_divide[i]->divide(p_ID, world, cart_topo); //Calls new parallel version of divide()
				p_ID = p_ID + 1; 																					//ID of cell to be created passed into divide()
		}

		for( int i=0; i < cells_ready_to_die.size(); i++ )
		{
			cells_ready_to_die[i]->die();
		}
		num_divisions_in_current_step+=  cells_ready_to_divide.size();
		num_deaths_in_current_step+=  cells_ready_to_die.size();

		cells_ready_to_die.clear();
		cells_ready_to_divide.clear();
		last_cell_cycle_time= t;
	}

	double time_since_last_mechanics= t- last_mechanics_time;

	// if( time_since_last_mechanics>= mechanics_dt || !initialzed)
	if( fabs(time_since_last_mechanics - mechanics_dt_) < mechanics_dt_tolerance || !initialzed)
	{
		if(!initialzed)
		{
			time_since_last_mechanics = mechanics_dt_;
		}

		// new February 2018
		// if we need gradients, compute them
		if( default_microenvironment_options.calculate_gradients )
		{
			microenvironment.compute_all_gradient_vectors();
		}
		// end of new in Feb 2018

		// Compute velocities

		pack_moore_info(world, cart_topo);
		#pragma omp parallel for
		for( int i=0; i < (*all_cells).size(); i++ )
		{

			/* For parallel version, use update_velocity_parallel instead of update_velocity */

			if(!(*all_cells)[i]->is_out_of_domain && (*all_cells)[i]->is_movable && (*all_cells)[i]->functions.update_velocity_parallel)
			{
				// update_velocity already includes the motility update
				//(*all_cells)[i]->phenotype.motility.update_motility_vector( (*all_cells)[i] ,(*all_cells)[i]->phenotype , time_since_last_mechanics );
				(*all_cells)[i]->functions.update_velocity_parallel( (*all_cells)[i], (*all_cells)[i]->phenotype, time_since_last_mechanics, world, cart_topo);
			}

			if( (*all_cells)[i]->functions.custom_cell_rule )
			{
				(*all_cells)[i]->functions.custom_cell_rule((*all_cells)[i], (*all_cells)[i]->phenotype, time_since_last_mechanics);
			}
		}

		// Calculate new positions
		#pragma omp parallel for
		for( int i=0; i < (*all_cells).size(); i++ )
		{
			if(!(*all_cells)[i]->is_out_of_domain && (*all_cells)[i]->is_movable)
			{
				//For serial call just remove the argument 'world'
				(*all_cells)[i]->update_position(time_since_last_mechanics, world, cart_topo);
			}
		}

		pack(all_cells, world, cart_topo);

		for(int i=0; i< (*all_cells).size(); i++)
    	if((*all_cells)[i]->crossed_to_left_subdomain || (*all_cells)[i]->crossed_to_right_subdomain)
    	{
    		//std::cout<<"CELL ID="<<(*all_cells)[i]->ID<<" removed from Rank="<<world.rank<<std::endl;

    		/* Since the current cell is exchanged with the last cell and the last cell is deleted */
    		/* a situation can arise where the current cell + last cell is to be deleted 					 */
    		/* In this case the last cell comes into the current position BUT since the loop 			 */
    		/* is incremented, the curent cell (i.e. the last cell) is NOT deleted 								 */
    		/* Thus, we need to decrement the counter IF we delete a cell 												 */

    		(*all_cells)[i]->remove_crossed_cell();
    		i = i - 1;
    	}

    /*=======================DELETE LATER===============================================*/
    //Checking for uniqueness of IDs - later delete

    // for(int temp_rank = 0; temp_rank < world.size; temp_rank++)
    // {
    // 	if(temp_rank == world.rank)
    // 	{
    // 		std::cout<<"IDs for Rank = "<<world.rank <<std::endl;
    // 		for(int i=0; i < (*all_cells).size(); i++)
    // 		{
    // 			for(int j=i+1; j < (*all_cells).size(); j++)
    // 			{
    // 				if((*all_cells)[i]->ID == (*all_cells)[j]->ID)
    // 					std::cout<<"Found a duplicate ID on Rank="<<world.rank<<std::endl;
    // 			}
    // 		}
    // 		std::cout<<std::endl;
    // 	}
    // 	MPI_Barrier(cart_topo.mpi_cart_comm);
    // }
    /*=======================DELETE LATER===============================================*/

    //Checking for uniqueness ends

		cart_topo.Send_and_Receive_Cells(no_cells_cross_left,  position_left,  snd_buf_left,
       															 no_cells_cross_right, position_right, snd_buf_right,
       															 no_of_cells_from_left,  rcv_buf_left,
       															 no_of_cells_from_right, rcv_buf_right,
       															 world
       															);

    unpack(world, cart_topo);

		/*------------------------------------------------------------------------*/
		/* Here I need to (1) Check if crossed_to_left_subdomain or crossed_to 		*/
		/* right_subdomain is true. If yes then call Cell packing function				*/
		/* send cell across and then run a loop on all cells again to delete			*/
		/* the cells whose crossed_to_left_subdomain or crossed_to_right_subdomain*/
		/* is true.																																*/
		/*------------------------------------------------------------------------*/

		// When somebody reviews this code, let's add proper braces for clarity!!!

		// Update cell indices in the container
		for( int i=0; i < (*all_cells).size(); i++ )
			if(!(*all_cells)[i]->is_out_of_domain && (*all_cells)[i]->is_movable)
				(*all_cells)[i]->update_voxel_in_container(world, cart_topo);					//Changed to parallel version
		last_mechanics_time=t;
	}

	initialzed=true;
	return;
}

void Cell_Container::pack_moore_info(mpi_Environment &world, mpi_Cartesian &cart_topo)
{


	/*--------------------------------------------------------------*/
	/* 4 variables below are data members of class Cell_Container 	*/
	/* These are used in sending packed cells also - use carefully	*/
	/*--------------------------------------------------------------*/

	position_left 			 		= 0;									//Must be initialized to 0
	position_right 			 		= 0;									//Must be initialized to 0
	snd_buf_left.resize(0); 											//When we enter function again, this is reset
	snd_buf_right.resize(0);											//Reset the others also

	int z_dim = underlying_mesh.z_coordinates.size();
	int y_dim = underlying_mesh.y_coordinates.size();

	int len_snd_buf_left  	 = 0;
	int len_snd_buf_right 	 = 0;

	/*---------------------*/
	/* Data Packing first  */
	/*---------------------*/

	/* Data Packing of the left boundary voxel cells of all processes except rank = 0 */

	if(world.rank > 0)
	{
		for(int k=0; k<z_dim; k++)
		{
			for(int j=0; j<y_dim; j++)
			{
				int local_vxl_inex 		 			= underlying_mesh.voxel_index(0, j, k);
				int global_mesh_index 			= underlying_mesh.voxels[local_vxl_inex].global_mesh_index;
				int no_of_cells_in_vxl 			= agent_grid[local_vxl_inex].size();
				std::vector<double> centers = underlying_mesh.voxels[local_vxl_inex].center;
				double max_i_dist						= max_cell_interactive_distance_in_voxel[local_vxl_inex];

				/* Pack the global voxel index, no. of cells in the voxel, 	 */
				/* voxel centers, and max_cell_interactive_distance_in_voxel */

				len_snd_buf_left = len_snd_buf_left + 2 * sizeof(int) + 4 * sizeof(double);
				snd_buf_left.resize(len_snd_buf_left);

				MPI_Pack(&global_mesh_index, 	1, MPI_INT, 	 &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
				MPI_Pack(&no_of_cells_in_vxl, 1, MPI_INT, 	 &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
				MPI_Pack(&centers[0], 				3, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
				MPI_Pack(&max_i_dist, 				1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);

				for(int vec_len=0; vec_len < no_of_cells_in_vxl; vec_len++)
				{
					Cell *pCell = agent_grid[local_vxl_inex][vec_len];

					len_snd_buf_left = len_snd_buf_left + 1 * sizeof(int) + 8 * sizeof(double);
					snd_buf_left.resize(len_snd_buf_left);

					MPI_Pack(&pCell->ID, 																							 			 1, MPI_INT, 		&snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->position[0], 																		 			 3, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.geometry.radius, 											 			 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.geometry.nuclear_radius, 							 			 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.mechanics.cell_cell_repulsion_strength, 			 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.mechanics.relative_maximum_adhesion_distance, 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.mechanics.cell_cell_adhesion_strength, 			 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
				}
			}
		}
	}

/* Data Packing of the right boundary voxel cells of all processes except rank = world.size-1 */

if(world.rank < world.size-1)
	{
		for(int k=0; k<z_dim; k++)
		{
			for(int j=0; j<y_dim; j++)
			{
				int local_vxl_inex 		 			= underlying_mesh.voxel_index(underlying_mesh.x_coordinates.size()-1, j, k);
				int global_mesh_index 			= underlying_mesh.voxels[local_vxl_inex].global_mesh_index;
				int no_of_cells_in_vxl 			= agent_grid[local_vxl_inex].size();
				std::vector<double> centers = underlying_mesh.voxels[local_vxl_inex].center;
				double max_i_dist						= max_cell_interactive_distance_in_voxel[local_vxl_inex];

				/* Pack the global voxel index, no. of cells in the voxel, 	 */
				/* voxel centers, and max_cell_interactive_distance_in_voxel */

				len_snd_buf_right = len_snd_buf_right + 2 * sizeof(int) + 4 * sizeof(double);
				snd_buf_right.resize(len_snd_buf_right);

				MPI_Pack(&global_mesh_index, 	1, MPI_INT, 	 &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
				MPI_Pack(&no_of_cells_in_vxl, 1, MPI_INT, 	 &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
				MPI_Pack(&centers[0], 				3, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
				MPI_Pack(&max_i_dist, 				1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

				for(int vec_len=0; vec_len < no_of_cells_in_vxl; vec_len++)
				{
					Cell *pCell = agent_grid[local_vxl_inex][vec_len];

					len_snd_buf_right = len_snd_buf_right + 1 * sizeof(int) + 8 * sizeof(double);
					snd_buf_right.resize(len_snd_buf_right);

					MPI_Pack(&pCell->ID, 																							 			 1, MPI_INT, 		&snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->position[0], 																		 			 3, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.geometry.radius, 											 			 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.geometry.nuclear_radius, 							 			 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.mechanics.cell_cell_repulsion_strength, 			 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.mechanics.relative_maximum_adhesion_distance, 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
					MPI_Pack(&pCell->phenotype.mechanics.cell_cell_adhesion_strength, 			 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
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

	MPI_Irecv(&rcv_buf_right[0], size_of_data_recvd_from_right_process , MPI_PACKED, cart_topo.X_RIGHT, 3333, cart_topo.mpi_cart_comm, &rcv_req[0]);
	MPI_Isend(&snd_buf_left[0], len_snd_buf_left, MPI_PACKED, cart_topo.X_LEFT,  3333, cart_topo.mpi_cart_comm, &snd_req[0]);

	MPI_Irecv(&rcv_buf_left[0], size_of_data_recvd_from_left_process, MPI_PACKED, cart_topo.X_LEFT, 4444, cart_topo.mpi_cart_comm, &rcv_req[1]);
	MPI_Isend(&snd_buf_right[0], len_snd_buf_right, MPI_PACKED, cart_topo.X_RIGHT, 4444, cart_topo.mpi_cart_comm, &snd_req[1]);

	MPI_Waitall(2, snd_req, MPI_STATUSES_IGNORE);
	MPI_Waitall(2, rcv_req, MPI_STATUSES_IGNORE);

	/*--------------------------------------------------------------------------*/
	/* Now tackle unpacking of cells in a special class called Moore_Voxel_Info */
	/*--------------------------------------------------------------------------*/

	/* Declare list of Moore_Voxel_Info - one element each for every boundary voxel i.e. y_dim * z_dim voxels */

	/* std::vector<Moore_Voxel_Info> mbfr (moore boundary from right [process]), similarly mbfl */
	/* is declared in PhysiCell_cell_container.h as a public data member 												*/

	mbfr.resize(y_dim * z_dim);
	mbfl.resize(y_dim * z_dim);

	/* Additionally, we have um_mbfl = unordered map moore boundary from left ---> 								 */
	/* maps the global mesh index of a voxel 																											 */
	/* to a specific Moore_Voxel_Info object. The idea is to take the global mesh index from 			 */
	/* moore_connected_voxel_global_indices_left/right and use this as a key to hash into the Moore*/
	/* Voxel_Info object associated with this global_mesh_index in O(1) time											 */
	/* No need to resize() um_mbfl and um_mbfr 																										 */


	position_right = 0;																			//Was used in Packing but re-used here, its a data member (remember)
	position_left = 0; 																			//Was used in Packing but r-used here, its a data member (remember)

	int size_right = size_of_data_recvd_from_right_process;	//For convenience
	int size_left =  size_of_data_recvd_from_left_process; 	//For convenience

	um_mbfr.clear();
	um_mbfl.clear();


	/* first unpack data coming FROM right process in 'rcv_buf_right' into 'mbfr' 					 */
	/* Outer loop runs till y_dim * z_dim voxels and each voxel can have any number of cells */
	/* IMPORTANT: rcv_buf_right will be unpacked ONLY on ranks < world.size-1 							 */
	/* REMEMBER: size_right = 0 on rank=world.size-1, hence error 													 */

	if(world.rank < world.size-1)
	{
		for(int vxl_ctr=0; vxl_ctr<y_dim*z_dim; vxl_ctr++)
		{
			MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &mbfr[vxl_ctr].global_mesh_index, 											1, MPI_INT, 		MPI_COMM_WORLD);
			MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &mbfr[vxl_ctr].no_of_cells_in_vxl, 										1, MPI_INT, 		MPI_COMM_WORLD);
			mbfr[vxl_ctr].center.resize(3);
			MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &mbfr[vxl_ctr].center[0], 															3, MPI_DOUBLE, 	MPI_COMM_WORLD);
			MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &mbfr[vxl_ctr].max_cell_interactive_distance_in_voxel, 1, MPI_DOUBLE, 	MPI_COMM_WORLD);

			/* Now resize vector moore_cells using no_of_cells_in_vxl unpacked above */

			int no_of_cells_in_vxl = mbfr[vxl_ctr].no_of_cells_in_vxl;
			mbfr[vxl_ctr].moore_cells.resize(no_of_cells_in_vxl);

			for(int cell_ctr=0; cell_ctr<no_of_cells_in_vxl; cell_ctr++)
			{
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &mbfr[vxl_ctr].moore_cells[cell_ctr].ID, 				 												1, MPI_INT, 	 MPI_COMM_WORLD);
				mbfr[vxl_ctr].moore_cells[cell_ctr].position.resize(3);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &mbfr[vxl_ctr].moore_cells[cell_ctr].position[0], 												3, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &mbfr[vxl_ctr].moore_cells[cell_ctr].radius, 		 												1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &mbfr[vxl_ctr].moore_cells[cell_ctr].nuclear_radius, 										1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &mbfr[vxl_ctr].moore_cells[cell_ctr].cell_cell_repulsion_strength, 			1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &mbfr[vxl_ctr].moore_cells[cell_ctr].relative_maximum_adhesion_distance, 1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_right[0], size_right, &position_right, &mbfr[vxl_ctr].moore_cells[cell_ctr].cell_cell_adhesion_strength, 				1, MPI_DOUBLE, MPI_COMM_WORLD);
			}

			  um_mbfr[mbfr[vxl_ctr].global_mesh_index]=mbfr[vxl_ctr];

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
			MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &mbfl[vxl_ctr].global_mesh_index, 											1, MPI_INT, 		MPI_COMM_WORLD);
			MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &mbfl[vxl_ctr].no_of_cells_in_vxl, 											1, MPI_INT, 		MPI_COMM_WORLD);
			mbfl[vxl_ctr].center.resize(3);
			MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &mbfl[vxl_ctr].center[0], 															3, MPI_DOUBLE, 	MPI_COMM_WORLD);
			MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &mbfl[vxl_ctr].max_cell_interactive_distance_in_voxel, 	1, MPI_DOUBLE, 	MPI_COMM_WORLD);

			/* Now resize vector moore_cells using no_of_cells_in_vxl unpacked above */

			int no_of_cells_in_vxl = mbfl[vxl_ctr].no_of_cells_in_vxl;
			mbfl[vxl_ctr].moore_cells.resize(no_of_cells_in_vxl);

			for(int cell_ctr=0; cell_ctr<no_of_cells_in_vxl; cell_ctr++)
			{
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &mbfl[vxl_ctr].moore_cells[cell_ctr].ID, 				 													1, MPI_INT, 	 MPI_COMM_WORLD);
				mbfl[vxl_ctr].moore_cells[cell_ctr].position.resize(3);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &mbfl[vxl_ctr].moore_cells[cell_ctr].position[0], 												3, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &mbfl[vxl_ctr].moore_cells[cell_ctr].radius, 		 													1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &mbfl[vxl_ctr].moore_cells[cell_ctr].nuclear_radius, 											1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &mbfl[vxl_ctr].moore_cells[cell_ctr].cell_cell_repulsion_strength, 				1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &mbfl[vxl_ctr].moore_cells[cell_ctr].relative_maximum_adhesion_distance, 	1, MPI_DOUBLE, MPI_COMM_WORLD);
				MPI_Unpack(&rcv_buf_left[0], size_left, &position_left, &mbfl[vxl_ctr].moore_cells[cell_ctr].cell_cell_adhesion_strength, 				1, MPI_DOUBLE, MPI_COMM_WORLD);
			}

			  um_mbfl[mbfl[vxl_ctr].global_mesh_index]=mbfl[vxl_ctr];
		}
	}

	/* Now one step is left: Need to map the global mesh index to a specific Moor_Voxel_Info object */
	/* i.e. declare map_glbl_indx_to_Moore_Voxel_Info_object[global_mesh_index]=mbfl[vxl_ctr] 			*/
	/* Similarly for mbfr as well 																																	*/
	/* Declare: std::unordered_map<int, Moore_Voxel_Info> in PhysiCell_cell_container.h 						*/


}

void Cell_Container::register_agent( Cell* agent )
{
	agent_grid[agent->get_current_mechanics_voxel_index()].push_back(agent);
	return;
}

void Cell_Container::remove_agent(Cell* agent )
{
	remove_agent_from_voxel(agent, agent->get_current_mechanics_voxel_index());
	return;
}

void Cell_Container::add_agent_to_outer_voxel(Cell* agent)
{
	int escaping_face= find_escaping_face_index(agent);
	agents_in_outer_voxels[escaping_face].push_back(agent);
	agent->is_out_of_domain=true;
	return;
}

void Cell_Container::remove_agent_from_voxel(Cell* agent, int voxel_index)
{
	int delete_index = 0;
	while( agent_grid[voxel_index][ delete_index ] != agent )
	{
		delete_index++;
	}
	// move last item to index location
	agent_grid[agent->get_current_mechanics_voxel_index()][delete_index] = agent_grid[agent->get_current_mechanics_voxel_index()][agent_grid[agent->get_current_mechanics_voxel_index()].size()-1 ];
	// shrink the vector
	agent_grid[agent->get_current_mechanics_voxel_index()].pop_back();
	return;
}

void Cell_Container::add_agent_to_voxel(Cell* agent, int voxel_index)
{
	agent_grid[voxel_index].push_back(agent);
	return;
}

bool Cell_Container::contain_any_cell(int voxel_index)
{
	// Let's replace this with clearer statements.
	return agent_grid[voxel_index].size()==0?false:true;
}

int find_escaping_face_index(Cell* agent)
{
	if(agent->position[0] <= agent->get_container()->underlying_mesh.bounding_box[PhysiCell_constants::mesh_min_x_index])
	{ return PhysiCell_constants::mesh_lx_face_index; }
	if(agent->position[0] >= agent->get_container()->underlying_mesh.bounding_box[PhysiCell_constants::mesh_max_x_index])
	{ return PhysiCell_constants::mesh_ux_face_index; }
	if(agent->position[1] <= agent->get_container()->underlying_mesh.bounding_box[PhysiCell_constants::mesh_min_y_index])
	{ return PhysiCell_constants::mesh_ly_face_index; }
	if(agent->position[1] >= agent->get_container()->underlying_mesh.bounding_box[PhysiCell_constants::mesh_max_y_index])
	{ return PhysiCell_constants::mesh_uy_face_index; }
	if(agent->position[2] <= agent->get_container()->underlying_mesh.bounding_box[PhysiCell_constants::mesh_min_z_index])
	{ return PhysiCell_constants::mesh_lz_face_index; }
	if(agent->position[2] >= agent->get_container()->underlying_mesh.bounding_box[PhysiCell_constants::mesh_max_z_index])
	{ return PhysiCell_constants::mesh_uz_face_index; }
	return -1;
}

void Cell_Container::flag_cell_for_division( Cell* pCell )
{
	#pragma omp critical
	{cells_ready_to_divide.push_back( pCell );}
	return;
}

void Cell_Container::flag_cell_for_removal( Cell* pCell )
{
	#pragma omp critical
	{cells_ready_to_die.push_back( pCell );}
	return;
}


Cell_Container* create_cell_container_for_microenvironment( BioFVM::Microenvironment& m , double mechanics_voxel_size )
{
	Cell_Container* cell_container = new Cell_Container;
	cell_container->initialize( m.mesh.bounding_box[0], m.mesh.bounding_box[3],
		m.mesh.bounding_box[1], m.mesh.bounding_box[4],
		m.mesh.bounding_box[2], m.mesh.bounding_box[5],  mechanics_voxel_size );
	m.agent_container = (Agent_Container*) cell_container;

	if( BioFVM::get_default_microenvironment() == NULL )
	{
		BioFVM::set_default_microenvironment( &m );
	}

	return cell_container;
}

/*------------------------------------------------------------------*/
/* Parallel version of create_cell_container_for_microenvironment() */
/*------------------------------------------------------------------*/

Cell_Container* create_cell_container_for_microenvironment( BioFVM::Microenvironment& m , double mechanics_voxel_size, mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	Cell_Container* cell_container = new Cell_Container;
	cell_container->initialize( m.mesh.bounding_box[0], m.mesh.bounding_box[3],
		m.mesh.bounding_box[1], m.mesh.bounding_box[4],
		m.mesh.bounding_box[2], m.mesh.bounding_box[5],  mechanics_voxel_size, world, cart_topo );
	m.agent_container = (Agent_Container*) cell_container;

	if( BioFVM::get_default_microenvironment() == NULL )
	{
		BioFVM::set_default_microenvironment( &m );
	}

	return cell_container;
}

};
