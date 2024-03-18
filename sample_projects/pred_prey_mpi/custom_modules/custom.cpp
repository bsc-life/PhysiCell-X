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

/*================================================================================
+ If you use PhysiCell-X in your project, we would really appreciate if you can  +
+																																							   +
+ [1] Cite the PhysiCell-X repository by giving its URL												   +
+																																							   +
+ [2] Cite BioFVM-X: 																													   +
+		Saxena, Gaurav, Miguel Ponce-de-Leon, Arnau Montagud, David Vicente Dorca,   +
+		and Alfonso Valencia. "BioFVM-X: An MPI+ OpenMP 3-D Simulator for Biological + 
+		Systems." In International Conference on Computational Methods in Systems    +
+		Biology, pp. 266-279. Springer, Cham, 2021. 																 +
=================================================================================*/

#include "./custom.h"
#include "../DistPhy/DistPhy_Utils.h"
#include "../DistPhy/DistPhy_Collective.h"

using namespace DistPhy::mpi;

// Declare cell definitions here 

Cell_Definition motile_cell; 

void create_cell_types( void )
{
	// Set the random seed 
	
	SeedRandom( parameters.ints("random_seed") );  
	
	/*----------------------------------------------------------------*/ 
	/*   Put any modifications to default cell definition here if you */
	/*   want to have "inherited" by other cell types. 								*/
	/*   This is a good place to set default functions. 							*/
	/*----------------------------------------------------------------*/ 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	
	/*-----------------------------------------------------------------------------------------------------------------------------*/
	/* For parallel settings, set the update velocity in parallel function pointer - this will call the new function which has the */
	/* same name but a different prototype and implementation (overloaded function)																								 */
	/* There may be no need to set it at this stage but setting it here for correctness and completeness. 												 */
	/*-----------------------------------------------------------------------------------------------------------------------------*/
	
	cell_defaults.functions.update_velocity_parallel = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; 			// update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	// The following parses the cell definitions in the XML config file
	
	initialize_cell_definitions_from_pugixml(); 
	
	/*-----------------------------------------------------------*/ 
	/* Put any modifications to individual cell definitions here.*/ 
	/*  This is a good place to set custom functions. 					 */
	/*-----------------------------------------------------------*/ 
	
	if( parameters.bools("predators_eat_prey") == true )
	{ 
		get_cell_definition("predator").functions.custom_cell_rule = predator_hunting_function; 
	}

	if( parameters.bools("predators_cycle_if_big") == true )
	{ 
		get_cell_definition("predator").functions.update_phenotype = predator_cycling_function; 
	}

	if( parameters.bools("prey_quorom_effect") == true )
	{ 
		get_cell_definition("prey").functions.update_phenotype = prey_cycling_function; 
	}
		
	// This builds the map of cell definitions and summarizes the setup. 
		
	build_cell_definitions_maps(); 
	
	//display_cell_definitions( std::cout ); <---- Will be printed out by all processes 
	
	/*---------------------------------------------------------------------------------------*/
	/* display_cell_definitions(...) has been disabled above as the printing will be done by */
	/* ALL processes. If you want to print (for checking, testing etc.) then create a new 	 */
	/* version of void create_cell_types(mpi_Environment &world, mpi_Cartesian &cart_topo)	 */
	/* and then use the IOProcessor(world) function to print using ONLY the root process  	 */
	/*---------------------------------------------------------------------------------------*/
	
	return; 
}

void setup_microenvironment( void )
{
	// Set domain parameters, put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here, initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

/*==============================================*/
/* Parallel version of setup_microenvironment() */
/*==============================================*/

void setup_microenvironment( mpi_Environment &world, mpi_Cartesian &cart_topo )
{	
	initialize_microenvironment(world, cart_topo); 	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// Create some of each type of cells 
	
	Cell* pC;
	
	// Place the preys 
	
	for( int n = 0 ; n < parameters.ints("number_of_prey") ; n++ )
	{
		std::vector<double> position = {0,0,0}; 
		position[0] = Xmin + UniformRandom()*Xrange; 
		position[1] = Ymin + UniformRandom()*Yrange; 
		position[2] = Zmin + UniformRandom()*Zrange; 
		
		pC = create_cell( get_cell_definition("prey") ); 
		pC->assign_position( position );
	}
	
	// Place the predators 
	
	for( int n = 0 ; n < parameters.ints("number_of_predators") ; n++ )
	{
		std::vector<double> position = {0,0,0}; 
		position[0] = Xmin + UniformRandom()*Xrange; 
		position[1] = Ymin + UniformRandom()*Yrange; 
		position[2] = Zmin + UniformRandom()*Zrange; 
		
		pC = create_cell( get_cell_definition("predator") ); 
		pC->assign_position( position );
	}	
	
	return; 
}

/*------------------------------------*/
/* Parallel version of setup_tissue() */
/*------------------------------------*/

void setup_tissue(Microenvironment &m, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	double Xmin = m.mesh.bounding_box[0]; 
	double Ymin = m.mesh.bounding_box[1]; 
	double Zmin = m.mesh.bounding_box[2]; 

	double Xmax = m.mesh.bounding_box[3]; 
	double Ymax = m.mesh.bounding_box[4]; 
	double Zmax = m.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// Create some of each type of cell 
	
	// The following 3 temporary variables are common to both Prey and Predator sections 

	Cell* pCell;
	std::vector<std::vector<double>> generated_positions_at_root;
	std::vector<double> position = {0,0,0};														//Temporary buffer
	
	/*------------------*/
	/* Prey Section 		*/
	/*------------------*/	

  mpi_CellPositions cp_prey;                		//To store cell positions, cell IDs, no. of cell IDs at root only (for all processes)
  mpi_MyCells       mc_prey;                		//To store cell positions, cell IDs, no. of cells at each process.
  

 // First Generate Prey positions on rank 0 
 	
 if(world.rank == 0)
 {
	for( int n = 0 ; n < parameters.ints("number_of_prey") ; n++ )
	{
		 
		position[0] = Xmin + UniformRandom()*Xrange; 
		position[1] = Ymin + UniformRandom()*Yrange; 
		position[2] = Zmin + UniformRandom()*Zrange;
		generated_positions_at_root.push_back(position); 
	}
	
	/*---------------------------------------------------------------------------------------------------*/
	/* (1) Obtain the current highest ID - rank 0 knows this as rank 0 generates all IDs 								 */
	/* (2) Distribute the positions and IDs to respective processes 										 								 */
	/* (3) Set the maximum ID as the current maximum ID would be needed if new cells are to be generated */
	/*---------------------------------------------------------------------------------------------------*/
	
	int strt_cell_ID = Basic_Agent::get_max_ID_in_parallel();                               //IDs for new cells (positions) will start from the current highest ID      
  cp_prey.positions_to_rank_list(generated_positions_at_root, 
                            		 m.mesh.bounding_box[0], m.mesh.bounding_box[3], 
                            		 m.mesh.bounding_box[1], m.mesh.bounding_box[4], 
                            		 m.mesh.bounding_box[2], m.mesh.bounding_box[5], 
                            		 m.mesh.dx, m.mesh.dy, m.mesh.dz, 
                            		 world, cart_topo, strt_cell_ID);
        
  Basic_Agent::set_max_ID_in_parallel(strt_cell_ID + generated_positions_at_root.size()); //Highest ID now is the starting ID + no. of generated coordinates !
 }
 
 distribute_cell_positions(cp_prey, mc_prey, world, cart_topo);                           //Distribute cell positions
 
 // Create preys at individual processes 
 
 for( int i=0; i < mc_prey.my_no_of_cell_IDs; i++ )
	{	
 		pCell = create_cell( get_cell_definition("prey"), mc_prey.my_cell_IDs[i] ); 
		pCell->assign_position(mc_prey.my_cell_coords[3*i],mc_prey.my_cell_coords[3*i+1],mc_prey.my_cell_coords[3*i+2],world, cart_topo); //pCell->assign_position( positions[i] );

	}

/*-------------------------------------------------------------------------------------------------------------------*/	
/* Clear off the generated_positions_at_root vector of vectors as now we need it for storing positions of Predators. */
/* Exactly the same procedure as above is repeated for Predators now 																								 */
/*-------------------------------------------------------------------------------------------------------------------*/	

	generated_positions_at_root.clear();

	/*------------------*/
	/* Predator Section */
	/*------------------*/	

  mpi_CellPositions cp_pred;                		//To store cell positions, cell IDs, no. of cell IDs at root only (for all processes)
  mpi_MyCells       mc_pred;                		//To store cell positions, cell IDs, no. of cells at each process. 
	
	if(world.rank == 0)
 	{
		for( int n = 0 ; n < parameters.ints("number_of_predators") ; n++ )
		{ 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange;
			generated_positions_at_root.push_back(position); 
		}
	
		int strt_cell_ID = Basic_Agent::get_max_ID_in_parallel();                               //IDs for new cells (positions) will start from the current highest ID      
  	cp_pred.positions_to_rank_list(generated_positions_at_root, 
                            			 m.mesh.bounding_box[0], m.mesh.bounding_box[3], 
                            			 m.mesh.bounding_box[1], m.mesh.bounding_box[4], 
                            			 m.mesh.bounding_box[2], m.mesh.bounding_box[5], 
                            			 m.mesh.dx, m.mesh.dy, m.mesh.dz, 
                            			 world, cart_topo, strt_cell_ID);
        
  	Basic_Agent::set_max_ID_in_parallel(strt_cell_ID + generated_positions_at_root.size()); //Highest ID now is the starting ID + no. of generated coordinates !
 	}
 
 distribute_cell_positions(cp_pred, mc_pred, world, cart_topo);                             //Distribute cell positions
 
 for( int i=0; i < mc_pred.my_no_of_cell_IDs; i++ )
	{	
 		pCell = create_cell( get_cell_definition("predator"), mc_pred.my_cell_IDs[i] ); 
		pCell->assign_position(mc_pred.my_cell_coords[3*i],mc_pred.my_cell_coords[3*i+1],mc_pred.my_cell_coords[3*i+2],world, cart_topo); //pCell->assign_position( positions[i] );
	}
	
	/*---------------------------------------------------------------------------------*/	
	/* Please note that the Predator cells and Prey cells will be displayed separately */
	/* in the output. For example if the initial number of preys is 20 and the no. of  */
	/* predators is 5, then you see an output like: 																	 */
	//	 	MPI Rank = 1 No of cells = 13 <---- Preys on MPI Rank 1
	//		MPI Rank = 0 No of cells = 7  <---- Preys on MPI Rank 0
	//		MPI Rank = 0 No of cells = 2  <---- Predators on Rank 0
	//		MPI Rank = 1 No of cells = 3  <---- Predators on Rank 1
	/*---------------------------------------------------------------------------------*/
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	static int prey_type = get_cell_definition( "prey" ).type; 
	static int predator_type = get_cell_definition( "predator" ).type; 
	
	// Start with flow cytometry coloring 
	
	std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 
	
	// Color live prey 
		
	if( pCell->phenotype.death.dead == false && pCell->type == prey_type )
	{
		 output[0] = parameters.strings("prey_color");  
		 output[2] = parameters.strings("prey_color");  
	}
	
	// Color live predators 

	if( pCell->phenotype.death.dead == false && pCell->type == predator_type )
	{
		 output[0] = parameters.strings("predator_color");  
		 output[2] = parameters.strings("predator_color");  
	}
	
	return output; 
}


void predator_hunting_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	Cell* pTestCell = NULL; 
	
	double sated_volume = pCell->parameters.pReference_live_phenotype->volume.total * 
		parameters.doubles("relative_sated_volume" ); 
	
	for( int n=0; n < pCell->cells_in_my_container().size() ; n++ )
	{
		pTestCell = pCell->cells_in_my_container()[n]; 
		
		// If it is not me, not dead, and not my type, eat it 
		
		if( pTestCell != pCell && pTestCell->type != pCell->type && pTestCell->phenotype.death.dead == false )
		{
			// Only eat if I am not full 
			
			if( phenotype.volume.total < sated_volume )
			{
				pCell->ingest_cell(pTestCell); 
				return; 
			}
	
		}
	}
	
	return; 
}

void predator_cycling_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	double sated_volume = pCell->parameters.pReference_live_phenotype->volume.total * 
		parameters.doubles("relative_sated_volume" ); 
	
	if( phenotype.volume.total > sated_volume )
	{ 
		phenotype.cycle.data.transition_rate(0,1) = get_cell_definition("prey").phenotype.cycle.data.transition_rate(0,1) * 0.01; 
	}
	else
	{ 
		phenotype.cycle.data.transition_rate(0,1) = 0; 
	}
	return; 
}

void prey_cycling_function( Cell* pCell , Phenotype& phenotype, double dt )
{
	static int signal_index = microenvironment.find_density_index( "prey signal" ); 
	
	double threshold = parameters.doubles("prey_quorom_threshold" ) + 1e-16 ; 
	double factor = (threshold - pCell->nearest_density_vector()[signal_index] )/threshold; 
	if( factor < 0 )
	{ 
		factor = 0.0; 
	} 
	
	phenotype.cycle.data.transition_rate(0,1) = get_cell_definition("prey").phenotype.cycle.data.transition_rate(0,1); 
	phenotype.cycle.data.transition_rate(0,1) *= factor; 
	
	return; 
}
int total_basic_agent_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	int local_count = all_basic_agents.size();
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}



int total_cell_agent_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	int local_count = (*all_cells).size();
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}


// Cell count functions
int total_live_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	int local_count = 0;
	#pragma omp parallel for reduction (+: local_count )
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell* pCell = (*all_cells)[i];
		if(pCell->phenotype.death.dead==false)
		{ local_count++;; }
	}
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}

int total_dead_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
    int local_count = 0;
    #pragma omp parallel for reduction (+: local_count )
    for (int i = 0; i < (*all_cells).size(); i++)
    {
        Cell* pCell = (*all_cells)[i];
        if(pCell->phenotype.death.dead==true)
        { local_count++; }
    }
    int global_count;
    MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
    return global_count;
}

int total_necrosis_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	int local_count = 0;
	#pragma omp parallel for reduction (+: local_count )
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell* pCell = (*all_cells)[i];
		if( pCell->phenotype.death.dead==false )
		{ continue; }
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )	
		{ local_count++; }
	}
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}

int total_apoptosis_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{	
	int local_count = 0;
	#pragma omp parallel for reduction (+: local_count )
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell* pCell = (*all_cells)[i];
		if( pCell->phenotype.death.dead==false )
		{ continue; }
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )
		{ local_count++; }
	}
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}
