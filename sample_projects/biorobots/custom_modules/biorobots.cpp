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

#include "./biorobots.h"

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

void setup_microenvironment( mpi_Environment& world, mpi_Cartesian& cart_topo )
{
	// set domain parameters
	
	initialize_microenvironment(world, cart_topo); 	
	
	if (IOProcessor(world))
		microenvironment.display_information( std::cout ); 

	return; 
}

void create_cell_types( mpi_Environment& world, mpi_Cartesian& cart_topo)
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	//cell_definitions_by_index.clear();
	std::cout << " Cell definitions size in 1: " << cell_definitions_by_index.size() << std::endl;
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	std::cout << " Cell definitions size in 2: " << cell_definitions_by_index.size() << std::endl;

	initialize_cell_definitions_from_pugixml(world, cart_topo); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
	std::cout << " Cell definitions size in 3: " << cell_definitions_by_index.size() << std::endl;

	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(world, cart_topo); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 

	Cell_Definition* pCD = find_cell_definition( "director cell"); 
	pCD->functions.update_phenotype = director_cell_rule; 

	pCD = find_cell_definition( "cargo cell");
	pCD->functions.update_phenotype = cargo_cell_rule; 
	pCD->functions.contact_function = standard_elastic_contact_function; 

	pCD = find_cell_definition( "worker cell");
	pCD->functions.update_phenotype = worker_cell_rule; 
	pCD->functions.contact_function = standard_elastic_contact_function; 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
	std::cout << " Cell definitions size in 4: " << cell_definitions_by_index.size() << std::endl;

	display_cell_definitions( std::cout , world, cart_topo); 
	
	std::cout << " Cell definitions size in 5: " << cell_definitions_by_index.size() << std::endl;
	return; 
}

void director_cell_rule( Cell* pCell , Phenotype& phenotype , double dt )
{
	return; 
	std::vector<Cell*> nearby = pCell->cells_in_my_container(); 
	
	// if at least 2 neighbors, turn off secretion 
		// if size >= 3, then we have "self" and at least two more 
	if( nearby.size() > 2 )
	{
		pCell->phenotype.secretion.set_all_secretion_to_zero(); 
		pCell->custom_data[ "secreting" ] = 0.0; 
		
		pCell->functions.update_phenotype = NULL; 
	}
	
	return; 
}

std::vector<std::string> robot_coloring_function( Cell* pCell )
{
	std::string color = "black"; 
	std::vector< std::string > output( 4 , color ); 
	
	// black cells if necrotic 
	if( pCell->phenotype.death.dead == true )
	{ return output; }

	output[3] = "none"; // no nuclear outline color 
	
	static std::string worker_color = parameters.strings( "worker_color" ); 
	static std::string cargo_color = parameters.strings( "cargo_color" ); 
	static std::string director_color = parameters.strings( "director_color" ); 

	static int worker_ID = find_cell_definition( "worker cell" )->type; 
	static int cargo_ID = find_cell_definition( "cargo cell" )->type; 
	static int director_ID = find_cell_definition( "director cell" )->type; 

	if( pCell->type == worker_ID )
	{ color = worker_color; }
	else if( pCell->type == cargo_ID )
	{ color = cargo_color; }
	else if( pCell->type == director_ID )
	{ color = director_color; }
	
	output[0] = color; 
	output[2] = color; 
	
	return output; 
}

void create_cargo_cluster_6( std::vector<double>& center , mpi_Environment& world, mpi_Cartesian& cart_topo )
{
	// create a hollow cluster at position, with random orientation 
	static Cell_Definition* pCargoDef = find_cell_definition("cargo cell");	
	
	static double spacing = 0.95 * pCargoDef->phenotype.geometry.radius * 2.0; 
	static double d_Theta = 1.047197551196598 ; // 2*pi / 6.0 
	
	double theta = 6.283185307179586 * UniformRandom(); 
	
	static std::vector<double> position(3,0.0); 
	
	Cell* pC; 
	for( int i=0; i < 6; i++ )
	{
		pC = create_cell( *pCargoDef ); 
		
		position[0] = center[0] + spacing*cos( theta ); 
		position[1] = center[1] + spacing*sin( theta ); 
		
		pC->assign_position( position, world, cart_topo ); 
		
		theta += d_Theta; 
	}
	
	return; 
}

void create_cargo_cluster_7( std::vector<double>& center , mpi_Environment& world, mpi_Cartesian& cart_topo )
{
	// create a filled cluster at position, with random orientation 
	static Cell_Definition* pCargoDef = find_cell_definition("cargo cell");	

	create_cargo_cluster_6( center, world, cart_topo );
	Cell* pC = create_cell( *pCargoDef ); 
	pC->assign_position( center , world, cart_topo); 
	
	return; 
}


void create_cargo_cluster_3( std::vector<double>& center )
{
	// create a small cluster at position, with random orientation 
	static Cell_Definition* pCargoDef = find_cell_definition("cargo cell");	
	
	static double spacing = 0.95 * pCargoDef->phenotype.geometry.radius * 1.0; 
	static double d_Theta = 2.094395102393195 ; // 2*pi / 3.0 
	
	double theta = 6.283185307179586 * UniformRandom(); 
	
	static std::vector<double> position(3,0.0); 
	
	Cell* pC; 
	for( int i=0; i < 3; i++ )
	{
		pC = create_cell( *pCargoDef ); 
		
		position[0] = center[0] + spacing*cos( theta ); 
		position[1] = center[1] + spacing*sin( theta ); 
		
		pC->assign_position( position ); 
		
		theta += d_Theta; 
	}
	
	return; 
}


void setup_tissue( DistPhy::mpi::mpi_Environment &world, DistPhy::mpi::mpi_Cartesian &cart_topo )
{
	Microenvironment microenvironment = (*default_microenvironment_options.pMicroenvironment);

	double Xmin = microenvironment.mesh.x_coordinates[0]; 
	double Ymin = microenvironment.mesh.y_coordinates[0]; 
	double Zmin = microenvironment.mesh.z_coordinates[0];

	int x_voxels = microenvironment.mesh.x_coordinates.size();
	int y_voxels = microenvironment.mesh.y_coordinates.size();
	int z_voxels = microenvironment.mesh.z_coordinates.size();

	double Xmax = microenvironment.mesh.x_coordinates[x_voxels - 1]; 
	double Ymax = microenvironment.mesh.y_coordinates[y_voxels - 1]; 
	double Zmax = microenvironment.mesh.z_coordinates[z_voxels - 1];

	int number_of_directors = parameters.ints("number_of_directors"); // 15; 
	int local_directors =  (number_of_directors/world.size) + ((number_of_directors % world.size) > world.rank);
	int number_of_cargo_clusters = parameters.ints("number_of_cargo_clusters"); // 100;  
	int local_cargo_clusters = (number_of_cargo_clusters/world.size) + ((number_of_cargo_clusters % world.size) > world.rank);
	int number_of_workers = parameters.ints("number_of_workers"); // 50; 
	int local_workers =  (number_of_workers/world.size) + ((number_of_workers % world.size) > world.rank);

	if (IOProcessor(world))
		std::cout << "Placing cells ... " << std::endl; 
	
	// randomly place seed cells 
	
	std::vector<double> position(3,0.0); 
	
	double x_range = Xmax - Xmin; 
	double y_range = Ymax - Ymin;
	double z_range = Zmax - Zmin; 

	double relative_margin = 0.2;  
	double relative_outer_margin = 0.02; 
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell* pC;
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") / world.size ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*x_range; 
			position[1] = Ymin + UniformRandom()*y_range; 
			position[2] = Zmin + UniformRandom()*z_range; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position, world, cart_topo );
		}
	}
	std::cout << std::endl; 
	
	if (IOProcessor(world))
		std::cout << "\tPlacing " << number_of_directors << " director cells ... " << std::endl; 
	
	static Cell_Definition* pCargoDef = find_cell_definition("cargo cell");	
	static Cell_Definition* pDirectorDef = find_cell_definition("director cell");	
	static Cell_Definition* pWorkerDef = find_cell_definition("worker cell");

	for( int i=0; i < local_directors ; i++ )
	{
		// pick a random location 
		position[0] = Xmin + x_range*( relative_margin + (1.0-2*relative_margin)*UniformRandom() ); 
		
		position[1] = Ymin + y_range*( relative_outer_margin + (1.0-2*relative_outer_margin)*UniformRandom() ); 

		position[2] = Zmin + z_range*( relative_outer_margin + (1.0-2*relative_outer_margin)*UniformRandom() ); 
		
		// place the cell
		Cell* pC;
		pC = create_cell( *pDirectorDef ); 
		pC->assign_position( position, world, cart_topo );
		pC->is_movable = false; 
	}
	
	// place cargo clusters on the fringes 

	
	
	if (IOProcessor(world))
		std::cout << "\tPlacing cargo cells ... " << std::endl; 
	
	for( int i=0; i <  local_cargo_clusters ; i++ )
	{
		// pick a random location 
		
		position[0] = Xmin + x_range*( relative_margin + (1.0-2*relative_margin)*UniformRandom() ); 
		
		position[1] = Ymin + y_range*( relative_outer_margin + (1.0-2*relative_outer_margin)*UniformRandom() ); 

		position[2] = Zmin + z_range*( relative_outer_margin + (1.0-2*relative_outer_margin)*UniformRandom() );  
		
		if( UniformRandom() < 0.5 )
		{
			Cell* pCell = create_cell( *pCargoDef ); 
			pCell->assign_position( position  , world , cart_topo); 
		}
		else
		{
			create_cargo_cluster_7( position, world, cart_topo ); 
		}
	}
	
	// place "workersworkers"
	if (IOProcessor(world))
		std::cout << "\tPlacing worker cells ... " << std::endl; 
	for( int i=0; i < local_workers ; i++ )
	{
		// pick a random location 
		
		position[0] = Xmin + x_range*( relative_margin + (1.0-2*relative_margin)*UniformRandom() ); 
		
		position[1] = Ymin + y_range*( relative_outer_margin + (1.0-2*relative_outer_margin)*UniformRandom() ); 

		position[2] = Zmin + z_range*( relative_outer_margin + (1.0-2*relative_outer_margin)*UniformRandom() ); 
		
		// place the cell
		Cell* pC;

		pC = create_cell( *pWorkerDef ); 
		pC->assign_position( position, world, cart_topo);
	}	
	
	if (IOProcessor(world))
		std::cout << "done!" << std::endl; 
	// make a plot 
	
	PhysiCell_SVG_options.length_bar = 200; 
	//SVG_plot_mpi( "initial.svg" , microenvironment, 0.0 , 0.0 , robot_coloring_function, world, cart_topo);	
	
	return; 
}


void cargo_cell_rule( Cell* pCell , Phenotype& phenotype , double dt )
{
	
	return; 
}


void attach_cells( Cell* pCell_1, Cell* pCell_2 )
{
	#pragma omp critical
	{
		
	bool already_attached = false; 
	for( int i=0 ; i < pCell_1->state.neighbors.size() ; i++ )
	{
		if( pCell_1->state.neighbors[i] == pCell_2 )
		{ already_attached = true; }
	}
	if( already_attached == false )
	{ pCell_1->state.neighbors.push_back( pCell_2 ); }
	
	already_attached = false; 
	for( int i=0 ; i < pCell_2->state.neighbors.size() ; i++ )
	{
		if( pCell_2->state.neighbors[i] == pCell_1 )
		{ already_attached = true; }
	}
	if( already_attached == false )
	{ pCell_2->state.neighbors.push_back( pCell_1 ); }

	}

	return; 
}

void dettach_cells( Cell* pCell_1 , Cell* pCell_2 )
{
	#pragma omp critical
	{
		bool found = false; 
		int i = 0; 
		while( !found && i < pCell_1->state.neighbors.size() )
		{
			// if cell 2 is in cell 1's list, remove it
			if( pCell_1->state.neighbors[i] == pCell_2 )
			{
				int n = pCell_1->state.neighbors.size(); 
				// copy last entry to current position 
				pCell_1->state.neighbors[i] = pCell_1->state.neighbors[n-1]; 
				// shrink by one 
				pCell_1->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}
	
		found = false; 
		i = 0; 
		while( !found && i < pCell_2->state.neighbors.size() )
		{
			// if cell 1 is in cell 2's list, remove it
			if( pCell_2->state.neighbors[i] == pCell_1 )
			{
				int n = pCell_2->state.neighbors.size(); 
				// copy last entry to current position 
				pCell_2->state.neighbors[i] = pCell_2->state.neighbors[n-1]; 
				// shrink by one 
				pCell_2->state.neighbors.pop_back(); 
				found = true; 
			}
			i++; 
		}

	}
	
	return; 
}

void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant )
{
	std::vector<double> displacement = pAttachedTo->position - pActingOn->position; 
	axpy( &(pActingOn->velocity) , elastic_constant , displacement ); 
	
	return; 
}

void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if I am 
	std::vector<double> velocity(3,0.0); 
	
	for( int i=0; i < pCell->state.neighbors.size() ; i++ )
	{
		add_elastic_velocity( pCell, pCell->state.neighbors[i], pCell->custom_data["elastic coefficient"] ); 
	}

	return; 
}	


void worker_cell_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
	static double threshold = parameters.doubles("drop_threshold"); // 0.4; 
	
	static int cargo_index = microenvironment.find_density_index( "cargo signal" ); // 1 
	static int director_index = microenvironment.find_density_index( "director signal" ); // 0 
	
	// have I arrived? If so, release my cargo 
	if( pCell->nearest_density_vector()[director_index] > threshold )
	{
		for( int i=0; i < pCell->state.neighbors.size(); i++ )
		{
			Cell* pTemp = pCell->state.neighbors[i]; 
			dettach_cells( pCell, pTemp ); 
			
			pTemp->custom_data[ "receptor" ] = 0.0; 
			pTemp->phenotype.cycle.data.transition_rate( 0,0 ) = 0; 
		}
	}
	
	// am I searching for cargo? if so, see if I've found it
	if( pCell->state.neighbors.size() == 0 )
	{
		std::vector<Cell*> nearby = pCell->cells_in_my_container(); 
		for( int i=0; i < nearby.size(); i++ )
		{
			// if it is expressing the receptor, dock with it 
			if( nearby[i]->custom_data["receptor"] > 0.5 )
			{
				PhysiCell::attach_cells( pCell, nearby[i] ); 
				nearby[i]->custom_data["receptor"] = 0.0; 
				nearby[i]->phenotype.secretion.set_all_secretion_to_zero(); 
			}
		}
		
	}
	
	return; 
}

void worker_cell_motility( Cell* pCell, Phenotype& phenotype, double dt )
{
	// if attached, biased motility towards director chemoattractant 
	// otherwise, biased motility towards cargo chemoattractant 
	
	static double attached_worker_migration_bias = 
		parameters.doubles("attached_worker_migration_bias"); 
	static double unattached_worker_migration_bias = 
		parameters.doubles("unattached_worker_migration_bias"); 
		
	static int cargo_index = microenvironment.find_density_index( "cargo signal" ); // 1 
	static int director_index = microenvironment.find_density_index( "director signal" ); // 0 
	
	if( pCell->state.neighbors.size() > 0 )
	{
		phenotype.motility.migration_bias = attached_worker_migration_bias; 

		phenotype.motility.migration_bias_direction = pCell->nearest_gradient(director_index);	
		normalize( &( phenotype.motility.migration_bias_direction ) );			
	}
	else
	{
		phenotype.motility.migration_bias = unattached_worker_migration_bias; 
		
		phenotype.motility.migration_bias_direction = pCell->nearest_gradient(cargo_index);	
		normalize( &( phenotype.motility.migration_bias_direction ) );			
	}
	
	return; 
}
