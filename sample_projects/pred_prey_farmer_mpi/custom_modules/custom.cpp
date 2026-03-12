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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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

#include "./custom.h"
#include "../DistPhy/DistPhy_Utils.h"
#include "../DistPhy/DistPhy_Collective.h"
#include <mpi.h>
#include <vector>

using namespace DistPhy::mpi;

namespace
{
	void create_cell_types_base( mpi_Environment* world, mpi_Cartesian* cart_topo )
	{
		if (parameters.ints.find_index("random_seed") != -1)
		{
			SeedRandom(parameters.ints("random_seed"));
		}
		
		initialize_default_cell_definition(); 
		
		cell_defaults.functions.volume_update_function = standard_volume_update_function;
		cell_defaults.functions.update_velocity = standard_update_cell_velocity;
		cell_defaults.functions.update_velocity_parallel = standard_update_cell_velocity;

		cell_defaults.functions.update_migration_bias = weighted_motility_function; 
		cell_defaults.functions.update_phenotype = NULL; 
		cell_defaults.functions.custom_cell_rule = NULL; 
		
		cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
		cell_defaults.functions.calculate_distance_to_membrane = NULL; 
		
		if( world && cart_topo )
		{
			initialize_cell_definitions_from_pugixml( *world, *cart_topo ); 
		}
		else
		{
			initialize_cell_definitions_from_pugixml(); 
		}

		cell_defaults.functions.update_phenotype = phenotype_function; 
		cell_defaults.functions.custom_cell_rule = custom_function; 
		
		Cell_Definition* pFarmerDef = find_cell_definition( "farmer" ); 
		Cell_Definition* pPreyDef = find_cell_definition( "prey" ); 
		Cell_Definition* pPredDef = find_cell_definition( "predator" ); 

		pFarmerDef->functions.custom_cell_rule = farmer_custom_function; 
		pPreyDef->functions.custom_cell_rule = prey_custom_function; 
		pPreyDef->functions.update_phenotype = prey_phenotype_function; 
		pPreyDef->functions.update_migration_bias = prey_motility_function; 

		pPredDef->functions.custom_cell_rule = predator_custom_function; 
		pPredDef->functions.update_phenotype = predator_phenotype_function; 
		pPredDef->functions.update_migration_bias = predator_motility_function; 

		build_cell_definitions_maps(); 

		if( world && cart_topo )
		{
			setup_signal_behavior_dictionaries( *world, *cart_topo ); 	
		}
		else
		{
			setup_signal_behavior_dictionaries(); 	
		}

		setup_cell_rules(); 

		if( world == nullptr || IOProcessor(*world) )
		{
			display_cell_definitions( std::cout ); 
		}
	}

	struct TypeCount
	{
		std::string def_name;
		int count;
	};

	std::vector<double> random_position(double Xmin,double Xmax,double Ymin,double Ymax,double Zmin,double Zmax)
	{
		std::vector<double> position(3,0.0); 
		position[0] = Xmin + UniformRandom()* (Xmax - Xmin); 
		position[1] = Ymin + UniformRandom()* (Ymax - Ymin); 
		position[2] = Zmin + UniformRandom()* (Zmax - Zmin); 
		return position;
	}

	void distribute_type(const TypeCount& type_info, Microenvironment &m, mpi_Environment &world, mpi_Cartesian &cart_topo, double Xmin, double Xmax, double Ymin, double Ymax, double Zmin, double Zmax)
	{
		mpi_CellPositions cp;                
	    mpi_MyCells       mc;                

		if(world.rank == 0)
		{
			std::vector<std::vector<double>> positions_root;
			positions_root.reserve(type_info.count);
			for( int n = 0 ; n < type_info.count ; n++ )
			{
				positions_root.push_back( random_position(Xmin, Xmax, Ymin, Ymax, Zmin, Zmax) );
			}

			int start_cell_ID = Basic_Agent::get_max_ID_in_parallel();                               
	        
	        cp.positions_to_rank_list(positions_root, 
	                                  m.mesh.bounding_box[0], m.mesh.bounding_box[3], 
	                                  m.mesh.bounding_box[1], m.mesh.bounding_box[4], 
	                                  m.mesh.bounding_box[2], m.mesh.bounding_box[5], 
	                                  m.mesh.dx, m.mesh.dy, m.mesh.dz, 
	                                  world, cart_topo, start_cell_ID);
	        
	        Basic_Agent::set_max_ID_in_parallel(start_cell_ID + positions_root.size()); 
		}
		else
		{
			cp.no_of_IDs_all_procs.resize(world.size,0);
			cp.cell_IDs_all_procs.resize(world.size);
			cp.cell_coords_all_procs.resize(world.size);
		}

		distribute_cell_positions(cp, mc, world, cart_topo);                                        

		Cell_Definition* pDef = find_cell_definition( type_info.def_name );
		for( int i=0; i < mc.my_no_of_cell_IDs; i++ )
		{
			Cell* pC = create_cell( *pDef, mc.my_cell_IDs[i] ); 
			pC->assign_position( mc.my_cell_coords[3*i], mc.my_cell_coords[3*i+1], mc.my_cell_coords[3*i+2], world, cart_topo );
			for( int k=0 ; k < pC->phenotype.death.rates.size() ; k++ )
			{ pC->phenotype.death.rates[k] = 0.0; }
		}
	}
}

void create_cell_types( void )
{
	create_cell_types_base( nullptr, nullptr );
	return; 
}

void create_cell_types( mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	create_cell_types_base( &world, &cart_topo );
	return;
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

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
	
	// create some of each type of cell 
	
	Cell* pC;
	
	// place farmers 
	Cell_Definition* pCD = find_cell_definition( "farmer"); 
	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
	for( int n = 0 ; n < parameters.ints("number_of_farmers") ; n++ )
	{
		std::vector<double> position = {0,0,0}; 
		position[0] = Xmin + UniformRandom()*Xrange; 
		position[1] = Ymin + UniformRandom()*Yrange; 
		position[2] = Zmin + UniformRandom()*Zrange; 
		
		pC = create_cell( *pCD ); 
		pC->assign_position( position );
	}

	// place prey 
	pCD = find_cell_definition( "prey"); 
	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
	for( int n = 0 ; n < parameters.ints("number_of_prey") ; n++ )
	{
		std::vector<double> position = {0,0,0}; 
		position[0] = Xmin + UniformRandom()*Xrange; 
		position[1] = Ymin + UniformRandom()*Yrange; 
		position[2] = Zmin + UniformRandom()*Zrange; 
		
		pC = create_cell( *pCD ); 
		pC->assign_position( position );
	}

	// place predators 
	pCD = find_cell_definition( "predator"); 
	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
	for( int n = 0 ; n < parameters.ints("number_of_predators") ; n++ )
	{
		std::vector<double> position = {0,0,0}; 
		position[0] = Xmin + UniformRandom()*Xrange; 
		position[1] = Ymin + UniformRandom()*Yrange; 
		position[2] = Zmin + UniformRandom()*Zrange; 
		
		pC = create_cell( *pCD ); 
		pC->assign_position( position );
	}

/*	
	std::cout << std::endl; 
	
	fill_rectangle( {0,0,0,200,100,0} , cell_definitions_by_index[0] , 0.9 ); 
	
	fill_circle( {-143,-100,0} , 150 , 1 , 0.95 ); 
	
	draw_line( {-137,132,0} , {217,194,0} , 2 , 0.9 ); 
	
	fill_annulus( {300,300,0} , 150, 90 , 0 , 0.75 ); 
*/	
	// XML-based cell placement as needed 
	load_cells_from_pugixml(); 

	return; 
}

void setup_tissue( Microenvironment &m, mpi_Environment &world, mpi_Cartesian &cart_topo )
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

	std::vector<TypeCount> types = {
		{"farmer", parameters.ints("number_of_farmers")},
		{"prey", parameters.ints("number_of_prey")},
		{"predator", parameters.ints("number_of_predators")}
	};

	for( const auto& t : types )
	{
		distribute_type( t, m, world, cart_topo, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax );
	}

	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	static Cell_Definition* pFarmerDef = find_cell_definition( "farmer" ); 
	static Cell_Definition* pPreyDef = find_cell_definition( "prey" ); 
	static Cell_Definition* pPredDef = find_cell_definition( "predator" ); 
	
	if( pCell->type == pFarmerDef->type )
	{ return { "grey", "black", "grey", "grey" }; } 
	
	if( pCell->type == pPreyDef->type )
	{ return { "blue", "black", "blue", "blue" }; } 

	if( pCell->type == pPredDef->type )
	{ return { "orange", "black", "orange", "orange" }; } 

	return paint_by_number_cell_coloring(pCell); 
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

/* some help */ 

void weighted_motility_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// find the indices for each major substrate 
	static int prey_index = microenvironment.find_density_index( "prey signal"); 
	static int predator_index = microenvironment.find_density_index( "predator signal"); 
	static int food_index = microenvironment.find_density_index( "food"); 
	
	// zero out the motility bias direction. use a pointer to make this easier 
	std::vector<double>* pV = &phenotype.motility.migration_bias_direction; 
	(*pV) = {0,0,0}; // pCell->position; 
	//*pV *= -0.00001; 
 	
	// v += prey_weight * grad(prey) 
	axpy( pV, pCell->custom_data["prey_weight"] , pCell->nearest_gradient(prey_index) );
	// v += predator_weight * grad(predator) 
	axpy( pV, pCell->custom_data["predator_weight"] , pCell->nearest_gradient(predator_index) );
	// v += food_weight * grad(food) 
	axpy( pV, pCell->custom_data["food_weight"] , pCell->nearest_gradient(food_index) );
	
	normalize( pV ); 
}

void avoid_boundaries( Cell* pCell )
{
	// add velocity to steer clear of the boundaries 
	static double Xmin = microenvironment.mesh.bounding_box[0]; 
	static double Ymin = microenvironment.mesh.bounding_box[1]; 
	static double Zmin = microenvironment.mesh.bounding_box[2]; 

	static double Xmax = microenvironment.mesh.bounding_box[3]; 
	static double Ymax = microenvironment.mesh.bounding_box[4]; 
	static double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	static double avoid_zone = 25; 
	static double avoid_speed = -0.5; // must be negative 
	
	// near edge: 
	bool near_edge = false; 
	if( pCell->position[0] < Xmin + avoid_zone || pCell->position[0] > Xmax - avoid_zone )
	{ near_edge = true; } 
	
	if( pCell->position[1] < Ymin + avoid_zone || pCell->position[1] > Ymax - avoid_zone )
	{ near_edge = true; } 
	
	if( default_microenvironment_options.simulate_2D == false )
	{
		if( pCell->position[2] < Zmin + avoid_zone || pCell->position[2] > Zmax - avoid_zone )
		{ near_edge = true; } 
	}
	
	if( near_edge )
	{
		pCell->velocity = pCell->position; // move towards origin 
		pCell->velocity *= avoid_speed; // move towards origin 
	}
	
	return; 
}

void wrap_boundaries( Cell* pCell )
{
	return avoid_boundaries( pCell ); 
	
	// add velocity to steer clear of the boundaries 
	static double Xmin = microenvironment.mesh.bounding_box[0]; 
	static double Ymin = microenvironment.mesh.bounding_box[1]; 
	static double Zmin = microenvironment.mesh.bounding_box[2]; 

	static double Xmax = microenvironment.mesh.bounding_box[3]; 
	static double Ymax = microenvironment.mesh.bounding_box[4]; 
	static double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	static double avoid_zone = 20; 

	static bool setup_done = false; 
	if( setup_done == false )
	{
		Xmax -= avoid_zone; 
		Xmin += avoid_zone; 
		
		Ymax -= avoid_zone;
		Ymin += avoid_zone; 
		
		setup_done = true; 
	}
	
	bool wrapped = false; 
	
	std::vector<double> p = pCell->position;
	double Delta;


	while( p[0] < Xmin )
	{
		Delta = Xmin - p[0]; 
		p[0] = Xmax - Delta; 
		wrapped = true; 
	}
	while( p[0] > Xmax )
	{
		Delta = p[0] - Xmax; 
		p[0] = Xmin + Delta; 
		wrapped = true; 
	}
	
	while( p[1] < Ymin )
	{
		Delta = Ymin - p[1]; 
		p[1] = Ymax - Delta; 
		wrapped = true; 
	}
	while( p[1] > Ymax )
	{
		Delta = p[1] - Ymax; 
		p[1] = Ymin + Delta; 
		wrapped = true; 
	}

	if( default_microenvironment_options.simulate_2D == false )
	{
		while( p[2] < Zmin )
		{
			Delta = Zmin - p[2]; 
			p[2] = Zmax - Delta; 
			wrapped = true; 
		}
		while( p[2] > Zmax )
		{
			Delta = p[2] - Zmax; 
			p[2] = Zmin + Delta; 
			wrapped = true; 
		}
	}

	if( wrapped == true ) 
	{ 
		#pragma omp critical 
		{ pCell->assign_position( p ); }
	} 
	return; 
}

std::vector<Cell*> get_possible_neighbors( Cell* pCell)
{
	std::vector<Cell*> neighbors = {};
	
	// First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end =
		pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
		
	for( neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{ neighbors.push_back( *neighbor ); }
	
	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end
		= pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();
	
	for( neighbor_voxel_index = pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin(); neighbor_voxel_index!= neighbor_voxel_index_end; ++neighbor_voxel_index)
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
		continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{ neighbors.push_back( *neighbor ); }
	}
	return neighbors;
}	

/* prey functions */ 

void prey_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// sample food
	static int nFood = microenvironment.find_density_index("food"); 
	double food = pCell->nearest_density_vector()[nFood]; 
	
	// death based on food
	static int nNecrosis = phenotype.death.find_death_model_index( "necrosis" );
	
	if( food < 0.1 )
	{
		pCell->start_death( nNecrosis ); 
		pCell->functions.update_phenotype = NULL; 
		return; 
	} 
	
	// division based on food
	static Cell_Definition* pCD = find_cell_definition( "prey" ); 
	phenotype.cycle.data.exit_rate(0) = pCD->phenotype.cycle.data.exit_rate(0); 
	double multiplier = (food-0.1)/0.9; 
	phenotype.cycle.data.exit_rate(0) *= multiplier; 
	
	return; 
}

void prey_custom_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	wrap_boundaries( pCell );
}

void prey_motility_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	return weighted_motility_function(pCell, phenotype, dt ); 
}

/* predator functions */ 

void predator_phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	static Cell_Definition* pFarmerDef = find_cell_definition( "farmer" ); 
	static Cell_Definition* pPreyDef = find_cell_definition( "prey" ); 
	static Cell_Definition* pPredDef = find_cell_definition( "predator" ); 	
	
	// hunting 

	static double max_detection_distance = 2; 
	
	// see who is nearby 
	
	std::vector<Cell*> nearby = get_possible_neighbors( pCell); 
	
	for( int i=0 ; i < nearby.size() ; i++ )
	{
		Cell* pC = nearby[i]; 
		// is it prey ? 
		
		if( pC->type == pPreyDef->type )
		{
			bool eat_it = true; 
			// in range? 
			std::vector<double> displacement = pC->position; 
			displacement -= pCell->position; 
			double distance = norm( displacement ); 
			if( distance > pCell->phenotype.geometry.radius + pC->phenotype.geometry.radius 
				+ max_detection_distance )
			{ eat_it = false; } 
			
			// am I hungry? 
			
			if( eat_it == true )
			{
				// eat it! 
				pCell->ingest_cell( pC ); 
				
				// increase energy 
				pCell->custom_data["energy"] += 100; 	
			}
		}
	}
	
	// update energy 
	
	static double decay_rate = 0.00025; 
	pCell->custom_data["energy"] /= (1.0 + dt*decay_rate); 
	
	// low energy kills
	
	// death based on food
	static int nNecrosis = phenotype.death.find_death_model_index( "necrosis" );
	
	if( pCell->custom_data["energy"] < 0.1 )
	{
		pCell->start_death( nNecrosis ); 
		pCell->functions.update_phenotype = NULL; 
		return; 
	} 	
	
	
	// need energy to reproduce 
	
	return; 
}

void predator_custom_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	static Cell_Definition* pFarmerDef = find_cell_definition( "farmer" ); 
	static Cell_Definition* pPreyDef = find_cell_definition( "prey" ); 
	static Cell_Definition* pPredDef = find_cell_definition( "predator" ); 	
	
	wrap_boundaries( pCell ); 
	
	return; 
}

void predator_motility_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	return weighted_motility_function(pCell,phenotype,dt); 	
}

/* farmer functions */ 

void farmer_custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{
	return wrap_boundaries( pCell ); 
} 
