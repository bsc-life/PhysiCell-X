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
# [1] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
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
		cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
		
		cell_defaults.functions.volume_update_function = standard_volume_update_function;
		cell_defaults.functions.update_velocity = standard_update_cell_velocity;
		cell_defaults.functions.update_velocity_parallel = standard_update_cell_velocity;

		cell_defaults.functions.update_migration_bias = NULL; 
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
		
		cell_defaults.functions.update_phenotype = NULL; 
		cell_defaults.functions.custom_cell_rule = custom_function; 
		cell_defaults.functions.contact_function = contact_function; 
		
		cell_defaults.phenotype.mechanics.attachment_elastic_constant = 
			parameters.doubles("attachment_elastic_constant"); 
			
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

		if( world == nullptr || IOProcessor( *world ) )
		{
			display_cell_definitions( std::cout ); 
		}
	}

	void seed_definition(Cell_Definition* pCD, int count, Microenvironment &m, mpi_Environment &world, mpi_Cartesian &cart_topo)
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

		std::vector<std::vector<double>> positions_root;

		mpi_CellPositions cp;                
	    mpi_MyCells       mc;                

		if(world.rank == 0)
		{
			positions_root.reserve( count );
			for( int n = 0 ; n < count ; n++ )
			{
				std::vector<double> position(3,0.0); 
				position[0] = Xmin + UniformRandom() * (Xmax - Xmin); 
				position[1] = Ymin + UniformRandom() * (Ymax - Ymin); 
				position[2] = Zmin + UniformRandom() * (Zmax - Zmin); 
				positions_root.push_back( position );
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

		for( int i=0; i < mc.my_no_of_cell_IDs; i++ )
		{
			Cell* pC = create_cell( *pCD, mc.my_cell_IDs[i] ); 
			pC->assign_position( mc.my_cell_coords[3*i], mc.my_cell_coords[3*i+1], mc.my_cell_coords[3*i+2], world, cart_topo );

			pC->custom_data["head"] = UniformRandom(); 
			pC->custom_data["head_initial"] = pC->custom_data["head"];
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
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	load_cells_from_pugixml(); 	
	
	for( int n=0; n < (*all_cells).size(); n++ )
	{
		Cell* pCell = (*all_cells)[n]; 
		pCell->custom_data["head"] = UniformRandom(); 
		pCell->custom_data["head_initial"] = pCell->custom_data["head"];
	}
	
	return; 
}

void setup_tissue( Microenvironment &m, mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		if( IOProcessor(world) )
		{
			std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		}
		seed_definition( pCD, parameters.ints("number_of_cells"), m, world, cart_topo );
	}

	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	if( pCell->state.number_of_attached_cells() == 0 )
	{ return { "grey", "black", "grey", "grey"}; }

	if( pCell->state.number_of_attached_cells() == 1 )
	{
		if( pCell->custom_data["head"] > pCell->state.attached_cells[0]->custom_data["head"] )
		{ return { "red", "black", "red", "red"};  } 
		
		return { "orange", "black", "orange", "orange"}; 
	}

	if( pCell->state.number_of_attached_cells() >= 2 )
	{
		int intensity = (int) floor( 255.0 * pCell->custom_data["head"] ); 
		std::string strColor = std::to_string(intensity); 
		std::string color = "rgb(" + strColor + "," + strColor + ",255)"; 
		
		if( pCell->state.number_of_attached_cells() > 2 )
		{ return { "yellow", "black" , color, color }; }

		return { color , "black", color, color}; 	
	}

	return { "yellow", "black", "yellow", "yellow" }; 
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{
	static int nSignal = microenvironment.find_density_index("signal");
	
	int number_of_attachments = pCell->state.number_of_attached_cells(); 
	std::vector<Cell*> nearby = pCell->nearby_interacting_cells(); 
	
	if( number_of_attachments == 0 )
	{
		int n = 0; 
		while( number_of_attachments < (int) pCell->custom_data["max_attachments"] && n < nearby.size() )
		{
			if( nearby[n]->state.number_of_attached_cells() < nearby[n]->custom_data["max_attachments"] )
			{
				attach_cells( nearby[n] , pCell ); 
				number_of_attachments++;
			}
			n++; 
		}
	}

	if( number_of_attachments == 0 )
	{ pCell->functions.update_migration_bias = chemotaxis_function; } 
	
	if( number_of_attachments == 1 )
	{
		pCell->custom_data["head"] = pCell->custom_data["head_initial"];

		bool head = false; 
		if( pCell->custom_data["head"] > pCell->state.attached_cells[0]->custom_data["head"] )
		{ head = true; } 
		
		if( head )
		{ pCell->functions.update_migration_bias = head_migration_direction; }
		else
		{ pCell->functions.update_migration_bias = tail_migration_direction; }
		phenotype.secretion.secretion_rates[nSignal] = 100; 
	} 
	
	if( number_of_attachments > 1 )
	{
		pCell->functions.update_migration_bias = middle_migration_direction;
		phenotype.secretion.secretion_rates[nSignal] = 1; 
	} 
	
	return; 
} 

void contact_function( Cell* pMe, Phenotype& phenoMe, 
	Cell* pOther, Phenotype& phenoOther, double dt )
{
	standard_elastic_contact_function(pMe,phenoMe,pOther,phenoOther,dt);

	if( pMe->state.number_of_attached_cells() > 0 )
	{
		double head_me = pMe->custom_data["head"];
		double head_other = pOther->custom_data["head"]; 
		
		if( head_me > head_other )
		{
			double amount_to_transfer = dt * pMe->custom_data["transfer_rate"] 
				* (head_me - head_other ); 
			pMe->custom_data["head"] -= amount_to_transfer; 
			#pragma omp critical
			{ pOther->custom_data["head"] += amount_to_transfer; }
		}	
	}

}

void head_migration_direction( Cell* pCell, Phenotype& phenotype, double dt )
{
	phenotype.motility.chemotaxis_direction = parameters.doubles("head_migration_direction"); 
	
	phenotype.motility.migration_speed = parameters.doubles("head_migration_speed"); 
	phenotype.motility.migration_bias = parameters.doubles("head_migration_bias");
	phenotype.motility.persistence_time = parameters.doubles("head_migration_persistence"); 
	
	return chemotaxis_function( pCell,phenotype,dt); 
}

void tail_migration_direction( Cell* pCell, Phenotype& phenotype, double dt )
{
	phenotype.motility.chemotaxis_direction = parameters.doubles("tail_migration_direction"); 
	
	phenotype.motility.migration_speed = parameters.doubles("tail_migration_speed");  
	phenotype.motility.migration_bias = parameters.doubles("tail_migration_bias"); 
	phenotype.motility.persistence_time = parameters.doubles("tail_migration_persistence"); 

	return chemotaxis_function( pCell,phenotype,dt); 
}

void middle_migration_direction( Cell* pCell, Phenotype& phenotype , double dt )
{
	Cell* pUpstream = pCell->state.attached_cells[0]; 

	if( pCell->state.attached_cells[1]->custom_data["head"] > 
		pCell->state.attached_cells[0]->custom_data["head"] )
	{ pUpstream = pCell->state.attached_cells[1]; }
	
	phenotype.motility.migration_speed = parameters.doubles("middle_migration_speed"); 
	phenotype.motility.migration_bias_direction = 
		pUpstream->phenotype.motility.migration_bias_direction;
		
	normalize( &(phenotype.motility.migration_bias_direction) ); 
		
	return; 
}

