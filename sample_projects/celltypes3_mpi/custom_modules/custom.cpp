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

		get_cell_definition("A").functions.update_phenotype = A_phenotype; 
		get_cell_definition("B").functions.update_phenotype = B_phenotype; 
		get_cell_definition("C").functions.update_phenotype = C_phenotype; 
			
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

	std::vector<double> random_position_within(double Xmin,double Xmax,double Ymin,double Ymax,double Zmin,double Zmax,double max_radius)
	{
		std::vector<double> position(3,0.0); 
		double r = max_radius + 1; 
		while( r > max_radius )
		{
			position[0] = Xmin + UniformRandom()* (Xmax - Xmin); 
			position[1] = Ymin + UniformRandom()* (Ymax - Ymin); 
			position[2] = Zmin + UniformRandom()* (Zmax - Zmin); 
			
			r = norm( position ); 
		}
		return position;
	}

	void distribute_type(const TypeCount& type_info, Microenvironment &m, mpi_Environment &world, mpi_Cartesian &cart_topo, double Xmin, double Xmax, double Ymin, double Ymax, double Zmin, double Zmax, double max_radius)
	{
		mpi_CellPositions cp;                
	    mpi_MyCells       mc;                

		if(world.rank == 0)
		{
			std::vector<std::vector<double>> positions_root;
			positions_root.reserve(type_info.count);
			for( int n = 0 ; n < type_info.count ; n++ )
			{
				positions_root.push_back( random_position_within(Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, max_radius) );
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
}

void create_cell_types( mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	create_cell_types_base( &world, &cart_topo );
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
	
	double max_radius = parameters.doubles("max_distance_from_origin");
	if( Xmax > max_radius )
	{ Xmax = max_radius; }
	if( Xmin < -max_radius )
	{ Xmin = -max_radius; }
	
	if( Ymax > max_radius )
	{ Ymax = max_radius; }
	if( Ymin < -max_radius )
	{ Ymin = -max_radius; }

	if( Zmax > max_radius )
	{ Zmax = max_radius; }
	if( Zmin < -max_radius )
	{ Zmin = -max_radius; }
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}

	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 

	Cell* pC;
	
	// place A
	for( int n = 0 ; n < parameters.ints("number_of_A") ; n++ )
	{
		std::vector<double> position = {0,0,0}; 
		
		double r = max_radius + 1; 
		while( r > max_radius )
		{
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			r = norm( position ); 
		}
		
		pC = create_cell( get_cell_definition("A") ); 
		pC->assign_position( position );
		for( int k=0 ; k < pC->phenotype.death.rates.size() ; k++ )
		{ pC->phenotype.death.rates[k] = 0.0; }
	}
	
	// place B
	for( int n = 0 ; n < parameters.ints("number_of_B") ; n++ )
	{
		std::vector<double> position = {0,0,0}; 
		
		double r = max_radius + 1; 
		while( r > max_radius )
		{
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			r = norm( position ); 
		}
		
		pC = create_cell( get_cell_definition("B") ); 
		pC->assign_position( position );
		for( int k=0 ; k < pC->phenotype.death.rates.size() ; k++ )
		{ pC->phenotype.death.rates[k] = 0.0; }
	}

	// place C
	for( int n = 0 ; n < parameters.ints("number_of_C") ; n++ )
	{
		std::vector<double> position = {0,0,0}; 

		double r = max_radius + 1; 
		while( r > max_radius )
		{
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			r = norm( position ); 
		}
		
		pC = create_cell( get_cell_definition("C") ); 
		pC->assign_position( position );
		for( int k=0 ; k < pC->phenotype.death.rates.size() ; k++ )
		{ pC->phenotype.death.rates[k] = 0.0; }
	}

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
	
	double max_radius = parameters.doubles("max_distance_from_origin");
	if( Xmax > max_radius )
	{ Xmax = max_radius; }
	if( Xmin < -max_radius )
	{ Xmin = -max_radius; }
	
	if( Ymax > max_radius )
	{ Ymax = max_radius; }
	if( Ymin < -max_radius )
	{ Ymin = -max_radius; }

	if( Zmax > max_radius )
	{ Zmax = max_radius; }
	if( Zmin < -max_radius )
	{ Zmin = -max_radius; }
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}

	std::vector<TypeCount> types = {
		{"A", parameters.ints("number_of_A")},
		{"B", parameters.ints("number_of_B")},
		{"C", parameters.ints("number_of_C")}
	};

	for( const auto& t : types )
	{
		distribute_type( t, m, world, cart_topo, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, max_radius );
	}

	return; 
}

std::vector<std::string> regular_colors( Cell* pCell )
{
	static int A_type = get_cell_definition( "A" ).type; 
	static int B_type = get_cell_definition( "B" ).type; 
	static int C_type = get_cell_definition( "C" ).type; 
	
	std::vector<std::string> output = {"black" , "black" , "black" , "black"} ;

	if( pCell->type == A_type )
	{
		 output[0] = parameters.strings("A_color");  
		 output[2] = parameters.strings("A_color");  
	}
	
	if( pCell->type == B_type )
	{
		 output[0] = parameters.strings("B_color");  
		 output[2] = parameters.strings("B_color");  
	}
	
	if( pCell->type == C_type )
	{
		 output[0] = parameters.strings("C_color");  
		 output[2] = parameters.strings("C_color");  
	}
	
	if( pCell->phenotype.death.dead == true )
	{
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
		{
			output[2] = "chocolate";
		}
		else
		{
			output[2] = "black"; 
		}
	}

	return output; 
}

std::vector<std::string> pseudo_fluorescence( Cell* pCell )
{
	static int A_type = get_cell_definition( "A" ).type; 
	static int B_type = get_cell_definition( "B" ).type; 
	static int C_type = get_cell_definition( "C" ).type; 
	
	static int nA = microenvironment.find_density_index( "signal A" ); 
	static int nB = microenvironment.find_density_index( "signal B" ); 
	static int nC = microenvironment.find_density_index( "signal C" ); 
	
	static Cell_Definition* pCD_A  = find_cell_definition("A");
	static Cell_Definition* pCD_B  = find_cell_definition("B");
	static Cell_Definition* pCD_C  = find_cell_definition("C");
	
	std::vector<std::string> output = {"black" , "black" , "black" , "black"} ;
	
	char color [32]; 
	double max_fluorescence = 1; 

	double value = 0.0; 

	if( pCell->type == A_type )
	{
		value = pCell->phenotype.secretion.secretion_rates[nA] 
			/ ( 0.001 + pCD_A->phenotype.secretion.secretion_rates[nA] ) ;
			
		value *= (1.0-pCell->phenotype.volume.fluid_fraction) * max_fluorescence;  
		if( pCell->phenotype.death.dead == true )
		{ value = (1.0-pCell->phenotype.volume.fluid_fraction) * max_fluorescence; }
		sprintf( color, "rgba(255,0,255,%f)", value ); 		
	}
	
	if( pCell->type == B_type )
	{
		value = pCell->phenotype.secretion.secretion_rates[nB] 
			/ ( 0.001 + pCD_B->phenotype.secretion.secretion_rates[nB] ); 
		value *= (1.0-pCell->phenotype.volume.fluid_fraction) * max_fluorescence;  
		if( pCell->phenotype.death.dead == true )
		{ value = (1.0-pCell->phenotype.volume.fluid_fraction) * max_fluorescence; }
		sprintf( color, "rgba(0,255,0,%f)", value ); 		
	}
	
	if( pCell->type == C_type )
	{
		value = pCell->phenotype.secretion.secretion_rates[nC] 
			/ ( 0.001 + pCD_C->phenotype.secretion.secretion_rates[nC] ); 
		value *= (1.0-pCell->phenotype.volume.fluid_fraction) * max_fluorescence;  
		if( pCell->phenotype.death.dead == true )
		{ value = (1.0-pCell->phenotype.volume.fluid_fraction) * max_fluorescence; }
		sprintf( color, "rgba(0,255,255,%f)", value ); 		
	}

	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{ sprintf(color,"rgba(0,0,0,%f)" , value ); } 

	output = { color, "none", color , "none" }; 
	return output; 
}

std::vector<std::string> nanohub_fluorescence( Cell* pCell )
{
	static int A_type = get_cell_definition( "A" ).type; 
	static int B_type = get_cell_definition( "B" ).type; 
	static int C_type = get_cell_definition( "C" ).type; 
	
	static int nA = microenvironment.find_density_index( "signal A" ); 
	static int nB = microenvironment.find_density_index( "signal B" ); 
	static int nC = microenvironment.find_density_index( "signal C" ); 
	
	static Cell_Definition* pCD_A  = find_cell_definition("A");
	static Cell_Definition* pCD_B  = find_cell_definition("B");
	static Cell_Definition* pCD_C  = find_cell_definition("C");
	
	std::vector<std::string> output = {"black" , "black" , "black" , "black"} ;
	
	char color [32]; 
	double max_fluorescence = 1; 

	double value = 0.0; 
	
	value = pCell->phenotype.secretion.secretion_rates[nA] 
			/ ( 0.001 + pCD_A->phenotype.secretion.secretion_rates[nA] ) ;
			
	value *= (1.0-pCell->phenotype.volume.fluid_fraction) * max_fluorescence;  
	if( pCell->phenotype.death.dead == true )
	{ value = (1.0-pCell->phenotype.volume.fluid_fraction) * max_fluorescence; }
	sprintf( color, "rgba(255,0,255,%f)", value ); 
	output[0].assign( color ); 	

	value = pCell->phenotype.secretion.secretion_rates[nB] 
			/ ( 0.001 + pCD_B->phenotype.secretion.secretion_rates[nB] ); 
	value *= (1.0-pCell->phenotype.volume.fluid_fraction) * max_fluorescence;  
	if( pCell->phenotype.death.dead == true )
	{ value = (1.0-pCell->phenotype.volume.fluid_fraction) * max_fluorescence; }
	sprintf( color, "rgba(0,255,0,%f)", value ); 		
	output[1].assign( color );

	value = pCell->phenotype.secretion.secretion_rates[nC] 
			/ ( 0.001 + pCD_C->phenotype.secretion.secretion_rates[nC] ); 
	value *= (1.0-pCell->phenotype.volume.fluid_fraction) * max_fluorescence;  
	if( pCell->phenotype.death.dead == true )
	{ value = (1.0-pCell->phenotype.volume.fluid_fraction) * max_fluorescence; }
	sprintf( color, "rgba(0,255,255,%f)", value ); 		
	output[2].assign( color );

	return output; 
}

void SVG_plot_dark( std::string filename , Microenvironment& M, double z_slice , double time, std::vector<std::string> (*cell_coloring_function)(Cell*) )
{
	PhysiCell_SVG_options.SVG_viewbox_min_x = 0;
	PhysiCell_SVG_options.SVG_viewbox_max_x = parameters.ints("plot_x"); 

	PhysiCell_SVG_options.SVG_viewbox_min_y = 0; 
	PhysiCell_SVG_options.SVG_viewbox_max_y = parameters.ints("plot_y"); 

	PhysiCell_SVG_options.x_range = parameters.ints("plot_x") - 20; 
	PhysiCell_SVG_options.y_range = parameters.ints("plot_y") - 20; 
	
	PhysiCell_SVG_options.SVG_plotting_shift_x = -0.5*parameters.ints("plot_x") + 10; 
	PhysiCell_SVG_options.SVG_plotting_shift_y = -0.5*parameters.ints("plot_y") + 10; 

	PhysiCell_SVG_options.background = "black"; 
	PhysiCell_SVG_options.foreground_color = "white"; 
	PhysiCell_SVG_options.marker_border_color = "white"; 
	PhysiCell_SVG_options.plot_nuclei = parameters.bools("plot_nuclei");  
	PhysiCell_SVG_options.plot_nuclei_black = true;  
	
	PhysiCell_SVG_options.fill_color = "rgb(0,0,0)"; 
	PhysiCell_SVG_options.stroke_color = "white"; 
	PhysiCell_SVG_options.frame_thickness = parameters.ints("frame_thickness");  
	
	PhysiCell_SVG_options.agents_as_rectangles = parameters.bools("agents_as_rectangles"); 
	PhysiCell_SVG_options.show_cell_type = parameters.bools("show_cell_type"); 
	PhysiCell_SVG_options.show_scale_bar = parameters.bools("show_scale_bar"); 
	
	PhysiCell_SVG_options.scale_bar_length = parameters.ints("scale_bar_length"); 
	PhysiCell_SVG_options.scale_bar_thickness = parameters.doubles("scale_bar_thickness"); 
	PhysiCell_SVG_options.length_bar = parameters.ints("length_bar"); 
	PhysiCell_SVG_options.show_time = parameters.bools("show_time");
	PhysiCell_SVG_options.show_runtime = parameters.bools("show_runtime"); 
	
	SVG_plot( filename , M, z_slice , time , cell_coloring_function ); 
	return; 
}

// Below are the up_down_signal utilities and phenotype functions (unchanged)

up_down_signal::up_down_signal()
{
    up = 0.0;
	down = 0.0; 
	
	base_parameter = 0.0;
	max_parameter = 1.0; 

    no_promoters = true; 
    no_inhibitors = true; 
	return; 
} 

void up_down_signal::add_effect( double factor, char factor_type )
{
	if( factor_type == 'N' || factor_type == 'n' )
	{ return; }

	if( factor_type == 'P' || factor_type == 'p' )
	{
		up += factor; 
		no_promoters = false; 
		return; 
	}

	if( factor_type == 'I' || factor_type == 'i' )
	{
		down += factor; 
		no_inhibitors = false; 
		return; 
	}

	std::cout << "Warning: unknown factor_type " << factor_type << " in up_down_signal::add_effect()" << std::endl; 
	return; 
} 

void up_down_signal::add_effect( double factor, std::string factor_type )
{
	if( factor_type == "N" || factor_type == "n" || factor_type == "neutral" )
	{ return; }

	if( factor_type == "P" || factor_type == "p" || factor_type == "promoter" )
	{
		up += factor; 
		no_promoters = false; 
		return; 
	}

	if( factor_type == "I" || factor_type == "i" || factor_type == "inhibitor" )
	{
		down += factor; 
		no_inhibitors = false; 
		return; 
	}

	std::cout << "Warning: unknown factor_type " << factor_type << " in up_down_signal::add_effect()" << std::endl; 
	return; 
}

void up_down_signal::display( void )
{
	std::cout << "Up: " << up << " down: " << down << " base: " << base_parameter << " max: " << max_parameter << std::endl; 
	return; 
}

double up_down_signal::compute_effect_hill( void )
{
	double output = base_parameter; 
	double promoter_effect = 0.0; 
	double inhibitor_effect = 0.0; 
	
	if( no_promoters == false )
	{
		promoter_effect = up / ( 1.0 + up ); 
		output += promoter_effect * (max_parameter - base_parameter); 
	}
	
	if( no_inhibitors == false )
	{
		inhibitor_effect = 1.0 / ( 1.0 + down ); 
		output *= inhibitor_effect; 
	}
	
	return output; 
}

double up_down_signal::compute_effect( void )
{
	double output = base_parameter; 
	if( no_promoters == false )
	{
		output += up; 
	}
	if( no_inhibitors == false )
	{
		output -= down; 
	}
	return output; 
}

void up_down_signal::reset( void )
{
	up = 0.0; 
	down = 0.0; 
	no_promoters = true; 
	no_inhibitors = true; 
	return; 
}

void A_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	double dt_phenotype = dt; 
	static int start_phase_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int end_phase_index = start_phase_index; 
	
	double prob = dt_phenotype * parameters.doubles("A_o2_birth_threshold") / ( parameters.doubles("A_o2_birth_threshold") + pCell->nearest_density_vector()[0] ); 
	if( UniformRandom() < prob )
	{
		phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = 0.0; 
	}
	else
	{
		phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = parameters.doubles("A_transition_rate"); 
	}
	
	return; 
}

void B_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	double dt_phenotype = dt; 
	static int start_phase_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int end_phase_index = start_phase_index; 
	
	double prob = dt_phenotype * parameters.doubles("B_o2_birth_threshold") / ( parameters.doubles("B_o2_birth_threshold") + pCell->nearest_density_vector()[0] ); 
	if( UniformRandom() < prob )
	{
		phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = 0.0; 
	}
	else
	{
		phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = parameters.doubles("B_transition_rate"); 
	}
	
	return; 
}

void C_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	double dt_phenotype = dt; 
	static int start_phase_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int end_phase_index = start_phase_index; 
	
	double prob = dt_phenotype * parameters.doubles("C_o2_birth_threshold") / ( parameters.doubles("C_o2_birth_threshold") + pCell->nearest_density_vector()[0] ); 
	if( UniformRandom() < prob )
	{
		phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = 0.0; 
	}
	else
	{
		phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = parameters.doubles("C_transition_rate"); 
	}
	
	return; 
}
