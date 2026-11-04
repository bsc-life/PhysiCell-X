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
#define _USE_MATH_DEFINES
#include <cmath>

#include "./custom.h"
#include "../DistPhy/DistPhy_Utils.h"
#include "../DistPhy/DistPhy_Collective.h"
#include <algorithm>  // ✅ Pour std::find
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace DistPhy::mpi;

bool injected = false;

double dynamic_cart_x_left = -450.0;
double dynamic_cart_x_right = 450.0;
double dynamic_cart_y_min = -1020.0;
double dynamic_cart_y_max = 1020.0;

namespace
{
double safe_radial_distance(double x, double y)
{
    // Several boundary rules divide by r; at the origin the radial force is zero.
    return std::max(std::hypot(x, y), 1e-12);
}
}


void create_cell_types( mpi_Environment& world, mpi_Cartesian& cart_topo )
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
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	cell_defaults.functions.update_velocity_parallel = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL;  
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	
	
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(world, cart_topo);

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(world, cart_topo);

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(world, cart_topo);

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	
	cell_defaults.functions.contact_function = contact_function;

	if (parameters.bools("square_real_dim_bc") == true)
	{
		find_cell_definition("FLS")->functions.custom_cell_rule = realdim_function_square;
		find_cell_definition("Osteoclast")->functions.custom_cell_rule = realdim_function_square;
		find_cell_definition("Chondrocytes")->functions.custom_cell_rule = realdim_function_square;
		find_cell_definition("M1")->functions.custom_cell_rule = realdim_function_square;
		find_cell_definition("M2")->functions.custom_cell_rule = realdim_function_square;
		find_cell_definition("LTCD4+")->functions.custom_cell_rule = realdim_function_square;
		find_cell_definition("LB")->functions.custom_cell_rule = realdim_function_square;
		find_cell_definition("Plasmocytes")->functions.custom_cell_rule = realdim_function_square;
		find_cell_definition("LTCD8+")->functions.custom_cell_rule = realdim_function_square;
		find_cell_definition("LTreg")->functions.custom_cell_rule = realdim_function_square;
		find_cell_definition("LTh17")->functions.custom_cell_rule = realdim_function_square;
		find_cell_definition("LTh1")->functions.custom_cell_rule = realdim_function_square;
		find_cell_definition("LTh2")->functions.custom_cell_rule = realdim_function_square;
	}
	
	else if (parameters.bools("circle_real_dim_bc") == true)
	{
		find_cell_definition("FLS")->functions.custom_cell_rule = realdim_function_circle;
		find_cell_definition("Osteoclast")->functions.custom_cell_rule = realdim_function_circle;
		find_cell_definition("Chondrocytes")->functions.custom_cell_rule = realdim_function_circle;
		find_cell_definition("M1")->functions.custom_cell_rule = realdim_function_circle;
		find_cell_definition("M2")->functions.custom_cell_rule = realdim_function_circle;
		find_cell_definition("LTCD4+")->functions.custom_cell_rule = realdim_function_circle;
		find_cell_definition("LB")->functions.custom_cell_rule = realdim_function_circle;
		find_cell_definition("Plasmocytes")->functions.custom_cell_rule = realdim_function_circle;
		find_cell_definition("LTCD8+")->functions.custom_cell_rule = realdim_function_circle;
		find_cell_definition("LTreg")->functions.custom_cell_rule = realdim_function_circle;
		find_cell_definition("LTh17")->functions.custom_cell_rule = realdim_function_circle;
		find_cell_definition("LTh1")->functions.custom_cell_rule = realdim_function_circle;
		find_cell_definition("LTh2")->functions.custom_cell_rule = realdim_function_circle;

	}
	
	else if (parameters.bools("init_circular_bc") == true)
	{
		find_cell_definition("FLS")->functions.custom_cell_rule = init_circular_function;
		find_cell_definition("Osteoclast")->functions.custom_cell_rule = init_circular_function;
		find_cell_definition("Chondrocytes")->functions.custom_cell_rule = init_circular_function;
		find_cell_definition("M1")->functions.custom_cell_rule = init_circular_function;
		find_cell_definition("M2")->functions.custom_cell_rule = init_circular_function;
		find_cell_definition("LTCD4+")->functions.custom_cell_rule = init_circular_function;
		find_cell_definition("LB")->functions.custom_cell_rule = init_circular_function;
		find_cell_definition("Plasmocytes")->functions.custom_cell_rule = init_circular_function;
		find_cell_definition("LTCD8+")->functions.custom_cell_rule = init_circular_function;
		find_cell_definition("LTreg")->functions.custom_cell_rule = init_circular_function;
		find_cell_definition("LTh17")->functions.custom_cell_rule = init_circular_function;
		find_cell_definition("LTh1")->functions.custom_cell_rule = init_circular_function;
		find_cell_definition("LTh2")->functions.custom_cell_rule = init_circular_function;

	}
	
	else if (parameters.bools("init_square_bc") == true)
	{
		find_cell_definition("FLS")->functions.custom_cell_rule = init_square_function;
		find_cell_definition("Osteoclast")->functions.custom_cell_rule = init_square_function;
		find_cell_definition("Chondrocytes")->functions.custom_cell_rule = init_square_function;
		find_cell_definition("M1")->functions.custom_cell_rule = init_square_function;
		find_cell_definition("M2")->functions.custom_cell_rule = init_square_function;
		find_cell_definition("LTCD4+")->functions.custom_cell_rule = init_square_function;
		find_cell_definition("LB")->functions.custom_cell_rule = init_square_function;
		find_cell_definition("Plasmocytes")->functions.custom_cell_rule = init_square_function;
		find_cell_definition("LTCD8+")->functions.custom_cell_rule = init_square_function;
		find_cell_definition("LTreg")->functions.custom_cell_rule = init_square_function;
		find_cell_definition("LTh17")->functions.custom_cell_rule = init_square_function;
		find_cell_definition("LTh1")->functions.custom_cell_rule = init_square_function;
		find_cell_definition("LTh2")->functions.custom_cell_rule = init_square_function;

	}
	
	else if (parameters.bools.find_index("final_square_bc") != -1 && parameters.bools("final_square_bc") == true)
{
	find_cell_definition("FLS")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("Osteoclast")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("Chondrocytes")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("M1")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("M2")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("LTCD4+")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("LB")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("Plasmocytes")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("LTCD8+")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("LTreg")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("LTh17")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("LTh1")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("LTh2")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("M0")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("Collagen")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("Joint")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("MoDC")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("Monocytes")->functions.custom_cell_rule = final_square_function;
	find_cell_definition("Neutrophils")->functions.custom_cell_rule = final_square_function;
}
	else if (parameters.bools.find_index("BCandinj") != -1 && parameters.bools("BCandinj") == true)
{
	find_cell_definition("FLS")->functions.custom_cell_rule = combined_function;
	find_cell_definition("Osteoclast")->functions.custom_cell_rule = combined_function;
	find_cell_definition("Chondrocytes")->functions.custom_cell_rule = combined_function;
	find_cell_definition("M1")->functions.custom_cell_rule = combined_function;
	find_cell_definition("M2")->functions.custom_cell_rule = combined_function;
	find_cell_definition("LTCD4+")->functions.custom_cell_rule = combined_function;
	find_cell_definition("LB")->functions.custom_cell_rule = combined_function;
	find_cell_definition("Plasmocytes")->functions.custom_cell_rule = combined_function;
	find_cell_definition("LTCD8+")->functions.custom_cell_rule = combined_function;
	find_cell_definition("LTreg")->functions.custom_cell_rule = combined_function;
	find_cell_definition("LTh17")->functions.custom_cell_rule = combined_function;
	find_cell_definition("LTh1")->functions.custom_cell_rule = combined_function;
	find_cell_definition("LTh2")->functions.custom_cell_rule = combined_function;
	find_cell_definition("M0")->functions.custom_cell_rule = combined_function;
	find_cell_definition("Collagen")->functions.custom_cell_rule = combined_function;
	find_cell_definition("Joint")->functions.custom_cell_rule = combined_function;
	find_cell_definition("MoDC")->functions.custom_cell_rule = combined_function;
	find_cell_definition("Monocytes")->functions.custom_cell_rule = combined_function;
	find_cell_definition("Neutrophils")->functions.custom_cell_rule = combined_function;
	find_cell_definition("bone")->functions.custom_cell_rule = combined_function;
	

}
	else
	{
		
		display_cell_definitions(std::cout, world, cart_topo);
		return;
	}

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions(std::cout, world, cart_topo);
	
	return; 
}




void setup_microenvironment( mpi_Environment& world, mpi_Cartesian& cart_topo )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(world, cart_topo);
	
	return; 
}

void setup_tissue( Microenvironment& m, mpi_Environment& world, mpi_Cartesian& cart_topo )
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
	
	for( int k=0; k < cell_definitions_by_index.size() && parameters.ints("number_of_cells") > 0 ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		mpi_CellPositions cp;
		mpi_MyCells mc;
		if (IOProcessor(world))
		{
			std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
			std::vector<std::vector<double>> positions;
			for (int n = 0; n < parameters.ints("number_of_cells"); ++n)
			{
				positions.push_back({Xmin + UniformRandom()*Xrange,
				                     Ymin + UniformRandom()*Yrange,
				                     Zmin + UniformRandom()*Zrange});
			}
			int start_id = Basic_Agent::get_max_ID_in_parallel();
			cp.positions_to_rank_list(positions,
				m.mesh.bounding_box[0], m.mesh.bounding_box[3],
				m.mesh.bounding_box[1], m.mesh.bounding_box[4],
				m.mesh.bounding_box[2], m.mesh.bounding_box[5],
				m.mesh.dx, m.mesh.dy, m.mesh.dz, world, cart_topo, start_id);
			Basic_Agent::set_max_ID_in_parallel(start_id + positions.size());
		}
		distribute_cell_positions(cp, mc, world, cart_topo);
		for (int i = 0; i < mc.my_no_of_cell_IDs; ++i)
		{
			Cell* cell = create_cell(*pCD, mc.my_cell_IDs[i]);
			cell->assign_position(mc.my_cell_coords[3*i], mc.my_cell_coords[3*i+1],
			                      mc.my_cell_coords[3*i+2], world, cart_topo);
		}
	}
	if (IOProcessor(world)) std::cout << std::endl;
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(world, cart_topo);

	// The MPI CSV loader constructs cells with rank-local IDs. Reassign all
	// initial cells, including generated ones, before halo exchange uses IDs.
	const int local_count = static_cast<int>(all_cells->size());
	int first_id = 0;
	int global_count = 0;
	MPI_Exscan(&local_count, &first_id, 1, MPI_INT, MPI_SUM, cart_topo.mpi_cart_comm);
	if (world.rank == 0) first_id = 0;
	MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, cart_topo.mpi_cart_comm);
	int invalid_local = 0;
	for (int i = 0; i < local_count; ++i)
	{
		Cell* cell = (*all_cells)[i];
		cell->ID = first_id + i;
		if (cell->position.size() != 3 ||
		    !std::isfinite(cell->position[0]) ||
		    !std::isfinite(cell->position[1]) ||
		    !std::isfinite(cell->position[2]))
		{
			if (invalid_local++ < 5)
				std::cerr << "Rank " << world.rank << ": invalid initial position for cell " << cell->ID << std::endl;
		}
	}
	Basic_Agent::set_max_ID_in_parallel(global_count);
	int invalid_global = 0;
	MPI_Allreduce(&invalid_local, &invalid_global, 1, MPI_INT, MPI_SUM, cart_topo.mpi_cart_comm);
	if (invalid_global > 0)
		throw std::runtime_error("nonfinite initial cell positions");
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }




void update_dynamic_boundaries(mpi_Cartesian& cart_topo)
{
    dynamic_cart_x_right = 450;
    dynamic_cart_x_left = -450;
    dynamic_cart_y_max = 1020;
    dynamic_cart_y_min = -1020;

    for (auto& cell : *all_cells)
    {
        if (
			(cell->type == get_cell_definition("Chondrocytes").type ||
			 cell->type == get_cell_definition("Collagen").type) &&
			 cell->phenotype.death.dead == false
		)

        {
            double x = cell->position[0];
            double y = cell->position[1];
			
			if (y > 0)
				{
					if (y > dynamic_cart_y_max) dynamic_cart_y_max = y;
				}
            
			if (y < 0)
				{
					if (y < dynamic_cart_y_min) dynamic_cart_y_min = y;
				}
            
			
			// side
			//left
			if (y >= -1015.0 && y <= 1015.0)
				{
					if (x < 0)
					{
						if (x > dynamic_cart_x_left) dynamic_cart_x_left = x;
					}
				}
			
			//right
			if (y >= -1015.0 && y <= 1015.0)
				{
					if (x > 0)
					{
						if (x < dynamic_cart_x_right) dynamic_cart_x_right = x;
					}
				}
			
            
        }
    }
	// The boundary is a property of the whole domain, including cells on other ranks.
	double local_max[2] = {dynamic_cart_x_left, dynamic_cart_y_max};
	double local_min[2] = {dynamic_cart_x_right, dynamic_cart_y_min};
	double global_max[2], global_min[2];
	MPI_Allreduce(local_max, global_max, 2, MPI_DOUBLE, MPI_MAX, cart_topo.mpi_cart_comm);
	MPI_Allreduce(local_min, global_min, 2, MPI_DOUBLE, MPI_MIN, cart_topo.mpi_cart_comm);
	dynamic_cart_x_left = global_max[0];
	dynamic_cart_y_max = global_max[1];
	dynamic_cart_x_right = global_min[0];
	dynamic_cart_y_min = global_min[1];
}

void update_model_environment(mpi_Environment& world, mpi_Cartesian& cart_topo)
{
	if (parameters.bools.find_index("BCandinj") != -1 && parameters.bools("BCandinj"))
	{
		if (!injected && PhysiCell_globals.current_time >= parameters.doubles("tofa_injection_time"))
		{
			const int tofa_index = microenvironment.find_density_index("tofacitinib");
			const double concentration = parameters.doubles("tofa_concentration");
			for (int n = 0; n < microenvironment.number_of_voxels(); ++n)
				microenvironment(n)[tofa_index] = concentration;
			injected = true;
			if (IOProcessor(world))
				std::cout << "Tofacitinib injected at t = " << PhysiCell_globals.current_time / 60.0 << " hours" << std::endl;
		}
		update_dynamic_boundaries(cart_topo);
	}
	treatment_function(world);
}
void trigger_apoptosis(Cell* cell)
{
	if (!cell->phenotype.death.dead)
	{
		int model = cell->phenotype.death.find_death_model_index("apoptosis");

		// Ne déclenche que si le taux est > 0, donc respectera la règle CSV dynamique
		if (model >= 0)
		{
			cell->start_death(model);
		}
	}
}


void combined_function(Cell* pCell, Phenotype& phenotype, double dt)
{
	if (parameters.bools("BCandinj") == false)
        return;
	
	double x = pCell->position[0];
    double y = pCell->position[1];
    double r = safe_radial_distance(x, y);
	
	if (pCell->type == get_cell_definition("Chondrocytes").type ||
        pCell->type == get_cell_definition("Collagen").type)
	
	{
		// side
		if (y >= dynamic_cart_y_min && y <= dynamic_cart_y_max && x <= dynamic_cart_x_right+10 && x >= dynamic_cart_x_right)
		{	
			trigger_apoptosis(pCell);
		}
					
		if (y >= dynamic_cart_y_min && y <= dynamic_cart_y_max && x >= dynamic_cart_x_left-10 && x <= dynamic_cart_x_left )
		{	
			trigger_apoptosis(pCell);
		}
        	
		
		//top
		if (x <= dynamic_cart_x_left)
		{
			if (y <= dynamic_cart_y_max && y >= dynamic_cart_y_max-20)
			{	
			trigger_apoptosis(pCell);
			}
		}
		
		if (x >= dynamic_cart_x_right)
		{
			if (y <= dynamic_cart_y_max && y >= dynamic_cart_y_max-20)
			{	
			trigger_apoptosis(pCell);
			}
		}
		
		//bottom
		if (x <= dynamic_cart_x_left)
		{
			if (y >= dynamic_cart_y_min && y <= dynamic_cart_y_min+20)
			{	
			trigger_apoptosis(pCell);
			}
		}
		
		if (x >= dynamic_cart_x_right)
		{
			if (y >= dynamic_cart_y_min && y <= dynamic_cart_y_min+20)
			{	
			trigger_apoptosis(pCell);
			}
		}
		
	}
	
		
    double push_strength = 1.0; // tu peux ajuster ce paramètre

    if (pCell->type == get_cell_definition("Osteoclast").type)
    {
		if (y > 1020.0)
		{
            pCell->velocity[1] -= push_strength;
		}
		if (y < -1020.0)
		{
            pCell->velocity[1] += push_strength;
		}
		
        if (x > 0.0)
        { 
			if (x < 450.0)
			{
            // Applique une force vers l'exterieur
            pCell->velocity[0] += push_strength;
			}
			
        }
		if (x < 0.0)
        { 
			if (x > -450.0)
			{
            pCell->velocity[0] -= push_strength;
			}
        }
    }
	
	if (pCell->type == get_cell_definition("Chondrocytes").type ||
        pCell->type == get_cell_definition("Collagen").type)
	
	{
		// global 
        if (r > 2020.0)
        {
            
            pCell->velocity[0] -= push_strength * (x / r);
            pCell->velocity[1] -= push_strength * (y / r);
        }
		
		// side
		
        if (y >= -1015.0 && y <= 1015.0)
			{
				if (x >= 450.0)
				{
				// Applique une force vers l'exterieur
				pCell->velocity[0] -= push_strength;
				}
			}
        
		
		if (y >= -1015.0 && y <= 1015.0)
			{
				if (x <= -450.0)
				{
				pCell->velocity[0] += push_strength;
				}
			}
        
		
		
		//top
		if (x < -468.0)
		{
			if (y > 1015.0 && y < 1020)
			{
				pCell->velocity[1] += push_strength;
			}
		}
		
		if (x > 468.0)
		{
			if (y > 1015.0 && y < 1020)
			{
				pCell->velocity[1] += push_strength;
			}
		}
		
		//bottom
		if (x < -468.0)
		{
			if (y > -1020 && y < -1015)
			{
				pCell->velocity[1] -= push_strength;
			}
		}
		
		if (x > 468.0)
		{
			if (y > -1020 && y < -1015)
			{
				pCell->velocity[1] -= push_strength;
			}
		}
		
		
	}
	
	if (pCell->type == get_cell_definition("FLS").type||
        pCell->type == get_cell_definition("M0").type)
    {
        
		if (r < 2025.0)
        {
            // Applique une force vers l'exterieur
            pCell->velocity[0] += push_strength * (x / r);
            pCell->velocity[1] += push_strength * (y / r);
        }
		
    }
	
	if (pCell->type == get_cell_definition("Joint").type)
    {
        if (r > 2058.0)
		{
			if (r < 2060.0)
			{
				// Applique une force vers l'exterieur
				pCell->velocity[0] += push_strength * (x / r);
				pCell->velocity[1] += push_strength * (y / r);
			}
		}
		
		// gauche haut
		
		if (x > 0 && y>0)
		{
			if (x < 1468 && y <=1420 ) {
				
			{
			pCell->velocity[0] += push_strength;
			}
		}
		}
		
		if (x > 0 && y>0)
		{
			if (x >= 1468 && y <1020 ) {
				
			{
			pCell->velocity[1] += push_strength;
			}
		}
		}
		
		// gauche bas
		if (x > 0 && y<0)
		{
			if (x < 1468 && y >=-1420 ) {
				
			{
			pCell->velocity[0] += push_strength;
			}
		}
		}
		
		if (x > 0 && y<0)
		{
			if (x >= 1468 && y >-1020 ) {
				
			{
			pCell->velocity[1] -= push_strength;
			}
		}
		}
		
		// droite haut
		
		if (x < 0 && y>0)
		{
			if (x > -1468 && y <=1420 ) {
				
			{
			pCell->velocity[0] -= push_strength;
			}
		}
		}
		
		if (x < 0 && y>0)
		{
			if (x <= -1468 && y <1020 ) {
				
			{
			pCell->velocity[1] += push_strength;
			}
		}
		}
		
		// droite bas
		if (x < 0 && y<0)
		{
			if (x > -1468 && y >=-1420 ) {
				
			{
			pCell->velocity[0] -= push_strength;
			}
		}
		}
		
		if (x < 0 && y<0)
		{
			if (x <= -1468 && y >-1020 ) {
				
			{
			pCell->velocity[1] -= push_strength;
			}
		}
		}
    }
	
    if (
        pCell->type == get_cell_definition("M1").type ||
        pCell->type == get_cell_definition("M2").type ||
        pCell->type == get_cell_definition("LTCD4+").type ||
        pCell->type == get_cell_definition("LB").type ||
        pCell->type == get_cell_definition("Plasmocytes").type ||
        pCell->type == get_cell_definition("LTCD8+").type ||
        pCell->type == get_cell_definition("LTreg").type ||
        pCell->type == get_cell_definition("LTh17").type ||
        pCell->type == get_cell_definition("LTh1").type ||
        pCell->type == get_cell_definition("LTh2").type ||
        pCell->type == get_cell_definition("MoDC").type ||
        pCell->type == get_cell_definition("Monocytes").type ||
        pCell->type == get_cell_definition("Neutrophils").type)
    {
        
		if (r > 2020.0)
        {
            // Applique une force vers le centre du grand cercle
            pCell->velocity[0] -= push_strength * (x / r);
            pCell->velocity[1] -= push_strength * (y / r);
        }
		
		//centre
		if (y >= dynamic_cart_y_min && y <= dynamic_cart_y_max && x >= dynamic_cart_x_right )
		{
			pCell->velocity[0] -= push_strength;
		}
		if (y >= dynamic_cart_y_min && y <= dynamic_cart_y_max && x <= dynamic_cart_x_left )
		{
			pCell->velocity[0] += push_strength;
		}
		
		//haut à gauche
		if (x >= -1468.0 && x <= dynamic_cart_x_left && y >0 )
		{ 
			if (y <= dynamic_cart_y_max)
			{
			pCell->velocity[1] += push_strength;
			}
		}
		
		// haut à droite
		if (x <= 1468.0 && x >= dynamic_cart_x_right && y >0 )
		{ 
			if (y <= dynamic_cart_y_max)
			{
			pCell->velocity[1] += push_strength;
			}
		}
	
		// bas à droite
		if (x <= 1468.0 && x >= dynamic_cart_x_right && y <0 )
		{ 
			if (y >= dynamic_cart_y_min)
			{
			pCell->velocity[1] -= push_strength;
			}
		}
		
		// bas à gauche
		if (x <= -1468.0 && x >= dynamic_cart_x_left && y <0 )
		{ 
			if (y >= dynamic_cart_y_min)
			{
			pCell->velocity[1] -= push_strength;
			}
		}
		
    }


	
}


void realdim_function_square( Cell* pCell, Phenotype& phenotype , double dt )
{
	if (parameters.bools("square_real_dim_bc") == false)
        return;

    double x = pCell->position[0];
    double y = pCell->position[1];
    double r = safe_radial_distance(x, y);

    double push_strength = 1.0; // tu peux ajuster ce paramètre

    if (pCell->type == get_cell_definition("Osteoclast").type)
    {
		if (y > 1020.0)
		{
            pCell->velocity[1] -= push_strength;
		}
		if (y < -1020.0)
		{
            pCell->velocity[1] += push_strength;
		}
		
        if (x > 0.0)
        { 
			if (x < 450.0)
			{
            // Applique une force vers l'exterieur
            pCell->velocity[0] += push_strength;
			}
			
        }
		if (x < 0.0)
        { 
			if (x > -450.0)
			{
            pCell->velocity[0] -= push_strength;
			}
        }
    }
	
	if (pCell->type == get_cell_definition("Chondrocytes").type)
	{
		// side
		
        if (x > 0.0)
        { 
			if (y >= -1020.0 && y <= 1020.0)
			{
				if (x > 450.0)
				{
				// Applique une force vers l'exterieur
				pCell->velocity[0] -= push_strength;
				}
			}
        }
		
		if (x < 0.0)
        { 
			if (y >= -1020.0 && y <= 1020.0)
			{
				if (x < -450.0)
				{
				pCell->velocity[0] += push_strength;
				}
			}
        }
		
		
		//top
		if (x < -468.0)
		{
			if (y > 0)
			{
				if (y < 1020)
			    {
				pCell->velocity[1] += push_strength;
			    }
			}
		}
		
		if (x > 468.0)
		{
			if (y > 0)
			{
				if (y < 1020)
				{
					pCell->velocity[1] += push_strength;
				}
			}
		}
		
		//bottom
		if (x < -468.0)
		{
			if (y < 0)
			{
				if (y > -1020)
				{
					pCell->velocity[1] -= push_strength;
				}
			}
		}
		
		if (x > 468.0)
		{
			if (y < 0)
			{
				if (y > -1020)
				{
					pCell->velocity[1] -= push_strength;
				}
			}
		}
		
		
	}
	
	
	if (pCell->type == get_cell_definition("FLS").type)
    {
        /* for v1
		if (r < 400.0)
        {
            // Applique une force vers l'exterieur
            pCell->velocity[0] += push_strength * (x / r);
            pCell->velocity[1] += push_strength * (y / r);
        }
		*/
		
		if (y > 0.0)
        { 
			if (y < 1460.0)
			{
            // Applique une force vers l'exterieur
            pCell->velocity[1] += push_strength;
			}
			
        }
		if (y < 0.0)
        { 
			if (y > -1460.0)
			{
            // Applique une force vers l'exterieur
            pCell->velocity[1] -= push_strength;
			}
			
        }
    }
	
    if (
        pCell->type == get_cell_definition("M1").type ||
        pCell->type == get_cell_definition("M2").type ||
        pCell->type == get_cell_definition("LTCD4+").type ||
        pCell->type == get_cell_definition("LB").type ||
        pCell->type == get_cell_definition("Plasmocytes").type ||
        pCell->type == get_cell_definition("LTCD8+").type ||
        pCell->type == get_cell_definition("LTreg").type ||
        pCell->type == get_cell_definition("LTh17").type ||
        pCell->type == get_cell_definition("LTh1").type ||
        pCell->type == get_cell_definition("LTh2").type )
    {
        /*
		if (r > 395.0)
        {
            // Applique une force vers le centre du grand cercle
            pCell->velocity[0] -= push_strength * (x / r);
            pCell->velocity[1] -= push_strength * (y / r);
        }
		*/
		//centre
		if (y >= -1420.0 && y <= 1420.0 && x > 50.0 )
		{
			pCell->velocity[0] -= push_strength;
		}
		if (y >= -1420.0 && y <= 1420.0 && x < -50.0 )
		{
			pCell->velocity[0] += push_strength;
		}
		
		//haut à gauche
		if (x >= -1468.0 && x <= -50.0 && y >0 )
		{ 
			if (y < 1420)
			{
			pCell->velocity[1] += push_strength;
			}
		}
		
		// haut à droite
		if (x <= 1468.0 && x >= 50.0 && y >0 )
		{ 
			if (y < 1420)
			{
			pCell->velocity[1] += push_strength;
			}
		}
	
		// bas à droite
		if (x <= 1468.0 && x >= 50.0 && y <0 )
		{ 
			if (y > -1420)
			{
			pCell->velocity[1] -= push_strength;
			}
		}
		
		// bas à gauche
		if (x <= -1468.0 && x >= -50.0 && y <0 )
		{ 
			if (y > -1420)
			{
			pCell->velocity[1] -= push_strength;
			}
		}
		
		//global
		if (y < -1460)
		{ 
			pCell->velocity[1] += push_strength;
		}
		if (y > 1460)
		{ 
			pCell->velocity[1] -= push_strength;
		}
		
    }
} 

void init_square_function( Cell* pCell, Phenotype& phenotype , double dt )
{
	if (parameters.bools("init_square_bc") == false)
        return;

    double x = pCell->position[0];
    double y = pCell->position[1];
    double r = safe_radial_distance(x, y);

    double push_strength = 1.0; // tu peux ajuster ce paramètre

    if (pCell->type == get_cell_definition("Osteoclast").type)
    {
		if (y > 200.0)
		{
            pCell->velocity[1] -= push_strength;
		}
		if (y < -200.0)
		{
            pCell->velocity[1] += push_strength;
		}
		
        if (x > 0.0)
        { 
			if (x < 190.0)
			{
            // Applique une force vers l'exterieur
            pCell->velocity[0] += push_strength;
			}
			
        }
		if (x < 0.0)
        { 
			if (x > -190.0)
			{
            pCell->velocity[0] -= push_strength;
			}
        }
    }
	
	if (pCell->type == get_cell_definition("Chondrocytes").type)
	{
		// global 
        if (r > 400.0)
        {
            
            pCell->velocity[0] -= push_strength * (x / r);
            pCell->velocity[1] -= push_strength * (y / r);
        }
		
		// side
		
        if (x > 0.0)
        { 
			if (y >= -220.0 && y <= 220.0)
			{
				if (x > 180.0)
				{
				// Applique une force vers l'exterieur
				pCell->velocity[0] -= push_strength;
				}
			}
        }
		
		if (x < 0.0)
        { 
			if (y >= -220.0 && y <= 220.0)
			{
				if (x < -180.0)
				{
				pCell->velocity[0] += push_strength;
				}
			}
        }
		
		
		//top
		if (x < -180.0)
		{
			if (y < 210)
			{
				pCell->velocity[1] += push_strength;
			}
		}
		
		if (x > 180.0)
		{
			if (y < 210)
			{
				pCell->velocity[1] += push_strength;
			}
		}
		
		//bottom
		if (x < -180.0)
		{
			if (y > -210)
			{
				pCell->velocity[1] -= push_strength;
			}
		}
		
		if (x > 180.0)
		{
			if (y > -210)
			{
				pCell->velocity[1] -= push_strength;
			}
		}
		
		
	}
	
	
	if (pCell->type == get_cell_definition("FLS").type)
    {
        
		if (r < 400.0)
        {
            // Applique une force vers l'exterieur
            pCell->velocity[0] += push_strength * (x / r);
            pCell->velocity[1] += push_strength * (y / r);
        }
		
    }
	
    if (
        pCell->type == get_cell_definition("M1").type ||
        pCell->type == get_cell_definition("M2").type ||
        pCell->type == get_cell_definition("LTCD4+").type ||
        pCell->type == get_cell_definition("LB").type ||
        pCell->type == get_cell_definition("Plasmocytes").type ||
        pCell->type == get_cell_definition("LTCD8+").type ||
        pCell->type == get_cell_definition("LTreg").type ||
        pCell->type == get_cell_definition("LTh17").type ||
        pCell->type == get_cell_definition("LTh1").type ||
        pCell->type == get_cell_definition("LTh2").type )
    {
        
		if (r > 395.0)
        {
            // Applique une force vers le centre du grand cercle
            pCell->velocity[0] -= push_strength * (x / r);
            pCell->velocity[1] -= push_strength * (y / r);
        }
		
		//centre
		if (y >= -220.0 && y <= 220.0 && x > 170.0 )
		{
			pCell->velocity[0] -= push_strength;
		}
		if (y >= -220.0 && y <= 220.0 && x < -170.0 )
		{
			pCell->velocity[0] += push_strength;
		}
		
		//haut à gauche
		if (x >= -328.0 && x <= -170.0 && y >0 )
		{ 
			if (y < 220)
			{
			pCell->velocity[1] += push_strength;
			}
		}
		
		// haut à droite
		if (x <= 328.0 && x >= 170.0 && y >0 )
		{ 
			if (y < 220)
			{
			pCell->velocity[1] += push_strength;
			}
		}
	
		// bas à droite
		if (x <= 328.0 && x >= 170.0 && y <0 )
		{ 
			if (y > -220)
			{
			pCell->velocity[1] -= push_strength;
			}
		}
		
		// bas à gauche
		if (x <= -328.0 && x >= -170.0 && y <0 )
		{ 
			if (y > -220)
			{
			pCell->velocity[1] -= push_strength;
			}
		}
		
    }
} 


void init_circular_function(Cell* pCell, Phenotype& phenotype, double dt)
{
    // Vérifie si on doit activer le placement circulaire
    if (parameters.bools("init_circular_bc") == false)
        return;

    double x = pCell->position[0];
    double y = pCell->position[1];
    double r = safe_radial_distance(x, y);

    double push_strength = 1; // tu peux ajuster ce paramètre

    // Chondrocytes r < 170
    if (pCell->type == get_cell_definition("Chondrocytes").type)
    {
        if (r < 110.0)
        {
            // Applique une force vers l'exterieur
            pCell->velocity[0] += push_strength * (x / r);
            pCell->velocity[1] += push_strength * (y / r);
        }
    }
	
	if (
		pCell->type == get_cell_definition("Osteoclast").type)

	{
        if (r > 100.0)
        {
            
            pCell->velocity[0] -= push_strength * (x / r);
            pCell->velocity[1] -= push_strength * (y / r);
        }	
		
	}
    
	if (pCell->type == get_cell_definition("FLS").type)
    {
        if (r < 450.0)
        {
            // Applique une force vers l'exterieur
            pCell->velocity[0] += push_strength * (x / r);
            pCell->velocity[1] += push_strength * (y / r);
        }
    }
	
    if (
        pCell->type == get_cell_definition("M1").type ||
        pCell->type == get_cell_definition("M2").type ||
        pCell->type == get_cell_definition("LTCD4+").type ||
        pCell->type == get_cell_definition("LB").type ||
        pCell->type == get_cell_definition("Plasmocytes").type ||
        pCell->type == get_cell_definition("LTCD8+").type ||
        pCell->type == get_cell_definition("LTreg").type ||
        pCell->type == get_cell_definition("LTh17").type ||
        pCell->type == get_cell_definition("LTh1").type ||
        pCell->type == get_cell_definition("LTh2").type )
    {
        if (r < 155.0)
        {
            // Applique une force vers l’extérieur du petit cercle
            pCell->velocity[0] += push_strength * (x / r);
            pCell->velocity[1] += push_strength * (y / r);
        }
        if (r > 445.0)
        {
            // Applique une force vers le centre du grand cercle
            pCell->velocity[0] -= push_strength * (x / r);
            pCell->velocity[1] -= push_strength * (y / r);
        }
    }
}


void realdim_function_circle(Cell* pCell, Phenotype& phenotype, double dt)
{
	
    if (parameters.bools("circle_real_dim_bc") == false)
        return;

    double x = pCell->position[0];
    double y = pCell->position[1];
    double r = safe_radial_distance(x, y);
    double push_strength = 1.0;
	
	// Variables globales pour l'évolution dynamique des rayons FLS
	double fls_inner_radius = 1504.0;
	double fls_outer_radius = 1544.0;
	bool fls_trigger_expand = false;


    // Chondrocytes confinement
    if (pCell->type == get_cell_definition("Chondrocytes").type)
    {
        if (r < 1062.0)
        {
            pCell->velocity[0] += push_strength * (x / r);
            pCell->velocity[1] += push_strength * (y / r);
        }
    }

    // Osteoclast confinement
    if (pCell->type == get_cell_definition("Osteoclast").type)
    {
        if (r > 1061.0)
        {
            pCell->velocity[0] -= push_strength * (x / r);
            pCell->velocity[1] -= push_strength * (y / r);
        }
    }

    // FLS confinement
    if (pCell->type == get_cell_definition("FLS").type)
    {
        if (r < fls_inner_radius)
        {
            pCell->velocity[0] += push_strength * (x / r);
            pCell->velocity[1] += push_strength * (y / r);
        }
    }

    // Cellules du liquide synovial
    if (
        pCell->type == get_cell_definition("M1").type ||
        pCell->type == get_cell_definition("M2").type ||
        pCell->type == get_cell_definition("LTCD4+").type ||
        pCell->type == get_cell_definition("LB").type ||
        pCell->type == get_cell_definition("Plasmocytes").type ||
        pCell->type == get_cell_definition("LTCD8+").type ||
        pCell->type == get_cell_definition("LTreg").type ||
        pCell->type == get_cell_definition("LTh17").type ||
        pCell->type == get_cell_definition("LTh1").type ||
        pCell->type == get_cell_definition("LTh2").type)
    {
        if (r < 1463)
        {
            pCell->velocity[0] += push_strength * (x / r);
            pCell->velocity[1] += push_strength * (y / r);
        }

        // Si cellule immunitaire dépasse la barrière → retour vers le centre
        if (r > fls_inner_radius)
        {
            pCell->velocity[0] -= 1.0 * (x / r);
            pCell->velocity[1] -= 1.0 * (y / r);
        }

         // Détection de proximité d'une FLS pour déclencher le déplacement global
        std::vector<Cell*> nearby = pCell->cells_in_my_container();
        for (Cell* neighbor : nearby)
        {
            if (neighbor == pCell)
                continue;

            if (neighbor->type == get_cell_definition("FLS").type)
            {
                double dx = neighbor->position[0] - x;
                double dy = neighbor->position[1] - y;
                double dist = sqrt(dx * dx + dy * dy);

                if (dist < 25.0)
                {
                    fls_trigger_expand = true;
                }
            }
        }
    }

    // 🔁 Étendre toutes les FLS vers l'extérieur si déclencheur activé
    if (fls_trigger_expand)
    {
        for (Cell* c : *all_cells)
        {
            if (c->type == get_cell_definition("FLS").type)
            {
                double rx = c->position[0];
                double ry = c->position[1];
                double rr = sqrt(rx * rx + ry * ry);

                double dx = rx / rr;
                double dy = ry / rr;
                c->velocity[0] += 0.1 * dx;
                c->velocity[1] += 0.1 * dy;
            }
        }
        fls_inner_radius += 1.0;
        fls_outer_radius += 1.0;
        fls_trigger_expand = false; // reset le déclencheur pour la frame suivante
    

    }
}

void final_square_function( Cell* pCell, Phenotype& phenotype , double dt )
{
	if (parameters.bools("final_square_bc") == false)
        return;

    double x = pCell->position[0];
    double y = pCell->position[1];
    double r = safe_radial_distance(x, y);

    double push_strength = 1.0; // tu peux ajuster ce paramètre

    if (pCell->type == get_cell_definition("Osteoclast").type)
    {
		if (y > 1020.0)
		{
            pCell->velocity[1] -= push_strength;
		}
		if (y < -1020.0)
		{
            pCell->velocity[1] += push_strength;
		}
		
        if (x > 0.0)
        { 
			if (x < 450.0)
			{
            // Applique une force vers l'exterieur
            pCell->velocity[0] += push_strength;
			}
			
        }
		if (x < 0.0)
        { 
			if (x > -450.0)
			{
            pCell->velocity[0] -= push_strength;
			}
        }
    }
	
	if (pCell->type == get_cell_definition("Chondrocytes").type ||
        pCell->type == get_cell_definition("Collagen").type)
	
	{
		// global 
        if (r > 2020.0)
        {
            
            pCell->velocity[0] -= push_strength * (x / r);
            pCell->velocity[1] -= push_strength * (y / r);
        }
		
		// side
		
        if (y >= -1015.0 && y <= 1015.0)
			{
				if (x > 450.0)
				{
				// Applique une force vers l'exterieur
				pCell->velocity[0] -= push_strength;
				}
			}
        
		
		if (y >= -1015.0 && y <= 1015.0)
			{
				if (x < -450.0)
				{
				pCell->velocity[0] += push_strength;
				}
			}
        
		
		
		//top
		if (x < -468.0)
		{
			if (y > 1015.0 && y < 1020)
			{
				pCell->velocity[1] += push_strength;
			}
		}
		
		if (x > 468.0)
		{
			if (y > 1015.0 && y < 1020)
			{
				pCell->velocity[1] += push_strength;
			}
		}
		
		//bottom
		if (x < -468.0)
		{
			if (y > -1020 && y < -1015)
			{
				pCell->velocity[1] -= push_strength;
			}
		}
		
		if (x > 468.0)
		{
			if (y > -1020 && y < -1015)
			{
				pCell->velocity[1] -= push_strength;
			}
		}
		
		
	}
	
	if (pCell->type == get_cell_definition("FLS").type||
        pCell->type == get_cell_definition("M0").type)
    {
        
		if (r < 2025.0)
        {
            // Applique une force vers l'exterieur
            pCell->velocity[0] += push_strength * (x / r);
            pCell->velocity[1] += push_strength * (y / r);
        }
		
    }
	
	if (pCell->type == get_cell_definition("Joint").type)
    {
        if (r > 2058.0)
		{
			if (r < 2060.0)
			{
				// Applique une force vers l'exterieur
				pCell->velocity[0] += push_strength * (x / r);
				pCell->velocity[1] += push_strength * (y / r);
			}
		}
		
		// gauche haut
		
		if (x > 0 && y>0)
		{
			if (x < 1468 && y <=1420 ) {
				
			{
			pCell->velocity[0] += push_strength;
			}
		}
		}
		
		if (x > 0 && y>0)
		{
			if (x >= 1468 && y <1020 ) {
				
			{
			pCell->velocity[1] += push_strength;
			}
		}
		}
		
		// gauche bas
		if (x > 0 && y<0)
		{
			if (x < 1468 && y >=-1420 ) {
				
			{
			pCell->velocity[0] += push_strength;
			}
		}
		}
		
		if (x > 0 && y<0)
		{
			if (x >= 1468 && y >-1020 ) {
				
			{
			pCell->velocity[1] -= push_strength;
			}
		}
		}
		
		// droite haut
		
		if (x < 0 && y>0)
		{
			if (x > -1468 && y <=1420 ) {
				
			{
			pCell->velocity[0] -= push_strength;
			}
		}
		}
		
		if (x < 0 && y>0)
		{
			if (x <= -1468 && y <1020 ) {
				
			{
			pCell->velocity[1] += push_strength;
			}
		}
		}
		
		// droite bas
		if (x < 0 && y<0)
		{
			if (x > -1468 && y >=-1420 ) {
				
			{
			pCell->velocity[0] -= push_strength;
			}
		}
		}
		
		if (x < 0 && y<0)
		{
			if (x <= -1468 && y >-1020 ) {
				
			{
			pCell->velocity[1] -= push_strength;
			}
		}
		}
    }
	
    if (
        pCell->type == get_cell_definition("M1").type ||
        pCell->type == get_cell_definition("M2").type ||
        pCell->type == get_cell_definition("LTCD4+").type ||
        pCell->type == get_cell_definition("LB").type ||
        pCell->type == get_cell_definition("Plasmocytes").type ||
        pCell->type == get_cell_definition("LTCD8+").type ||
        pCell->type == get_cell_definition("LTreg").type ||
        pCell->type == get_cell_definition("LTh17").type ||
        pCell->type == get_cell_definition("LTh1").type ||
        pCell->type == get_cell_definition("LTh2").type ||
        pCell->type == get_cell_definition("MoDC").type ||
        pCell->type == get_cell_definition("Monocytes").type ||
        pCell->type == get_cell_definition("Neutrophils").type)
    {
        
		if (r > 2020.0)
        {
            // Applique une force vers le centre du grand cercle
            pCell->velocity[0] -= push_strength * (x / r);
            pCell->velocity[1] -= push_strength * (y / r);
        }
		
		//centre
		if (y >= -1420.0 && y <= 1420.0 && x > 50.0 )
		{
			pCell->velocity[0] -= push_strength;
		}
		if (y >= -1420.0 && y <= 1420.0 && x < -50.0 )
		{
			pCell->velocity[0] += push_strength;
		}
		
		//haut à gauche
		if (x >= -1468.0 && x <= -50.0 && y >0 )
		{ 
			if (y < 1420)
			{
			pCell->velocity[1] += push_strength;
			}
		}
		
		// haut à droite
		if (x <= 1468.0 && x >= 50.0 && y >0 )
		{ 
			if (y < 1420)
			{
			pCell->velocity[1] += push_strength;
			}
		}
	
		// bas à droite
		if (x <= 1468.0 && x >= 50.0 && y <0 )
		{ 
			if (y > -1420)
			{
			pCell->velocity[1] -= push_strength;
			}
		}
		
		// bas à gauche
		if (x <= -1468.0 && x >= -50.0 && y <0 )
		{ 
			if (y > -1420)
			{
			pCell->velocity[1] -= push_strength;
			}
		}
		
    }
} 


 
void contact_function(Cell* pMe, Phenotype& phenoMe, Cell* pOther, Phenotype& phenoOther, double dt)
{return;}



void treatment_function (mpi_Environment& world)
{
	if (PhysiCell::parameters.bools.find_index("treatment") != -1) 
	{
		int treatment_substrate_index = BioFVM::microenvironment.find_density_index(PhysiCell::parameters.strings("treatment_substrate"));

		if (PhysiCell::parameters.bools("treatment")){
		
			if (
				(((int)PhysiCell::PhysiCell_globals.current_time) % PhysiCell::parameters.ints("treatment_period")) == 0 
				&& !BioFVM::microenvironment.get_substrate_dirichlet_activation(treatment_substrate_index)
			)
			{
				if (IOProcessor(world)) std::cout << PhysiCell::parameters.strings("treatment_substrate") << " activation at t=" << PhysiCell::PhysiCell_globals.current_time << std::endl;
				BioFVM::microenvironment.set_substrate_dirichlet_activation(treatment_substrate_index, true);	
			}

			if (
				(((int)PhysiCell::PhysiCell_globals.current_time) % PhysiCell::parameters.ints("treatment_period")) == PhysiCell::parameters.ints("treatment_duration") 
				&& BioFVM::microenvironment.get_substrate_dirichlet_activation(treatment_substrate_index)
			)
			{
				if (IOProcessor(world)) std::cout << PhysiCell::parameters.strings("treatment_substrate") << " inactivation at t=" << PhysiCell::PhysiCell_globals.current_time << std::endl;
				BioFVM::microenvironment.set_substrate_dirichlet_activation(treatment_substrate_index, false);	
			}
			
		} else if ( BioFVM::microenvironment.get_substrate_dirichlet_activation(treatment_substrate_index) ){
			if (IOProcessor(world)) std::cout << PhysiCell::parameters.strings("treatment_substrate") << " inactivation (NO TREATMENT) at t=" << PhysiCell::PhysiCell_globals.current_time << std::endl;
			BioFVM::microenvironment.set_substrate_dirichlet_activation(treatment_substrate_index, false);	
		}
	}
}
