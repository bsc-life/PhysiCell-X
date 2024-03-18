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

// declare cell definitions here
void create_cell_types(void)
{
	// use the same random seed so that future experiments have the
	// same initial histogram of oncoprotein, even if threading means
	// that future division and other events are still not identical
	// for all runs

	SeedRandom(parameters.ints("random_seed")); // or specify a seed here

	// housekeeping
	// std::cout << cell_defaults.name << std::endl; <---- will print using all processes

	
	/*-----------------------------------------------------------------------------------------------------------------------------*/
	/* For parallel settings, set the update velocity in parallel function pointer - this will call the new function which has the */
	/* same name but a different prototype and implementation (overloaded function)																								 */
	/* There may be no need to set it at this stage but setting it here for correctness and completeness. 												 */
	/*-----------------------------------------------------------------------------------------------------------------------------*/
	
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	cell_defaults.functions.update_velocity_parallel = standard_update_cell_velocity;
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	// no migration_bias needed
	cell_defaults.functions.update_migration_bias = NULL;
	// No custom rule needed
	cell_defaults.functions.custom_cell_rule = NULL;

	/*--------------------------------------------------------------------------------------
		update_pc_parameters_O2_based flag controls how is the update_phenotypes
	 	if update_pc_parameters_O2_based is set to true	the model will call
		tumor_cell_phenotype_with_signaling which just call two functions 
		sequentially: 1) update_cell_and_death_parameters_O2_based and 
		2) tnf_bm_interface_main; the former updates growth and death rates 
		based on oxygen while the second is the	function that update the boolean model. 
		If the flag is false then only the tnf_bm_interface_main is invoked
	 ---------------------------------------------------------------------------------------*/
	if (parameters.bools("update_pc_parameters_O2_based"))
	{
		cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_signaling; 
	}
	else
	{
		cell_defaults.functions.update_phenotype = tnf_bm_interface_main;	
	}

	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
	cell_defaults.functions.calculate_distance_to_membrane = NULL;
	cell_defaults.functions.set_orientation = NULL;

	/* This parses the cell definitions in the XML config file. */
	
	initialize_cell_definitions_from_pugixml();			//Only serial version exists and is sufficient.

	// initialize tnf

	tnf_receptor_model_setup();											
	tnf_boolean_model_interface_setup();						
	
	//submodel_registry.display(std::cout); <------ will print using all processes

	// Needs to initialize one of the receptor state to the total receptor value
	cell_defaults.custom_data["unbound_external_TNFR"] = cell_defaults.custom_data["TNFR_receptors_per_cell"];
	cell_defaults.custom_data["bound_external_TNFR"] = 0;
	cell_defaults.custom_data["bound_internal_TNFR"] = 0;

	build_cell_definitions_maps();
	
	//display_cell_definitions(std::cout); <------ will print using all processes
	
	/*------------------------------------------------------------------------------------*/
	/* Multiple print statements have been disabled above as the printing will be done by */
	/* ALL processes. If you want to print (for checking, testing etc.) then create a new */
	/* version of void create_cell_types(mpi_Environment &world, mpi_Cartesian &cart_topo)*/
	/* and then use the IOProcessor(world) function to print using ONLY the root process  */
	/*------------------------------------------------------------------------------------*/

	return;
}

void setup_microenvironment(void)
{
	// make sure to override and go back to 3D
	if (default_microenvironment_options.simulate_2D == true)
	{
		std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl;
		default_microenvironment_options.simulate_2D = false;
	}

	// initialize BioFVM
	initialize_microenvironment();

	return;
}

/*==============================================*/
/* Parallel version of setup_microenvironment() */
/*==============================================*/

void setup_microenvironment( mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	// make sure to override and go back to 3D
	if (default_microenvironment_options.simulate_2D == true)
	{
		std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl;
		default_microenvironment_options.simulate_2D = false;
	}
	
	//initialize BioFVM	
	initialize_microenvironment(world, cart_topo); 	
	
	return; 
}

void setup_tissue(void)
{

	std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);

	for (int i = 0; i < cells.size(); i++)
	{
		Cell *pCell;
		float x = cells[i].x;
		float y = cells[i].y;
		float z = cells[i].z;
		double elapsed_time = cells[i].elapsed_time;
		//check following funcionality
		//double elapsed_time = UniformRandom();
		pCell = create_cell(get_cell_definition("default"));
		pCell->phenotype.cycle.data.elapsed_time_in_phase = elapsed_time;
		pCell->assign_position(x, y, z);

		update_monitor_variables(pCell);
	}

	return;
}

/*---------------------------------------*/
/* Parallel version of setup_tissue(...) */
/*---------------------------------------*/

void setup_tissue(Microenvironment &m, mpi_Environment &world, mpi_Cartesian &cart_topo)
{

	double Xmin = m.mesh.bounding_box[0]; 
	double Ymin = m.mesh.bounding_box[1]; 
	double Zmin = m.mesh.bounding_box[2]; 

	double Xmax = m.mesh.bounding_box[3]; 
	double Ymax = m.mesh.bounding_box[4]; 
	double Zmax = m.mesh.bounding_box[5]; 
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin;
	
	/* The following 3 are temporary variables */

	Cell* pCell;
	std::vector<std::vector<double>> generated_positions_at_root;
	std::vector<double> position = {0,0,0};												//Temporary buffer
	
	mpi_CellPositions cp;   //To store cell positions, cell IDs, no. of cell IDs at root only (for all processes)
	mpi_MyCells       mc;   //To store cell positions, cell IDs, no. of cells at each process.
	
	/* First Generate Cell positions on rank 0 by reading from file */
	if(world.rank == 0)
	{
		std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);
	
		for (int i = 0; i < cells.size(); i++)
		{ 
			position[0] = cells[i].x; 
			position[1] = cells[i].y; 
			position[2] = cells[i].z;
			//double elapsed_time = cells[i].elapsed_time; ---> This needs to be made random and set when creating cell
			generated_positions_at_root.push_back(position); 
		}
			
		/*---------------------------------------------------------------------------------------------------*/
		/* (1) Obtain the current highest ID - rank 0 knows this as rank 0 generates all IDs 				 */
		/* (2) Distribute the positions and IDs to respective processes 									 */
		/* (3) Set the maximum ID as the current maximum ID would be needed if new cells are to be generated */
		/*---------------------------------------------------------------------------------------------------*/
	
		//IDs for new cells (positions) will start from the current highest ID 
		int strt_cell_ID = Basic_Agent::get_max_ID_in_parallel();                                    
		cp.positions_to_rank_list(generated_positions_at_root, 
								m.mesh.bounding_box[0], m.mesh.bounding_box[3], 
								m.mesh.bounding_box[1], m.mesh.bounding_box[4], 
								m.mesh.bounding_box[2], m.mesh.bounding_box[5], 
								m.mesh.dx, m.mesh.dy, m.mesh.dz, 
								world, cart_topo, strt_cell_ID);
			
		Basic_Agent::set_max_ID_in_parallel(strt_cell_ID + generated_positions_at_root.size()); //Highest ID now is the starting ID + no. of generated coordinates !
	} 
	
	distribute_cell_positions(cp, mc, world, cart_topo);                           //Distribute cell positions
 
 /* Create cells at individual processes */
 
 for( int i=0; i < mc.my_no_of_cell_IDs; i++ )
	{	
 		pCell = create_cell( get_cell_definition("default"), mc.my_cell_IDs[i] );
 		pCell->phenotype.cycle.data.elapsed_time_in_phase = UniformRandom(); 
		pCell->assign_position(mc.my_cell_coords[3*i],mc.my_cell_coords[3*i+1],mc.my_cell_coords[3*i+2],world, cart_topo); //pCell->assign_position( positions[i] );
		update_monitor_variables(pCell);				//No need to parallelize
	}
	 
} 


// custom cell phenotype function to run PhysiBoSS when is needed
void tumor_cell_phenotype_with_signaling(Cell *pCell, Phenotype &phenotype, double dt)
{
	if (phenotype.death.dead == true)
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);  //No parallelization needed
	tnf_bm_interface_main(pCell, phenotype, dt);

}


std::vector<std::string> my_coloring_function(Cell *pCell)
{
	// start with live coloring
	std::vector<std::string> output = false_cell_coloring_live_dead(pCell);

	// dead cells
	if (pCell->phenotype.death.dead == false)
	{
		static int nR_EB = pCell->custom_data.find_variable_index("bound_external_TNFR");
		float activation_threshold = pCell->custom_data.find_variable_index("TNFR_activation_threshold");

		int bounded_tnf = (int)round((pCell->custom_data[nR_EB] / activation_threshold) * 255.0);
		if (bounded_tnf > 0)
		{
			char szTempString[128];
			sprintf(szTempString, "rgb(%u,%u,%u)", bounded_tnf, bounded_tnf, 255 - bounded_tnf);
			output[0].assign("black");
			output[1].assign(szTempString);
			output[2].assign("black");
			output[3].assign(szTempString);
		}
	}

	return output;
}

// Function to read init files created with PhysiBoSSv2

std::vector<init_record> read_init_file(std::string filename, char delimiter, bool header)
{
	// File pointer
	std::fstream fin;
	std::vector<init_record> result;

	// Open an existing file
	fin.open(filename, std::ios::in);

	// Read the Data from the file
	// as String Vector
	std::vector<std::string> row;
	std::string line, word;

	if (header)
		getline(fin, line);

	do
	{
		row.clear();

		// read an entire row and
		// store it in a string variable 'line'
		getline(fin, line);

		// used for breaking words
		std::stringstream s(line);

		// read every column data of a row and
		// store it in a string variable, 'word'
		while (getline(s, word, delimiter))
		{

			// add all the column data
			// of a row to a vector
			row.push_back(word);
		}

		init_record record;
		record.x = std::stof(row[2]);
		record.y = std::stof(row[3]);
		record.z = std::stof(row[4]);
		record.radius = std::stof(row[5]);
		record.phase = std::stoi(row[13]);
		record.elapsed_time = std::stod(row[14]);

		result.push_back(record);
	} while (!fin.eof());

	return result;
}


void set_input_nodes(Cell* pCell) {}

void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt) {}

void color_node(Cell* pCell)
{
	std::string node_name = parameters.strings("node_to_visualize");
	pCell->custom_data[node_name] = pCell->phenotype.intracellular->get_boolean_variable_value(node_name);
}

void inject_density_sphere(int density_index, double concentration, double membrane_lenght)
{

// Inject given concentration on the extremities only
#pragma omp parallel for
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	{
		auto current_voxel = microenvironment.voxels(n);
		std::vector<double> cent = {current_voxel.center[0], current_voxel.center[1], current_voxel.center[2]};

		if ((membrane_lenght - norm(cent)) <= 0)
			microenvironment.density_vector(n)[density_index] = concentration;
	}
}

void remove_density(int density_index)
{
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
		microenvironment.density_vector(n)[density_index] = 0;
	std::cout << "Removal done" << std::endl;
}


double total_live_cell_count()
{
	double out = 0.0;

	for (int i = 0; i < (*all_cells).size(); i++)
	{
		if ((*all_cells)[i]->phenotype.death.dead == false && (*all_cells)[i]->type == 0)
		{
			out += 1.0;
		}
	}

	return out;
}

double total_live_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	double out = 0.0, global_out;

	for (int i = 0; i < (*all_cells).size(); i++)
	{
		if ((*all_cells)[i]->phenotype.death.dead == false && (*all_cells)[i]->type == 0)
		{
			out += 1.0;
		}
	}
	
	MPI_Reduce(&out, &global_out, 1, MPI_DOUBLE, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_out;
}

double total_dead_cell_count()
{
	double out = 0.0;

	for (int i = 0; i < (*all_cells).size(); i++)
	{
		if ((*all_cells)[i]->phenotype.death.dead == true && (*all_cells)[i]->phenotype.death.current_death_model_index == 0)
		{
			out += 1.0;
		}
	}

	return out;
}

double total_dead_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	double out = 0.0, global_out;

	for (int i = 0; i < (*all_cells).size(); i++)
	{
		if ((*all_cells)[i]->phenotype.death.dead == true && (*all_cells)[i]->phenotype.death.current_death_model_index == 0)
		{
			out += 1.0;
		}
	}

	MPI_Reduce(&out, &global_out, 1, MPI_DOUBLE, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_out;
}

double total_necrosis_cell_count()
{
	double out = 0.0;

	for (int i = 0; i < (*all_cells).size(); i++)
	{
		if ((*all_cells)[i]->phenotype.death.dead == true && (*all_cells)[i]->phenotype.death.current_death_model_index == 1)
		{
			out += 1.0;
		}
	}

	return out;
}

double total_necrosis_cell_count(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	double out = 0.0, global_out;

	for (int i = 0; i < (*all_cells).size(); i++)
	{
		if ((*all_cells)[i]->phenotype.death.dead == true && (*all_cells)[i]->phenotype.death.current_death_model_index == 1)
		{
			out += 1.0;
		}
	}

	MPI_Reduce(&out, &global_out, 1, MPI_DOUBLE, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_out;
}
