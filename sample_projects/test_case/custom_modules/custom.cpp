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
#include "../modules/PhysiCell_settings.h"
#include "../DistPhy/DistPhy_Utils.h"
#include "../DistPhy/DistPhy_Collective.h"
#include <sstream>

using namespace DistPhy::mpi;

void create_cell_types( void )
{
	// Use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints( "random_seed" ) ); 

	// Housekeeping 
	initialize_default_cell_definition();
	//same
	cell_defaults.parameters.o2_proliferation_saturation = 38.0;  
	cell_defaults.parameters.o2_reference = 38.0; 


	cell_defaults.functions.update_phenotype = NULL;
	cell_defaults.functions.update_phenotype_parallel = dummy_phenotype_update;
  
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;	
	// Turn the default cycle model to live, 
	// so it's easier to turn-off proliferation
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	//thsi we keep because they are not default.
	//Now lets chek physicell
	// commented the lines below because of consistency with vanilla version
	// cell_defaults.phenotype.geometry.polarity = 1.0; 
	// cell_defaults.phenotype.motility.restrict_to_2D = true; 
	// std::cout<<"After "<<cell_defaults.phenotype.geometry.polarity <<cell_defaults.phenotype.motility.restrict_to_2D<<std::endl;
	//These were not set in the vanilla version 
	// int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	// int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	// int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	
	
	// cell_defaults.functions.set_orientation = up_orientation;



	
	// Use default proliferation and death 
	



	initialize_cell_definitions_from_pugixml(); 
	
	

	// Set the default cell type to no phenotype updates 
	
	
	// cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_oncoprotein; 
 
	// cell_defaults.name = "cancer cell"; 
	// cell_defaults.type = 0; 
	
	// // Add custom data 
	
	cell_defaults.custom_data.add_variable( "oncoprotein" , "dimensionless", 1.0 ); 
	
	build_cell_definitions_maps(); //Uncommenting it
	
	// display_cell_definitions( world, cart_topo );
	//  <------ will print using all processes, thus disabled
	
	/*---------------------------------------------------------------------------------------*/
	/* display_cell_definitions(...) has been disabled above as the printing will be done by */
	/* ALL processes. If you want to print (for checking, testing etc.) then create a new 	 */
	/* version of void create_cell_types(mpi_Environment &world, mpi_Cartesian &cart_topo)	 */
	/* and then use the IOProcessor(world) function to print using ONLY the root process  	 */
	/*---------------------------------------------------------------------------------------*/

	return; 
}



/*==============================================*/
/* Parallel version of setup_microenvironment() */
/*==============================================*/

void setup_microenvironment(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	
	//PhysiCell-X ONLY supports 3-D problems, hence MUST set 3-D option
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		if(IOProcessor(world))
            std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl; 
		default_microenvironment_options.simulate_2D = false; 
	}
			
	initialize_microenvironment(world, cart_topo); 	

	return; 
}	


/*-------------------------------------------------------------------*/
/* Miguel Ponce-de-Leon's function for generating positions of cells */
/*-------------------------------------------------------------------*/

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
  	double y_spacing= cell_radius*2;
  	double z_spacing= cell_radius*sqrt(3);
	
	std::vector<double> tempPoint(3,0.0);

	for(double z=-sphere_radius;z<sphere_radius;z+=z_spacing, zc++)
	{
		for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
		{
			for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
				tempPoint[1]=y + (xc%2) * cell_radius;
				tempPoint[2]=z;

				if(sqrt(norm_squared(tempPoint))< sphere_radius)
				{ cells.push_back(tempPoint); }
			}

		}
	}
	return cells;

}
/*-------------------------------------------------------------------*/
/* Thalia Diniaco function for generating cell position by reading a csv */
/*-------------------------------------------------------------------*/
std::vector<std::vector<double>> read_cell_positions(std::string filename, char delimiter, bool header)
{
	std::fstream fin;
	std::vector<std::vector<double>> positions;
	
	fin.open(filename, std::ios::in);

	// Read the Data from the file
	// as String Vector
	std::vector<std::string> row;
	std::string line, word;

	if (header)
	{ getline(fin, line); }

	do
	{
		row.clear();

		// read an entire row and
		// store it in a string variable 'line'
		getline(fin, line);

		// used for breaking words
		std::stringstream s(line);

		while (getline(s, word, delimiter))
		{ 
			row.push_back(word); 
		}

		std::vector<double> tempPoint(3,0.0);
		tempPoint[0]= std::stof(row[0]);
		tempPoint[1]= std::stof(row[1]);
		tempPoint[2]= std::stof(row[2]);

		positions.push_back(tempPoint);
	} while (!fin.eof());

	return positions;

}
/*------------------------------------------------------------------------*/
/* Parallel version of setup_tissue(), replacing this function completely */
/* by Miguel's version of setup_tissue and then parallelizing             */
/*------------------------------------------------------------------------*/

void setup_tissue(Microenvironment &m, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
  // Place a cluster of tumor cells at the center
   
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	

	Cell* pCell = NULL; 
    
    std::vector<std::vector<double>> positions;		 
    std::vector<std::vector<double>> generated_positions_at_root;
    
    /*----------------------------------------------------------------------------------------------------*/
    /* Object of mpi_CellPositions must be declared for all processes because distribute_cell_positions() */
    /* function will pass 2 objects of the kind mpi_CellPositions and mpi_MyCells                         */
    /*----------------------------------------------------------------------------------------------------*/
    
    mpi_CellPositions cp;                //To store cell positions, cell IDs, no. of cell IDs at root only (for all processes)
    mpi_MyCells       mc;                //To store cell positions, cell IDs, no. of cells at each process.
	
    if(world.rank == 0) //Only the MPI Rank 0 process will generate positions
    {
		if ( parameters.bools("read_init") )
		{
			std::string csv_fname = parameters.strings("init_cells_filename");
			std::cout << "init file " << parameters.strings("init_cells_filename") << std::endl; 
			generated_positions_at_root = read_cell_positions(csv_fname, ',', false);
			std::cout << "finished reading init" << parameters.strings("init_cells_filename") << std::endl; 

		}
		else
		{
			double tumor_radius = parameters.doubles( "tumor_radius" ); // 250.0; 
			generated_positions_at_root = create_cell_sphere_positions(cell_radius,tumor_radius);

		}

        
        int strt_cell_ID = Basic_Agent::get_max_ID_in_parallel();                               //IDs for new cells (positions) will start from the current highest ID
        
        
        cp.positions_to_rank_list(generated_positions_at_root, 
                                  m.mesh.bounding_box[0], m.mesh.bounding_box[3], m.mesh.bounding_box[1], m.mesh.bounding_box[4], m.mesh.bounding_box[2], m.mesh.bounding_box[5], 
                                  m.mesh.dx, m.mesh.dy, m.mesh.dz, 
                                  world, cart_topo, strt_cell_ID);
        
        Basic_Agent::set_max_ID_in_parallel(strt_cell_ID + generated_positions_at_root.size()); //Highest ID now is the starting ID + no. of generated coordinates ! 
    }
    
    distribute_cell_positions(cp, mc, world, cart_topo);                                        //Distribute cell positions to individual processes
    if(IOProcessor(world))
        std::cout << "creating " << generated_positions_at_root.size() << " closely-packed tumor cells ... " << std::endl;

	double x = 0.0; 
	// double x_outer = tumor_radius; 
	double y = 0.0; 
	
	double p_mean = parameters.doubles( "oncoprotein_mean" ); 
	double p_sd 	= parameters.doubles( "oncoprotein_sd" ); 
	double p_min 	= parameters.doubles( "oncoprotein_min" ); 
	double p_max 	= parameters.doubles( "oncoprotein_max" ); 
	

		for( int i=0; i < mc.my_no_of_cell_IDs; i++ )
		{
		  
			pCell = create_cell(mc.my_cell_IDs[i]); // tumor cell --> This has to be replaced by create_cell(mc.my_cell_IDs[i])
	  	  		
			pCell->assign_position(mc.my_cell_coords[3*i],mc.my_cell_coords[3*i+1],mc.my_cell_coords[3*i+2],world, cart_topo); //pCell->assign_position( positions[i] );
				 	
		 	pCell->custom_data[0] = NormalRandom( p_mean, p_sd );
		

			if( pCell->custom_data[0] < p_min )
			{ 
				pCell->custom_data[0] = p_min; 
			}
			if( pCell->custom_data[0] > p_max )
			{ 
				pCell->custom_data[0] = p_max; 
			}
		} 
 
	
	
	double local_sum = 0.0, global_sum = 0.0; 
	double local_min = 9e9, global_min = 0.0; 
	double local_max = -9e9, global_max = 0.0; 
	int local_cells, global_cells; 

/*------------------------------------------------------------------------------------------------------*/
/* The global_min/max are only for display purposes. We can just calculate local_min then do MPI_Reduce */
/* at the root process to display it. Since the mean is needed by all processes to calculate squared 	  */
/* sum of differences, I should do MPI_Allreduce(). Right now everything is globally distributed			  */
/* later we can decide to use MPI_Reduce selectively. 																									*/
/*------------------------------------------------------------------------------------------------------*/	
	
	for( int i=0; i < all_cells->size() ; i++ )
	{
		double r = (*all_cells)[i]->custom_data[0]; 
		local_sum += r;
		if( r < local_min )
		{ 
			local_min = r; 
		} 
		if( r > local_max )
		{ 
			local_max = r; 
		}
	}
	
	local_cells 	= mc.my_no_of_cell_IDs; 											//or all_cells->size()
	
	global_sum 		= distribute_global_sum(local_sum, cart_topo);
	global_cells 	= distribute_global_sum(local_cells, cart_topo); 
	global_max		= distribute_global_max(local_max, cart_topo); 
	global_min	  = distribute_global_min(local_min, cart_topo); 
	 
	double mean = global_sum / ( global_cells + 1e-15 ); 
	
	// Compute standard deviation 
	
	local_sum = 0.0;
	 
	for( int i=0; i < all_cells->size(); i++ )
	{
		local_sum +=  ( (*all_cells)[i]->custom_data[0] - mean )*( (*all_cells)[i]->custom_data[0] - mean ); 
	}
	global_sum = distribute_global_sum(local_sum, cart_topo);
		
	double standard_deviation = sqrt( global_sum / ( global_cells - 1.0 + 1e-15 ) ); 
	
	if(IOProcessor(world))
	{
		std::cout << std::endl << "Oncoprotein summary: " << std::endl<< "===================" << std::endl; 
		std::cout << "mean: " << mean << std::endl; 
		std::cout << "standard deviation: " << standard_deviation << std::endl; 
		std::cout << "[min max]: [" << global_min << " " << global_max << "]" << std::endl << std::endl;
		std::cout << "end of setup tissue " << generated_positions_at_root.size() << " closely-packed tumor cells ... " << std::endl;
	} 
	
	return; 
}

// Custom cell phenotype function to scale immunostimulatory factor with hypoxia 

void tumor_cell_phenotype_with_oncoprotein( Cell* pCell, Phenotype& phenotype, double dt )
{
	// update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	
	// If cell is dead, don't bother with future phenotype changes. 
	
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 		
		return; 
	}

	// Multiply proliferation rate by the oncoprotein 
	
	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" ); 
	// this is stohastic
	// phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ; 
	
	return; 
}
void tumor_cell_phenotype_with_oncoprotein( Cell* pCell, Phenotype& phenotype, double dt,mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	
	// If cell is dead, don't bother with future phenotype changes. 
	
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 		
		return; 
	}

	// Multiply proliferation rate by the oncoprotein 
	
	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" ); 

	phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ; 
	
	return; 
}
void dummy_phenotype_update( Cell* pCell, Phenotype& phenotype, double dt,mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	
	
	// If cell is dead, don't bother with future phenotype changes. 
	
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 	
		pCell->functions.update_phenotype_parallel=NULL;	
		return; 
	}

	
	return; 
}

std::vector<std::string> heterogeneity_coloring_function( Cell* pCell )
{
	static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" ); 
	
	static double p_min = parameters.doubles( "oncoprotein_min" ); 
	static double p_max = parameters.doubles( "oncoprotein_max" ); 
	
	// Immune are black
	
	std::vector< std::string > output( 4, "black" ); 
	
	if( pCell->type == 1 )
	{ 
		return output; 
	} 
	
	// Live cells are green, but shaded by oncoprotein value
	 
	if( pCell->phenotype.death.dead == false )
	{
		int oncoprotein = (int) round( (1.0/(p_max-p_min)) * (pCell->custom_data[oncoprotein_i]-p_min) * 255.0 ); 
		char szTempString [128];
		sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein );
		output[0].assign( szTempString );
		output[1].assign( szTempString );

		sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/p_max) , (int)round(output[0][1]/p_max) , (int)round(output[0][2]/p_max) );
		output[2].assign( szTempString );
		
		return output; 
	}

	// If not, dead colors 
	
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )  // Apoptotic - Red
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
	}
	
	// Necrotic - Brown
	
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
	}	
	
	return output; 
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
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic)
		{ local_count++; }
	}
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}


int total_dead_but_ghost(mpi_Environment &world, mpi_Cartesian &cart_topo)
{	
	int local_count = 0;
	#pragma omp parallel for reduction (+: local_count )
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell* pCell = (*all_cells)[i];
		if( pCell->phenotype.death.dead==false )
		{ continue; }
		if( pCell->phenotype.cycle.current_phase().code != PhysiCell_constants::apoptotic &&
			pCell->phenotype.cycle.current_phase().code != PhysiCell_constants::necrotic_swelling && 
			pCell->phenotype.cycle.current_phase().code != PhysiCell_constants::necrotic_lysed && 
			pCell->phenotype.cycle.current_phase().code != PhysiCell_constants::necrotic )
		{ local_count++; 
			// std::cout<<pCell->phenotype.cycle.current_phase().code<<" Line written by Rank:"<<world.rank<<std::endl;
		}
	}
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}


int total_apoptosis_not_codealive(mpi_Environment &world, mpi_Cartesian &cart_topo)
{	
	int local_count = 0;
	#pragma omp parallel for reduction (+: local_count )
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell* pCell = (*all_cells)[i];
		if( pCell->phenotype.death.dead==false )
		{ continue; }
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic && pCell->phenotype.cycle.current_phase().code != 14 )
		{ local_count++; 
			// std::cout<<pCell->phenotype.cycle.current_phase().code<<" Line written by Rank:"<<world.rank<<std::endl;
		}
	}
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}

int total_apoptosis_and_codealive(mpi_Environment &world, mpi_Cartesian &cart_topo)
{	
	int local_count = 0;
	#pragma omp parallel for reduction (+: local_count )
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell* pCell = (*all_cells)[i];
		if( pCell->phenotype.death.dead==false )
		{ continue; }
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic && pCell->phenotype.cycle.current_phase().code == 14 )
		{ local_count++; 
			// std::cout<<pCell->phenotype.cycle.current_phase().code<<" Line written by Rank:"<<world.rank<<std::endl;
		}
	}
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;
}
int total_necrotic_ghost(mpi_Environment &world, mpi_Cartesian &cart_topo){
	int local_count = 0;
	#pragma omp parallel for reduction (+: local_count )
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell* pCell = (*all_cells)[i];
		if( pCell->phenotype.death.dead==false )
		{ continue; }
		if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling && 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed && 
			pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic)
			{
				if(pCell->phenotype.cycle.current_phase().code == 14)
					{local_count++; }
			}
		
	}
	int global_count;
	MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	return global_count;


}