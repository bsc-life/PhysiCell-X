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

#include "custom.h"
#include "../BioFVM/BioFVM.h"  
#include "../DistPhy/DistPhy_Utils.h"
#include "../DistPhy/DistPhy_Collective.h"
using namespace DistPhy::mpi;

// declare cell definitions here 

std::vector<bool> nodes;

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = NULL;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL;
	cell_defaults.functions.update_phenotype_parallel = tumor_cell_phenotype_with_signaling; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	cell_defaults.custom_data.add_variable(parameters.strings("node_to_visualize"), "dimensionless", 0.0 ); //for paraview visualization

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/

	build_cell_definitions_maps(); 
	
	display_cell_definitions( std::cout ); 
	
	return; 
}

/*==============================================*/
/* Parallel version of setup_microenvironment() */
/*==============================================*/

void setup_microenvironment(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		if(IOProcessor(world))
            std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = false; 
	}
			
	initialize_microenvironment(world, cart_topo); 	

	return; 
}	


/*------------------------------------------------------------------------*/
/* Parallel version of setup_tissue(), replacing this function completely */
/* by Miguel's version of setup_tissue and then parallelizing             */
/*------------------------------------------------------------------------*/
/*
void setup_tissue(Microenvironment &m, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
    
	Cell* pCell = NULL; 
    
    std::vector<std::vector<double>> positions;		//What is this variable for ?  
    std::vector<std::vector<double>> generated_positions_at_root;
    
    mpi_CellPositions cp;                //To store cell positions, cell IDs, no. of cell IDs at root only (for all processes)
    mpi_MyCells       mc;                //To store cell positions, cell IDs, no. of cells at each process.
	
   
	for (auto t_cell_definition: cell_definitions_by_index) {

		if(world.rank == 0) //Only the MPI Rank 0 process will generate positions
		{
			generated_positions_at_root.clear();
			std::vector<double> t_position;
			for (int i=0; i < 90; i+= 10)
				for (int j=0; j < 90; j+= 10){
					
					t_position.clear();
					switch (t_cell_definition->type) {
						case 0://get_cell_definition("default").type:	
							t_position.push_back(-i-10);
							t_position.push_back(-j-10);
							break;
						case 1://get_cell_definition("other").type:	
							t_position.push_back(+i+10);
							t_position.push_back(-j-10);
							break;
						case 2://get_cell_definition("another").type:	
							t_position.push_back(-i-10);
							t_position.push_back(+j+10);
							break;
						case 3://get_cell_definition("yet_another").type:	
							t_position.push_back(+i+10);
							t_position.push_back(+j+10);
							break;
						case 4://get_cell_definition("yet_yet_another").type:	
							t_position.push_back(+i+110);
							t_position.push_back(-j-10);
							break;
						case 5://get_cell_definition("last_one").type:	
							t_position.push_back(+i+110);
							t_position.push_back(+j+10);
							break;
							
							
					}
					t_position.push_back(0);
					generated_positions_at_root.push_back(t_position);
				}
				
				
			int strt_cell_ID = Basic_Agent::get_max_ID_in_parallel();                               //IDs for new cells (positions) will start from the current highest ID
        
        
			cp.positions_to_rank_list(generated_positions_at_root, 
									m.mesh.bounding_box[0], m.mesh.bounding_box[3], m.mesh.bounding_box[1], m.mesh.bounding_box[4], m.mesh.bounding_box[2], m.mesh.bounding_box[5], 
									m.mesh.dx, m.mesh.dy, m.mesh.dz, 
									world, cart_topo, strt_cell_ID);
			
			Basic_Agent::set_max_ID_in_parallel(strt_cell_ID + generated_positions_at_root.size()); //Highest ID now is the starting ID + no. of generated coordinates ! 
		}
	
		distribute_cell_positions(cp, mc, world, cart_topo);                                        //Distribute cell positions to individual processes
		
		
		if(IOProcessor(world)){
			std::cout << "creating " << generated_positions_at_root.size() << " cells ... " << std::endl;
			
		}	
	
	}
		
	for( int i=0; i < mc.my_no_of_cell_IDs; i++ )
	{	
		// Here it's simple to get the cell type, because of the ordering of the ids
		// but should it be part of the MyCells structure ?
		int cell_type = mc.my_cell_IDs[i]/81;
		Cell_Definition t_cell_definition = *(cell_definitions_by_index[cell_type]);
		pCell = create_cell( t_cell_definition, mc.my_cell_IDs[i] ); 
		pCell->assign_position(mc.my_cell_coords[3*i],mc.my_cell_coords[3*i+1],mc.my_cell_coords[3*i+2],world, cart_topo); //pCell->assign_position( positions[i] );
	}
	
	return; 
}*/

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius, double x_min, double x_max, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
   std::vector<std::vector<double>> cells;
   int xc=0,yc=0,zc=0;
   double x_spacing= cell_radius*sqrt(3);
   double y_spacing= cell_radius*2;
   double z_spacing= cell_radius*sqrt(3);
	
	std::vector<double> tempPoint(3,0.0);

	for(double z=-sphere_radius;z<sphere_radius;z+=z_spacing, zc++)
	{
		for(double x=x_min;x<x_max;x+=x_spacing, xc++)
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

void setup_tissue(Microenvironment &m, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
    
	Cell* pCell = NULL; 
    
    std::vector<std::vector<double>> positions;		//What is this variable for ?  
    std::vector<std::vector<double>> generated_positions_at_root;
	double tumor_radius = parameters.doubles( "tumor_radius" );
	double cell_radius = cell_defaults.phenotype.geometry.radius;
	double cell_spacing = 0.95 * 2.0 * cell_radius;
    
    /*----------------------------------------------------------------------------------------------------*/
    /* Object of mpi_CellPositions must be declared for all processes because distribute_cell_positions() */
    /* function will pass 2 objects of the kind mpi_CellPositions and mpi_MyCells                         */
    /*----------------------------------------------------------------------------------------------------*/
    
    mpi_CellPositions cp;                //To store cell positions, cell IDs, no. of cell IDs at root only (for all processes)
    mpi_MyCells       mc;                //To store cell positions, cell IDs, no. of cells at each process.
	
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

	int c_defs = cell_definitions_by_index.size();
	int n_cells = parameters.ints("number_of_cells");
	generated_positions_at_root.clear();
	
	/*
	
	
	if (world.rank == 0) {
		for( int k=0; k < c_defs ; k++ )
		{
			Cell_Definition* pCD = cell_definitions_by_index[k]; 
			std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
			for( int n = 0 ; n < n_cells ; n++ )
			{
				std::vector<double> position = {0,0,0}; 
				position[0] = Xmin + UniformRandom()*Xrange; 
				position[1] = Ymin + UniformRandom()*Yrange; 
				position[2] = Zmin + UniformRandom()*Zrange; 
				
				generated_positions_at_root.push_back(position);
			}
		}
		int strt_cell_ID = Basic_Agent::get_max_ID_in_parallel(); 
		cp.positions_to_rank_list(generated_positions_at_root, 
							m.mesh.bounding_box[0], m.mesh.bounding_box[3], m.mesh.bounding_box[1], m.mesh.bounding_box[4], m.mesh.bounding_box[2], m.mesh.bounding_box[5], 
							m.mesh.dx, m.mesh.dy, m.mesh.dz, 
							world, cart_topo, strt_cell_ID);
	
		Basic_Agent::set_max_ID_in_parallel(strt_cell_ID + generated_positions_at_root.size()); //Highest ID now is the starting ID + no. of generated coordinates ! 
	}
	
	distribute_cell_positions(cp, mc, world, cart_topo); */                                       //Distribute cell positions to individual processes

		double x = 0.0; 
		double x_outer = tumor_radius; 
		double y = 0.0; 

		double p_mean = parameters.doubles( "oncoprotein_mean" ); 
		double p_sd = parameters.doubles( "oncoprotein_sd" ); 
		double p_min = parameters.doubles( "oncoprotein_min" ); 
		double p_max = parameters.doubles( "oncoprotein_max" ); 
		Cell_Definition* pCD = cell_definitions_by_index[0]; //default
		

		int x_min = m.mesh.x_coordinates[0];
		int x_max = m.mesh.x_coordinates[m.mesh.x_coordinates.size() -1 ];
		vector<vector<double>> local_positions = create_cell_sphere_positions(cell_radius,tumor_radius, x_min, x_max ,world, cart_topo);
		unsigned long int local_cells = local_positions.size();

		// Allocate space for every rank’s count
		unsigned long int* remote_cells = new unsigned long int[world.size];

		// All ranks gather local_cells into remote_cells
		MPI_Allgather(&local_cells, 1, MPI_UNSIGNED_LONG,
					remote_cells,   1, MPI_UNSIGNED_LONG,
					cart_topo.mpi_cart_comm);
		
		if (world.rank == 0) {
			unsigned long int count = 0;
			for (int i = 0; i < world.size; i++) count += remote_cells[i];
			std::cout << "Initial distribution of " << count << " cells:" << std::endl;
		}
			 
		
		int start_id = 0;
		for (int i = 0; i < world.rank; i++) start_id += remote_cells[i];
	
		std::cout << "\t[Rank" << world.rank <<"]: " << local_positions.size() << std::endl;

		for( int i=0; i < local_positions.size(); i++ )
		{
		  
			pCell = create_cell(*pCD, start_id + i); // tumor cell --> This has to be replaced by create_cell(mc.my_cell_IDs[i])
	  	  		
			pCell->assign_position(local_positions[i][0],local_positions[i][1],local_positions[i][2],world, cart_topo); //pCell->assign_position( positions[i] );
				 	
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

		
        
        

		
		
	if(IOProcessor(world)){
		std::cout << "creating " << generated_positions_at_root.size() << " cells ... " << std::endl;
		
	}	
	
	/*
	for( int i=0; i < mc.my_no_of_cell_IDs; i++ )
	{	
		// Here it's simple to get the cell type, because of the ordering of the ids
		// but should it be part of the MyCells structure ?
		int cell_type = mc.my_cell_IDs[i]/81;
		Cell_Definition t_cell_definition = *(cell_definitions_by_index[cell_type]);
		pCell = create_cell( t_cell_definition, mc.my_cell_IDs[i] ); 
		pCell->assign_position(mc.my_cell_coords[3*i],mc.my_cell_coords[3*i+1],mc.my_cell_coords[3*i+2],world, cart_topo); //pCell->assign_position( positions[i] );
	}

	*/
	
	return; 
}


void tumor_cell_phenotype_with_signaling( Cell* pCell, Phenotype& phenotype, double dt, mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	// std::cout << "Updating tumor cell phenotype in parallel" << std::endl;
	if (pCell->phenotype.intracellular->need_update())
	{	
		if (
			pCell->type == get_cell_definition("last_one").type
			&& PhysiCell::PhysiCell_globals.current_time >= 100.0 
			&& pCell->phenotype.intracellular->get_parameter_value("$time_scale") == 0.0
		){
			pCell->phenotype.intracellular->set_parameter_value("$time_scale", 1);
			pCell->assign_position(pCell->position[0]-300.0, pCell->position[1], 0.0, world, cart_topo);
		}
		// set_input_nodes(pCell);

		pCell->phenotype.intracellular->update();
		
		// from_nodes_to_cell(pCell, phenotype, dt);
		color_node(pCell);
	}	
}


void set_input_nodes(Cell* pCell) {}

void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt) {}


std::vector<std::string> my_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4 , "rgb(0,0,0)" );
	
	if ( !pCell->phenotype.intracellular->get_boolean_variable_value( parameters.strings("node_to_visualize") ) )
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
		
	}
	else{
		output[0] = "rgb(0, 255,0)";
		output[2] = "rgb(0, 125,0)";
	}
	
	return output;
}

void color_node(Cell* pCell){
	std::string node_name = parameters.strings("node_to_visualize");
	pCell->custom_data[node_name] = pCell->phenotype.intracellular->get_boolean_variable_value(node_name);
}