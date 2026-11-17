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

/*=======================================*/
/* Include mpi.h in the parallel version */
/*=======================================*/

#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <filesystem>
#include <string> 

/*=============================*/
/* Include ./DistPhy/DistPhy.h */
/*=============================*/

#include "../../../modules/PhysiCell_standard_modules.h" 

#include "./custom_modules/heterogeneity.h"


using namespace BioFVM;
using namespace PhysiCell;

/*============================*/
/* Use DistPhy::mpi namespace */
/*============================*/

using namespace DistPhy::mpi; 


int main( int argc, char* argv[] )
{
    
/*=======================================================================================*/
/* Create mpi_Environment object, initialize it, then create Cartesian Topology          */
/*=======================================================================================*/
	
	mpi_Environment world;                         //object contains size of communicator, rank of process
    world.Initialize();  //OK                          //Initialize using MPI_THREAD_FUNNELED, find comm. size and comm. rank
    mpi_Cartesian cart_topo;                       //Contains dims[3], coords[3] array and MPI_Comm mpi_cart_comm
    cart_topo.Build_Cartesian_Topology(world);     //Create 1-D X decomposition by setting dims[1]=size. 
    cart_topo.Find_Cartesian_Coordinates(world);   //Find Cartesian Topology coordinates of each process
    cart_topo.Find_Left_Right_Neighbours(world); 	 //Finds left/right immediate neighbour processes ranks and stores in X_LEFT/X_RIGHT
    
    bool XML_status = false;
    
/*====================================================================================*/    
/* Call parallel version of load_ function, all processes must load the complete file */
/*====================================================================================*/    

    
	const std::string config_path = argc > 1 ? argv[1] : "./config/PhysiCell_settings.xml";
	const std::string output_folder = argc > 2 ? argv[2] : "";
	const bool validate = argc > 3 && std::string( argv[3] ) == "--validate";
	if( argc > 1 )
	{ 
	        XML_status = load_PhysiCell_config_file( config_path, world ); //Jose: por aqui
        
  }
	else
	{ 
        XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml", world ); 
        
  }
    
	if( !XML_status )
	{ 
        exit(-1); 
  }

	if( !output_folder.empty() )
	{
		PhysiCell_settings.folder = output_folder;
	}
	if( validate )
	{
		PhysiCell_settings.omp_num_threads = 1;
		default_microenvironment_options.X_range = { -80.0, 80.0 };
		default_microenvironment_options.Y_range = { -60.0, 60.0 };
		default_microenvironment_options.Z_range = { -60.0, 60.0 };
		default_microenvironment_options.dx = 20.0;
		default_microenvironment_options.dy = 20.0;
		default_microenvironment_options.dz = 20.0;
		parameters.doubles( "tumor_radius" ) = 30.0;
		parameters.doubles( "oncoprotein_mean" ) = 1.0;
		parameters.doubles( "oncoprotein_sd" ) = 0.0;
		parameters.doubles( "oncoprotein_min" ) = 1.0;
		parameters.doubles( "oncoprotein_max" ) = 1.0;
	}

	if( world.rank == 0 )
	{
		std::filesystem::create_directories( PhysiCell_settings.folder );
	}
	MPI_Barrier( cart_topo.mpi_cart_comm );

	//Setting OpenMP threads, this takes precedence over OMP_NUM_THREADS so be careful
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);
	
	//PNRG setup 
	SeedRandom(); 
	
	//Time units setup 
	std::string time_units = "min"; 

/*================================================================================================*/
/* Objects of both mpi_Environment and mpi_Cartesian are needed in setup microenvironment(...) 		*/
/* These objects are passed to initialize_microenvironment(...) which calls resize_uniform(...)   */
/*================================================================================================*/
	
	setup_microenvironment(world, cart_topo); //Jose: por aqui

	//PhysiCell setup
 	
	//Set mechanics voxel size, must be >= Diffusion voxel size
	double mechanics_voxel_size = 20;
    
/*=========================================================*/
/* Calling the parallel version of Cell Container creation */
/*=========================================================*/
    
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size, world, cart_topo);
	
	create_cell_types( world, cart_topo );
    
	setup_tissue(microenvironment, world, cart_topo);      //Custom parallel function 

	int validation_failed = 0;
	if( validate )
	{
		const int local_cells = static_cast<int>( all_cells->size() );
		int global_cells = 0;
		MPI_Allreduce( &local_cells, &global_cells, 1, MPI_INT, MPI_SUM, cart_topo.mpi_cart_comm );
		const int expected_cells = static_cast<int>( create_cell_sphere_positions(
			cell_defaults.phenotype.geometry.radius, parameters.doubles( "tumor_radius" ) ).size() );
		validation_failed = global_cells == expected_cells && global_cells > 0 ? 0 : 1;
		if( world.rank == 0 )
		{
			std::cout << ( validation_failed == 0 ? "PASS" : "FAIL" )
				<< ": distributed initial population count=" << global_cells
				<< ", expected=" << expected_cells << std::endl;
		}
	}
	
	/* Users typically start modifying here. START USERMODS */ 
	
	/*......................................................*/
	
	/* Users typically stop modifying here. END USERMODS 		*/ 
	
	//Set MultiCellDS save options, these MUST all be true 
	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true ); 
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true );
	
	//Save a simulation snapshot 
	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() ); 
	
	//Use the parallel version of the function for XML file writing
	//save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time, world, cart_topo );  //Jose: por aqui
	save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time, world, cart_topo );

	
	//Gracefully shut-down MPI 
  	world.Finalize(); 

	return validation_failed;
}
