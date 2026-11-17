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
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <iomanip>
#include <string> 

/*=============================*/
/* Include ./DistPhy/DistPhy.h */
/*=============================*/

#include "../../../modules/PhysiCell_standard_modules.h" 

#include "./custom_modules/custom.h"


using namespace BioFVM;
using namespace PhysiCell;

/*============================*/
/* Use DistPhy::mpi namespace */
/*============================*/

using namespace DistPhy::mpi; 

namespace
{

void configure_steady_state_validation_problem( void )
{
	default_microenvironment_options.X_range = { -120.0, 120.0 };
	default_microenvironment_options.Y_range = { -40.0, 40.0 };
	default_microenvironment_options.Z_range = { -40.0, 40.0 };
	default_microenvironment_options.dx = 20.0;
	default_microenvironment_options.dy = 20.0;
	default_microenvironment_options.dz = 20.0;
	default_microenvironment_options.simulate_2D = false;

	const size_t substrates = default_microenvironment_options.initial_condition_vector.size();
	default_microenvironment_options.outer_Dirichlet_conditions = true;
	default_microenvironment_options.Dirichlet_activation_vector.assign( substrates, true );
	default_microenvironment_options.Dirichlet_all.assign( substrates, false );
	default_microenvironment_options.Dirichlet_xmin.assign( substrates, true );
	default_microenvironment_options.Dirichlet_xmax.assign( substrates, true );
	default_microenvironment_options.Dirichlet_ymin.assign( substrates, false );
	default_microenvironment_options.Dirichlet_ymax.assign( substrates, false );
	default_microenvironment_options.Dirichlet_zmin.assign( substrates, false );
	default_microenvironment_options.Dirichlet_zmax.assign( substrates, false );
	default_microenvironment_options.Dirichlet_xmin_values.resize( substrates );
	default_microenvironment_options.Dirichlet_xmax_values.resize( substrates );
	for( size_t substrate = 0; substrate < substrates; ++substrate )
	{
		default_microenvironment_options.Dirichlet_xmin_values[substrate] = 1.0 + substrate;
		default_microenvironment_options.Dirichlet_xmax_values[substrate] = 1.0 + substrate;
	}
}

double expected_steady_state_value( size_t substrate )
{
	return 1.0 + substrate;
}

int validate_steady_state_solution( mpi_Environment& world, mpi_Cartesian& cart_topo, int iterations )
{
	std::fill( microenvironment.decay_rates.begin(), microenvironment.decay_rates.end(), 0.0 );
	for( unsigned int voxel = 0; voxel < microenvironment.number_of_voxels(); ++voxel )
	{
		for( unsigned int substrate = 0; substrate < microenvironment.number_of_densities(); ++substrate )
		{
			microenvironment.density_vector( voxel )[substrate] = expected_steady_state_value( substrate );
		}
	}

	for( int iteration = 0; iteration < iterations; ++iteration )
	{
		microenvironment.simulate_diffusion_decay( diffusion_dt, world, cart_topo );
	}

	double local_max_error = 0.0;
	int local_non_finite = 0;
	for( unsigned int voxel = 0; voxel < microenvironment.number_of_voxels(); ++voxel )
	{
		for( unsigned int substrate = 0; substrate < microenvironment.number_of_densities(); ++substrate )
		{
			const double value = microenvironment.density_vector( voxel )[substrate];
			local_non_finite += std::isfinite( value ) ? 0 : 1;
			local_max_error = std::max( local_max_error,
				std::fabs( value - expected_steady_state_value( substrate ) ) );
		}
	}

	double global_max_error = 0.0;
	int global_non_finite = 0;
	MPI_Allreduce( &local_max_error, &global_max_error, 1, MPI_DOUBLE, MPI_MAX,
		cart_topo.mpi_cart_comm );
	MPI_Allreduce( &local_non_finite, &global_non_finite, 1, MPI_INT, MPI_SUM,
		cart_topo.mpi_cart_comm );

	const bool passed = global_non_finite == 0 && global_max_error <= 1e-8;
	if( world.rank == 0 )
	{
		std::cout << ( passed ? "PASS" : "FAIL" )
			<< ": constant steady-state diffusion across " << world.size
			<< " MPI rank(s), max error=" << std::setprecision( 12 ) << global_max_error
			<< ", non-finite values=" << global_non_finite << std::endl;
	}
	return passed ? 0 : 1;
}

}


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
	const bool validate = argc > 2 && std::string( argv[2] ) == "--validate";
	const int iterations = validate
		? ( argc > 3 ? std::stoi( argv[3] ) : 20 )
		: ( argc > 2 ? std::stoi( argv[2] ) : 1000000 );

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

	if( validate )
	{
		configure_steady_state_validation_problem();
		PhysiCell_settings.omp_num_threads = 1;
	}

	//Setting OpenMP threads, this takes precedence over OMP_NUM_THREADS so be careful
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);
	
	//PNRG setup 
	//SeedRandom(); 
	
	//Time units setup 
	std::string time_units = "min"; 

/*================================================================================================*/
/* Objects of both mpi_Environment and mpi_Cartesian are needed in setup microenvironment(...) 		*/
/* These objects are passed to initialize_microenvironment(...) which calls resize_uniform(...)   */
/*================================================================================================*/
	
	setup_microenvironment(world, cart_topo); //Jose: por aqui
	if( validate )
	{
		const int status = validate_steady_state_solution( world, cart_topo, iterations );
		world.Finalize();
		return status;
	}

	//PhysiCell setup
 	
	//Set mechanics voxel size, must be >= Diffusion voxel size
	double mechanics_voxel_size = 20;
    
/*=========================================================*/
/* Calling the parallel version of Cell Container creation */
/*=========================================================*
    
	
	
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

  	if (world.rank == 0) {
		std::cout << "Granmurality was: " << microenvironment.granurality << std::endl;
		std::cout << "Granurality will be: " << world.size << std::endl;
	}
  	
	microenvironment.granurality = world.size; 

	//MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	
	for (int i = 0; i<iterations; ++i) {
		auto start = std::chrono::high_resolution_clock::now();
		microenvironment.simulate_diffusion_decay( diffusion_dt, world, cart_topo );
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> duracion = end - start;

		if (world.rank == 0) { 
			std::cout << "	Iteration " << i << ":  " << duracion.count() << " ms | Median bandwith: ";
			double s = duracion.count()/1000.0;
			double size = ((double)(microenvironment.thomas_i_jump * 8.0)) / (1024.0*1024.0);
			double messages = (microenvironment.granurality )* ((world.size-1)) * 2.0;
			//cout << size << " " << messages << " " << s << std::endl;
			std::cout << (size*messages) / s << " MB/s" << std::endl;
		}	
	}
	//Gracefully shut-down MPI 
	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
	//save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time, world, cart_topo );

  	world.Finalize(); 

	return 0; 
}
