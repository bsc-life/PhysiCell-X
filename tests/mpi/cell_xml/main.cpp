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
#include <iomanip>
#include <sstream>
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

namespace
{

unsigned long long definition_fingerprint( void )
{
	std::ostringstream stream;
	stream << std::setprecision( 17 );
	for( Cell_Definition* definition : cell_definitions_by_index )
	{
		stream << definition->name << "|" << definition->type << "|"
			<< definition->phenotype.cycle.model().code << "|";
		for( const auto& row : definition->phenotype.cycle.data.transition_rates )
		{
			for( double value : row )
			{
				stream << value << ",";
			}
		}
		for( double value : definition->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities )
		{
			stream << value << ",";
		}
		for( double value : definition->phenotype.death.rates )
		{
			stream << value << ",";
		}
		stream << ";";
	}

	unsigned long long hash = 1469598103934665603ULL;
	for( unsigned char value : stream.str() )
	{
		hash ^= value;
		hash *= 1099511628211ULL;
	}
	return hash;
}

bool close_enough( double left, double right )
{
	return std::fabs( left - right ) <= 1e-12;
}

bool validate_asymmetric_division_config( void )
{
	if( cell_definitions_by_index.size() != 3 )
	{
		return false;
	}

	Cell_Definition* stem = find_cell_definition( "stem" );
	Cell_Definition* progenitor_1 = find_cell_definition( "progenitor_1" );
	Cell_Definition* progenitor_2 = find_cell_definition( "progenitor_2" );
	if( stem == NULL || progenitor_1 == NULL || progenitor_2 == NULL )
	{
		return false;
	}

	const std::vector<double>& probabilities =
		stem->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities;
	return probabilities.size() == 3
		&& close_enough( probabilities[0], 0.0 )
		&& close_enough( probabilities[1], 1.0 )
		&& close_enough( probabilities[2], 0.0 );
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

    
	const bool validate = argc > 2 && std::string( argv[2] ) == "--validate";
	if( argc > 1 )
	{ 
        XML_status = load_PhysiCell_config_file( argv[1], world ); //Jose: por aqui
        
  }
	else
	{ 
        XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml", world ); 
        
  }
    
	if( !XML_status )
	{ 
        exit(-1); 
  	}

	initialize_cell_definitions_from_pugixml( world, cart_topo );
		build_cell_definitions_maps();

		int validation_failed = 0;
		if( validate )
		{
			validation_failed = validate_asymmetric_division_config() ? 0 : 1;

			const unsigned long long local_fingerprint = definition_fingerprint();
			unsigned long long minimum_fingerprint = 0;
			unsigned long long maximum_fingerprint = 0;
			MPI_Allreduce( &local_fingerprint, &minimum_fingerprint, 1,
				MPI_UNSIGNED_LONG_LONG, MPI_MIN, cart_topo.mpi_cart_comm );
			MPI_Allreduce( &local_fingerprint, &maximum_fingerprint, 1,
				MPI_UNSIGNED_LONG_LONG, MPI_MAX, cart_topo.mpi_cart_comm );
			if( minimum_fingerprint != maximum_fingerprint )
			{
				validation_failed = 1;
			}

			int any_validation_failed = 0;
			MPI_Allreduce( &validation_failed, &any_validation_failed, 1,
				MPI_INT, MPI_MAX, cart_topo.mpi_cart_comm );
			validation_failed = any_validation_failed;
			if( world.rank == 0 )
			{
				std::cout << ( validation_failed == 0 ? "PASS" : "FAIL" )
					<< ": XML cell definitions and asymmetric-division probabilities"
					<< ( validation_failed == 0 ? " agree across all MPI ranks" : " are invalid or rank-dependent" )
					<< std::endl;
			}
		}

	display_cell_definitions( std::cout, world, cart_topo );
	
	//Gracefully shut-down MPI 
  	world.Finalize(); 

	return validation_failed;
}
