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

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h" 
#include "./addons/PhysiBoSS/src/maboss_intracellular.h"	
// put custom code modules here! 
#include "./custom_modules/custom.h" 
// #include <cppkafka/utils/buffered_producer.h>
// #include <cppkafka/configuration.h>

using namespace BioFVM;
using namespace PhysiCell;
// using namespace cppkafka;

/*====================================================*/
/* Use DistPhy::mpi namespace in the parallel version */
/*====================================================*/
using namespace DistPhy::mpi;

int main( int argc, char* argv[] )
{
	// load and parse settings file(s)

/*=======================================================================================*/
/* Create mpi_Environment object, initialize it, then create Cartesian Topology          */
/*=======================================================================================*/
	
	mpi_Environment world;                         //object contains size of communicator, rank of process
  world.Initialize();                            //Initialize using MPI_THREAD_MULTIPLE, find comm. size and comm. rank
  mpi_Cartesian cart_topo;                       //Contains dims[3], ccoords[3] array and MPI_Comm mpi_cart_comm
  cart_topo.Build_Cartesian_Topology(world);     //Create 1-D X decomposition by setting dims[1]=size. 
  cart_topo.Find_Cartesian_Coordinates(world);   //Find Cartesian Topology coordinates of each process
  cart_topo.Find_Left_Right_Neighbours(world); 	 //Finds left/right immediate neighbour processes ranks and stores in X_LEFT/X_RIGHT	
	
	bool XML_status = false; 

/*=========================================================================================================*/    
/* Call parallel version of load_PhysiCell_config_file function, all processes must load the complete file */
/*=========================================================================================================*/
	
	if( argc > 1 )
	{ 
		XML_status = load_PhysiCell_config_file( argv[1], world ); 
	}
	else
	{ 
		XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml", world ); 
	}
	
	if( !XML_status )
	{ 
		exit(-1); 
	}

	// OpenMP setup, commented out for debugging
	//omp_set_num_threads(PhysiCell_settings.omp_num_threads);
	
	// PNRG setup 
	SeedRandom( parameters.ints("random_seed") ); // Or a seed can be specified here
	
	// time setup 
	std::string time_units = "min"; 

	/* Microenvironment setup */
	
/*================================================================================================*/
/* Objects of both DistPhy_Environment and DistPhy_Cartesian are needed in setup microenvironment */
/* These objects are passed to initialize_microenvironment() which calls resize_uniform()         */
/*================================================================================================*/	
	
	setup_microenvironment(world, cart_topo); // modify this in the custom code 
	
	// User parameters
	
	double time_add_tnf = parameters.ints("time_add_tnf");
	double time_put_tnf = 0;
	double duration_add_tnf = parameters.ints("duration_add_tnf");
	double time_tnf_next = 0;
	double time_remove_tnf = parameters.ints("time_remove_tnf");
	
	// change concentration units too match voxel volume
	// double concentration_tnf = parameters.doubles("concentration_tnf") * microenvironment.voxels(0).volume * 0.000001;
	
	double concentration_tnf = parameters.doubles("concentration_tnf") * 0.1;
	// radius around which the tnf pulse is injected
	double membrane_lenght = parameters.ints("membrane_length");
	// tnf density index
	static int tnf_idx = microenvironment.find_density_index("tnf");	

	// this is to emulate PhysiBoSSv1 TNF experiment
	
	/*==============================================================================*/
	/* seed_tnf = false in the next statement so the if block (compound statement) 	*/
	/* is NOT executed. If seed_tnf becomes true THEN we need to check if the 			*/
	/* the function called within the if conditions i.e. inject_density_sphere			*/
	/* is to be parallelized or not. I don't think it needs parallelization so only */
	/* changing the simulate_diffusion_decay(...) call to parallel version 					*/
	/*==============================================================================*/
	
	bool seed_tnf = false;
	// do small diffusion steps alone to initialize densities
	
	if ( seed_tnf )
	{
		inject_density_sphere(tnf_idx, concentration_tnf, membrane_lenght);
		for ( int i = 0; i < 25; i ++ )
			microenvironment.simulate_diffusion_decay( diffusion_dt, world, cart_topo );
	}
	
	/* PhysiCell setup */ 
 	
	// set mechanics voxel size, and match the data structure to BioFVM
	//GS changed it to 25 otherwise it was 30 because 500 - (-500) / 25 = 40 
	double mechanics_voxel_size = 20; 
	
/*---------------------------------------------------------*/
/* Calling the parallel version of Cell Container creation */
/*---------------------------------------------------------*/

	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size, world, cart_topo );
	
	/* Users typically start modifying here. START USERMODS */ 
	
	create_cell_types();
	
	/* Calling the parallel version of setup_tissue(...) */
	
	setup_tissue(microenvironment, world, cart_topo);

	/* Users typically stop modifying here. END USERMODS */ 
	
	// set MultiCellDS save options 

	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true ); 
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true );
	
	// save a simulation snapshot 
	
	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() ); 
	/* Use the parallel version of the function now */
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time, world, cart_topo ); 
	
	// save a quick SVG cross section through z = 0, after setting its 
	// length bar to 200 microns 

	PhysiCell_SVG_options.length_bar = 200; 

	// for simplicity, set a pathology coloring function 
	
	std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function;
	
	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, world, cart_topo );
	
	display_citations(); 
	
	// set the performance timers 

	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	
	std::ofstream report_file;
	if( PhysiCell_settings.enable_legacy_saves == true )
	{	
		sprintf( filename , "%s/simulation_report.txt" , PhysiCell_settings.folder.c_str() ); 
		
		report_file.open(filename); 	// create the data log file 
		report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;
	}
	
	// main loop
	try 
	{		
		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt )
		{
			// save data if it's time. 
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
			{
				std::cout << "Time to save" << std::endl;
				display_simulation_status( std::cout ); 
				if( PhysiCell_settings.enable_legacy_saves == true )
				{	
					//GS commented out the next statement
					//log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);
					//Count Necrotic Apoptotic Alive cells
					// Producer
					
					std::string message;
					std::string topic_name = "cells";
					double timepoint = PhysiCell_globals.current_time;
					int alive_no,necrotic_no,apoptotic_no;
					
					/* Call the parallel versions of the function now which use MPI_Reduce at rank 0 */
					alive_no 			= total_live_cell_count(world, cart_topo);
					necrotic_no 	= total_necrosis_cell_count(world, cart_topo);
					apoptotic_no 	= total_dead_cell_count(world, cart_topo);
					pid_t pid_var = getpid();
					
					// MessageBuilder builder(topic_name);
					
					if(IOProcessor(world)) //This is rank 0 which will gather no. of cell types from all processes
						message = std::to_string(pid_var) + ';' + std::to_string(timepoint) + ';' + std::to_string(alive_no) + ';' + std::to_string(apoptotic_no) + ';' + std::to_string(necrotic_no) + ';';
					
					// message = '{' + 'process_id' + std::to_string(pid_var) + ',' + 'timepoint' + std::to_string(timepoint) + ',' 
					// + 'alive' + std::to_string(alive_no) + ',' + 'apoptotic' + std::to_string(apoptotic_no) + ',' 
					// + 'necrotic' + std::to_string(necrotic_no) + '}';
					// Define the configuration structure
					// Configuration config = { { "metadata.broker.list", "localhost:9092" } };
				    // Create the producer
				    // BufferedProducer<std::string> producer(config);
					//Produce a message
					// The message that will be sent
					// std::cout << "Message to Kafka: " << message << std::endl;
					// std::string str(message);
					// builder.partition(0).payload(str);
				    // producer.add_message(builder);
				    // producer.flush();
				}
				
				if( PhysiCell_settings.enable_full_saves == true )
				{	
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
					
					/* Use parallel version of function */
					save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time, world, cart_topo ); 
					
					//GS commented out next statement
					//MaBoSSIntracellular::save( filename, *PhysiCell::all_cells ); <--- NEEDS to be parallelized
				}
				
				PhysiCell_globals.full_output_index++; 
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
			}
			
			// save SVG plot if it's time
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 0.01 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{	
					//Gaurav Saxena is debugging SVG individually	
					//sprintf( filename , "%s/snapshot%08u_RANK%u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index, world.rank ); 
					//SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function);
					
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index ); 
					SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, world, cart_topo );
					
					PhysiCell_globals.SVG_output_index++; 
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
				}
			}

			/*
			  Custom add-ons could potentially go here. 
			*/			
			if ( PhysiCell_globals.current_time >= time_put_tnf )
			{
				time_tnf_next = PhysiCell_globals.current_time + duration_add_tnf;
				time_put_tnf += time_add_tnf;
			}

			if ( PhysiCell_globals.current_time >= time_remove_tnf )
			{
				remove_density(tnf_idx);													//I think it does not need parallelization
				time_remove_tnf += PhysiCell_settings.max_time;
			}

			if ( PhysiCell_globals.current_time <= time_tnf_next )
			{
				inject_density_sphere(tnf_idx, concentration_tnf, membrane_lenght); //I think it does not need parallelization
			}

			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt, world, cart_topo );
			
			// update te TNF receptor model of each cell
			tnf_receptor_model_main( diffusion_dt );		//I think it does not need parallelization
			
			// run PhysiCell 
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time, world, cart_topo );
			
			PhysiCell_globals.current_time += diffusion_dt;
		}

		if( PhysiCell_settings.enable_legacy_saves == true )
		{		
			//GS commented out the next statement	
			//log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
			report_file.close();
		}
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}
	
	// save a final simulation snapshot 
	
	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time, world, cart_topo ); 
	
	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, world, cart_topo );

	
	// timer 
	std::cout << std::endl << "Total simulation runtime: " << std::endl; 
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 

	world.Finalize(); 
	
	return 0; 
}


