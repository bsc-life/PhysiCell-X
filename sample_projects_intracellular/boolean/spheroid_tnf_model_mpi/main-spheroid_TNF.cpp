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
+																				 +
+ [1] Cite the PhysiCell-X repository by giving its URL							 +
+																				 +
+ [2] Cite BioFVM-X: 															 +
+		Saxena, Gaurav, Miguel Ponce-de-Leon, Arnau Montagud, David Vicente Dorca,   +
+		and Alfonso Valencia. "BioFVM-X: An MPI+ OpenMP 3-D Simulator for Biological + 
+		Systems." In International Conference on Computational Methods in Systems    +
+		Biology, pp. 266-279. Springer, Cham, 2021. 							 +
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

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h" 
#include "./addons/PhysiBoSS/src/maboss_intracellular.h"	
#include "./custom_modules/custom.h" 


using namespace BioFVM;
using namespace PhysiCell;

/*====================================================*/
/* Use DistPhy::mpi namespace in the parallel version */
/*====================================================*/
using namespace DistPhy::mpi;

int main( int argc, char* argv[] )
{

	/*=======================================================================================*/
	/* Create mpi_Environment object, initialize it, then create Cartesian Topology          */
	/*=======================================================================================*/
		
	mpi_Environment world;                         //Object contains size of communicator, rank of process
	world.Initialize();                            //Initialize using MPI_THREAD_FUNNELED, find comm. size and comm. rank
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

/*=======================================================================================*/ 
/* Setting the number of threads here will take precedence over the threads specified in */ 
/* PhysiCell_settings.xml and the environment variable OMP_NUM_THREADS									 */
/*=======================================================================================*/
	
	//omp_set_num_threads(PhysiCell_settings.omp_num_threads); <--- We use OMP_NUM_THREADS
	
	SeedRandom( parameters.ints("random_seed") ); // Or a seed can be specified here
	std::string time_units = "min"; 

	/*========================*/
	/* Microenvironment setup */
	/*========================*/
		
	setup_microenvironment(world, cart_topo);  
	
	//User parameters
	double time_add_tnf = parameters.ints("time_add_tnf");
	double time_put_tnf = 0;
	double duration_add_tnf = parameters.ints("duration_add_tnf");
	double time_tnf_next = 0;
	double time_remove_tnf = parameters.ints("time_remove_tnf");
	
	
	double concentration_tnf = parameters.doubles("concentration_tnf") * 0.1;
	
	//Radius around which the tnf pulse is injected
	double membrane_lenght = parameters.ints("membrane_length");
	
	//TNF density index
	static int tnf_idx = microenvironment.find_density_index("tnf");	
	
	bool seed_tnf = false;
	
	//Do small diffusion steps to initialize densities in case seed_tnf is true
	if ( seed_tnf )
	{
		//inject_density_sphere(tnf_idx, concentration_tnf, membrane_lenght);
		inject_density_sphere(tnf_idx, concentration_tnf, membrane_lenght, world, cart_topo);
		for ( int i = 0; i < 25; i ++ )
			microenvironment.simulate_diffusion_decay( diffusion_dt, world, cart_topo );
	}

	/*=============================*/	
	/* PhysiCell/PhysiCell-X setup */ 
	/*=============================*/	
 	
	// Mechanical voxel size must be >= Diffusion voxel size (check with PhysiCell_settings.xml)
	double mechanics_voxel_size = 20; 
	

	// Calling the parallel version of Cell Container creation 
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size, world, cart_topo );
	
	//----> Users typically start modifying here. START USERMODS 
	
	create_cell_types();
	
	//Calling the parallel version of setup_tissue(...) 
	setup_tissue(microenvironment, world, cart_topo);

	//----> Users typically stop modifying here. END USERMODS  
	
	//Set MultiCellDS save options <--- These MUST all be 'true' here 
	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true ); 
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true );
	
	//Save a simulation snapshot 
	
	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() ); 
	
	//Use the parallel version of the function for XML file writing
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time, world, cart_topo ); 
	
	//Save SVG cross section through z = 0, after setting its length bar to 200 microns 
	PhysiCell_SVG_options.length_bar = 200; 

	//For simplicity, set a pathology coloring function 
	std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function;
	
	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() );
	
	//Use the parallel version of the function for SVG plot file
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, world, cart_topo );
		
	//Set the performance timers 
	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();

	std::ofstream report_file;
	if( world.rank == 0 )
	{
		sprintf(filename , "%s/simulation_report.txt" , PhysiCell_settings.folder.c_str() ); 
		report_file.open(filename); 	// create the data log file 
		report_file << "timepoint\talive\tapoptotic\tnectortic";
		report_file << "\ttotal_free_tnfr\ttotal_active_tnfr\ttotal_int_TNF\ttotal_tnf"<<std::endl;
	}

	if(IOProcessor(world))
		std::cout << PhysiCell_settings.enable_legacy_saves << std::endl;
	//Main loop of the program 
	try 
	{		
		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt )
		{
			//Save data if it's time. 
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
			{
				if(IOProcessor(world))
					std::cout << "Time to save" << std::endl;
					
				//Use the parallel version of the function	
				display_simulation_status( std::cout, world, cart_topo );
				
				double timepoint = PhysiCell_globals.current_time;
				
				//Count Necrotic, Apoptotic and Alive cells
				int alive_no, necrotic_no, apoptotic_no;
				float total_tnf, total_free_tnfr, total_active_tnfr, total_int_TNF;
				float total_active_TNF, total_active_FADD, total_active_NFKb;
				
				//Call the parallel versions of the function now which use MPI_Reduce at rank 0 
				alive_no 		  = total_live_cell_count(world, cart_topo);
				necrotic_no 	  = total_necrosis_cell_count(world, cart_topo);
				apoptotic_no 	  = total_dead_cell_count(world, cart_topo);
				total_tnf         = get_total_tnf(world, cart_topo);
				total_free_tnfr   = total_free_TNF_receptor(world, cart_topo);
				total_active_tnfr = total_active_TNF_receptor(world, cart_topo);
				total_int_TNF     = total_internalized_TNF_receptor(world, cart_topo);
				total_active_TNF  = total_active_TNF_node(world, cart_topo);
				total_active_FADD = total_active_FADD_node(world, cart_topo);
				total_active_NFKb = total_active_NFKb_node(world, cart_topo);

				if( world.rank == 0) 
				{
					report_file << PhysiCell_globals.current_time << "\t" << alive_no << "\t" << necrotic_no << "\t" << apoptotic_no;
					report_file << "\t" << total_free_tnfr << "\t" << total_active_tnfr << "\t" << total_int_TNF;
					report_file << "\t" << total_active_TNF << "\t" << total_active_FADD << "\t" << total_active_NFKb;
					report_file << "\t" << total_tnf  <<std::endl;
				}
				
				if( PhysiCell_settings.enable_full_saves == true )
				{	
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
					
					// Use parallel version of function 
					save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time, world, cart_topo ); 
					// Use serial version
					// save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time);
				}
				
				PhysiCell_globals.full_output_index++; 
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
			}
			
			//Save SVG plot if it's time
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 0.01 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{	
					
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index ); 
					SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, world, cart_topo );
					
					PhysiCell_globals.SVG_output_index++; 
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
				}
			}

			//Custom add-ons could potentially go here. 
		
			if ( PhysiCell_globals.current_time >= time_put_tnf )
			{
				time_tnf_next = PhysiCell_globals.current_time + duration_add_tnf;
				time_put_tnf += time_add_tnf;
			}

			if ( PhysiCell_globals.current_time >= time_remove_tnf )
			{
				remove_density(tnf_idx);													
				time_remove_tnf += PhysiCell_settings.max_time;
			}

			if ( PhysiCell_globals.current_time <= time_tnf_next )
			{
				inject_density_sphere(tnf_idx, concentration_tnf, membrane_lenght, world, cart_topo);
			}

			//Update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt, world, cart_topo );
			
			//Update the TNF receptor model of each cell
			tnf_receptor_model_main( diffusion_dt );		
			
			//Run PhysiCell/PhysiCell-X 
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time, world, cart_topo );
			
			PhysiCell_globals.current_time += diffusion_dt;
		}
		if (world.rank == 0)
		{
			report_file.close();
		}
	}
	catch( const std::exception& e )
	{ 
		//Reference to the base of a polymorphic object, information from length_error printed
		std::cout << e.what();  
	}
	
	//Save a final simulation snapshot 
	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time, world, cart_topo );
	// save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time); 
	
	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, world, cart_topo );
	// SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function);
	
	
	//Timer
	if(IOProcessor(world)) 
	{
		std::cout << std::endl << "Total simulation runtime: " << std::endl;
		BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 
 	 }
  
  	//Gracefully shut-down MPI i.e. distributed parallelization
	world.Finalize(); 
	
	return 0; 
}


