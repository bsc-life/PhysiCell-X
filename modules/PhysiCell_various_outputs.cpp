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
# Copyright (c) 2015-2025, Paul Macklin and the PhysiCell Project             #
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

#include "../core/PhysiCell.h"
#include "./PhysiCell_various_outputs.h"

#include "./PhysiCell_settings.h"

#include <iostream>
#include <fstream>
#include <unistd.h>

namespace PhysiCell{

int writePov(std::vector<Cell*> all_cells, double timepoint, double scale)
{
	static int TUMOR_TYPE=0; 
	static int VESSEL_TYPE=1; 
	
	std::string filename; 
	filename.resize( 1024 ); 
//	sprintf( (char*) filename.c_str() , "output//cells_%i.pov" , (int)round(timepoint) ); 
	sprintf( (char*) filename.c_str() , "%s/cells_%i.pov" , PhysiCell_settings.folder.c_str() ,  (int)round(timepoint) ); 
	std::ofstream povFile (filename.c_str(), std::ofstream::out);
	povFile<<"#include \"colors.inc\" \n";
	povFile<<"#include \"header.inc\" \n";
	
	for(int i=0;i<all_cells.size();i++)
	{
		std::string _nameCore;

		if( all_cells[i]->phenotype.cycle.pCycle_Model )
		{
			int code= all_cells[i]->phenotype.cycle.current_phase().code;
			if (code ==PhysiCell_constants::Ki67_positive_premitotic || code==PhysiCell_constants::Ki67_positive_postmitotic || code==PhysiCell_constants::Ki67_positive || code==PhysiCell_constants::Ki67_negative || code==PhysiCell_constants::live)
				_nameCore="LIVE";
			else if (code==PhysiCell_constants::apoptotic)
				_nameCore="APOP";
			else if (code==PhysiCell_constants::necrotic_swelling || code==PhysiCell_constants::necrotic_lysed || code==PhysiCell_constants::necrotic)
				_nameCore="NEC";
			else if (code==PhysiCell_constants::debris)
				_nameCore="DEBR";
			else
				_nameCore="MISC";
		}
		else if(all_cells[i]->type==TUMOR_TYPE)
			_nameCore="LIVE";
		else if(all_cells[i]->type==VESSEL_TYPE)
			_nameCore="ENDO";
		else
			_nameCore="MISC";
		std::string center= "<" + std::to_string(all_cells[i]->position[0]/scale) + "," + std::to_string(all_cells[i]->position[1]/scale) +","+ std::to_string(all_cells[i]->position[2]/scale) +">";
		std::string core = "sphere {\n\t" + center + "\n\t " + std::to_string( all_cells[i]->phenotype.geometry.radius/scale) + "\n\t FinishMacro ( " + center +","+ _nameCore+ "Finish,"+ _nameCore + "*1)\n}\n";
		povFile<< core;		
	}
	
	povFile<<"#include \"footer.inc\" \n";
	povFile.close();
	return 0;
}

int writePov(std::vector<Cell*> all_cells, double timepoint, double scale,mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	static int TUMOR_TYPE=0; 
	static int VESSEL_TYPE=1; 
	
	std::string filename; 
	filename.resize( 1024 ); 
//	sprintf( (char*) filename.c_str() , "output//cells_%i.pov" , (int)round(timepoint) ); 
	sprintf( (char*) filename.c_str() , "%s/cells_%i.pov" , PhysiCell_settings.folder.c_str() ,  (int)round(timepoint) ); 
	std::ofstream povFile (filename.c_str(), std::ofstream::out);
	if (world.rank == 0) {
		povFile<<"#include \"colors.inc\" \n";
		povFile<<"#include \"header.inc\" \n";
	}
	
	for (int ser_ctr = 0; ser_ctr < world.size; ser_ctr++) {
		if (world.rank == ser_ctr) {
			for(int i=0;i<all_cells.size();i++)
			{
				std::string _nameCore;

				if( all_cells[i]->phenotype.cycle.pCycle_Model )
				{
					int code= all_cells[i]->phenotype.cycle.current_phase().code;
					if (code ==PhysiCell_constants::Ki67_positive_premitotic || code==PhysiCell_constants::Ki67_positive_postmitotic || code==PhysiCell_constants::Ki67_positive || code==PhysiCell_constants::Ki67_negative || code==PhysiCell_constants::live)
						_nameCore="LIVE";
					else if (code==PhysiCell_constants::apoptotic)
						_nameCore="APOP";
					else if (code==PhysiCell_constants::necrotic_swelling || code==PhysiCell_constants::necrotic_lysed || code==PhysiCell_constants::necrotic)
						_nameCore="NEC";
					else if (code==PhysiCell_constants::debris)
						_nameCore="DEBR";
					else
						_nameCore="MISC";
				}
				else if(all_cells[i]->type==TUMOR_TYPE)
					_nameCore="LIVE";
				else if(all_cells[i]->type==VESSEL_TYPE)
					_nameCore="ENDO";
				else
					_nameCore="MISC";
				std::string center= "<" + std::to_string(all_cells[i]->position[0]/scale) + "," + std::to_string(all_cells[i]->position[1]/scale) +","+ std::to_string(all_cells[i]->position[2]/scale) +">";
				std::string core = "sphere {\n\t" + center + "\n\t " + std::to_string( all_cells[i]->phenotype.geometry.radius/scale) + "\n\t FinishMacro ( " + center +","+ _nameCore+ "Finish,"+ _nameCore + "*1)\n}\n";
				povFile<< core;		
			}
		}
		MPI_Barrier(cart_topo.mpi_cart_comm);
	}
	
	if (world.rank == world.size - 1) {
		povFile<<"#include \"footer.inc\" \n";
		povFile.close();
	}	
	return 0;
}

int writeCellReport(std::vector<Cell*> all_cells, double timepoint)
{
	std::string filename; 
	filename.resize( 1024 ); 
//	sprintf( (char*) filename.c_str() , "output//cells_%i.txt" , (int)round(timepoint) ); 
	sprintf( (char*) filename.c_str() , "%s/cells_%i.txt" , PhysiCell_settings.folder.c_str() , (int)round(timepoint) ); 
	std::ofstream povFile (filename.c_str(), std::ofstream::out);
	povFile<<"\tID\tx\ty\tz\tradius\tvolume_total\tvolume_nuclear_fluid\tvolume_nuclear_solid\tvolume_cytoplasmic_fluid\tvolume_cytoplasmic_solid\tvolume_calcified_fraction\tphenotype\telapsed_time\n";
	int phenotype_code;
	for(int i=0;i<all_cells.size();i++)
	{
		phenotype_code = all_cells[i]->phenotype.cycle.current_phase().code;
		// phenotype_code = phases.size()>0?all_cells[i]->phenotype.cycle.phases[all_cells[i]->phenotype.current_phase_index].code:-1;
		povFile<<i<<"\t"<<all_cells[i]->ID<<"\t"<<all_cells[i]->position[0]<<"\t" << all_cells[i]->position[1] <<"\t"<< all_cells[i]->position[2]<<"\t";
		povFile<<all_cells[i]->phenotype.geometry.radius<<"\t"<<all_cells[i]->phenotype.volume.total<<"\t"<<all_cells[i]->phenotype.volume.nuclear_fluid
		<<"\t"<<all_cells[i]->phenotype.volume.nuclear_solid<<"\t"<<all_cells[i]->phenotype.volume.cytoplasmic_fluid<<"\t"<<
		all_cells[i]->phenotype.volume.cytoplasmic_solid<<"\t"<<all_cells[i]->phenotype.volume.calcified_fraction<<"\t"<<phenotype_code<< 
		// "\t"<< all_cells[i]->phenotype.cycle.phases[all_cells[i]->phenotype.current_phase_index].elapsed_time <<std::endl;		
		"\t"<< all_cells[i]->phenotype.cycle.data.elapsed_time_in_phase <<std::endl;		
	}
	povFile.close();
	return 0;
}

int writeCellReport(std::vector<Cell*> all_cells, double timepoint, mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	std::string filename; 
	filename.resize( 1024 ); 
//	sprintf( (char*) filename.c_str() , "output//cells_%i.txt" , (int)round(timepoint) ); 
	sprintf( (char*) filename.c_str() , "%s/cells_%i.txt" , PhysiCell_settings.folder.c_str() , (int)round(timepoint) ); 
	std::ofstream povFile (filename.c_str(), std::ofstream::out);
	if (IOProcessor(world)) {
		povFile<<"\tID\tx\ty\tz\tradius\tvolume_total\tvolume_nuclear_fluid\tvolume_nuclear_solid\tvolume_cytoplasmic_fluid\tvolume_cytoplasmic_solid\tvolume_calcified_fraction\tphenotype\telapsed_time\n";
	}
	int phenotype_code;
	
	for (int ser_ctr = 0; ser_ctr < world.size; ser_ctr++) {
		if (world.rank == ser_ctr) {
			for(int i=0;i<all_cells.size();i++)
			{
				phenotype_code = all_cells[i]->phenotype.cycle.current_phase().code;
				// phenotype_code = phases.size()>0?all_cells[i]->phenotype.cycle.phases[all_cells[i]->phenotype.current_phase_index].code:-1;
				povFile<<i<<"\t"<<all_cells[i]->ID<<"\t"<<all_cells[i]->position[0]<<"\t" << all_cells[i]->position[1] <<"\t"<< all_cells[i]->position[2]<<"\t";
				povFile<<all_cells[i]->phenotype.geometry.radius<<"\t"<<all_cells[i]->phenotype.volume.total<<"\t"<<all_cells[i]->phenotype.volume.nuclear_fluid
				<<"\t"<<all_cells[i]->phenotype.volume.nuclear_solid<<"\t"<<all_cells[i]->phenotype.volume.cytoplasmic_fluid<<"\t"<<
				all_cells[i]->phenotype.volume.cytoplasmic_solid<<"\t"<<all_cells[i]->phenotype.volume.calcified_fraction<<"\t"<<phenotype_code<< 
				// "\t"<< all_cells[i]->phenotype.cycle.phases[all_cells[i]->phenotype.current_phase_index].elapsed_time <<std::endl;		
				"\t"<< all_cells[i]->phenotype.cycle.data.elapsed_time_in_phase <<std::endl;		
			}
		}
		MPI_Barrier(cart_topo.mpi_cart_comm);
	}
	povFile.close();
	return 0;
}

void display_simulation_status( std::ostream& os )
{
	os << "current simulated time: " << PhysiCell_globals.current_time << " " << 
		PhysiCell_settings.time_units << " (max: " << 
		PhysiCell_settings.max_time << " " << 
		PhysiCell_settings.time_units << ")" << std::endl; 
		
	os << "total agents: " << all_cells->size() << std::endl; 
	
	os << "interval wall time: ";
	BioFVM::TOC();
	BioFVM::display_stopwatch_value( os , BioFVM::stopwatch_value() ); 
	os << std::endl; 
	BioFVM::TIC(); 

	os << "MaBoSS updates: " << count_maboss_updates << " | avg wall time/update: ";
	if( count_maboss_updates > 0 )
	{
		os << ( time_maboss_update / count_maboss_updates ) / 1000.0 << " ms";
	}
	else
	{
		os << "n/a";
	}
	os << std::endl;

#ifdef MECHS_TIME
	if( count_mechanics_timing_steps > 0 )
	{
		const double avg_scale = 1.0 / static_cast<double>( count_mechanics_timing_steps ) / 1000.0;
		os << "MECHS_TIME avg/step: "
			<< "gradients " << time_mechs_gradient_total * avg_scale << " ms | "
			<< "halos " << time_mechs_halo_total * avg_scale << " ms | "
			<< "potentials " << time_mechs_potential_total * avg_scale << " ms" << std::endl;
		os << "                    "
			<< "interactions " << time_mechs_interactions_total * avg_scale << " ms | "
			<< "clear dummy " << time_mechs_clear_dummy_total * avg_scale << " ms | "
			<< "positions " << time_mechs_update_position_total * avg_scale << " ms" << std::endl;
		os << "                    "
			<< "pack " << time_mechs_pack_total * avg_scale << " ms | "
			<< "transfer " << time_mechs_transfer_total * avg_scale << " ms | "
			<< "unpack " << time_mechs_unpack_total * avg_scale << " ms | "
			<< "voxel update " << time_mechs_voxels_update_total * avg_scale << " ms" << std::endl;
	}
#endif
		
	os << "total wall time: "; 
	BioFVM::RUNTIME_TOC();
	BioFVM::display_stopwatch_value( os , BioFVM::runtime_stopwatch_value() ); 
	os << std::endl << std::endl; 

	time_maboss_update = 0.0;
	count_maboss_updates = 0;
#ifdef MECHS_TIME
	time_mechs_gradient_total = 0.0;
	time_mechs_halo_total = 0.0;
	time_mechs_potential_total = 0.0;
	time_mechs_interactions_total = 0.0;
	time_mechs_clear_dummy_total = 0.0;
	time_mechs_update_position_total = 0.0;
	time_mechs_pack_total = 0.0;
	time_mechs_transfer_total = 0.0;
	time_mechs_unpack_total = 0.0;
	time_mechs_voxels_update_total = 0.0;
	count_mechanics_timing_steps = 0;
#endif
		
	return;
}

/*----------------------------------------*/
/* Parallel version of the function above */
/*----------------------------------------*/

double getProcessMemoryGB() {
    std::ifstream statm("/proc/self/statm");
    long rssPages = 0;
    statm >> rssPages >> rssPages; // second value is RSS (resident set size)
    long pageSize = sysconf(_SC_PAGESIZE); // in bytes
    double memUsedGB = (rssPages * pageSize) / (1024.0 * 1024.0 * 1024.0);
    return memUsedGB;
}

void display_simulation_status( std::ostream& os, mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	if(IOProcessor(world))
	{
		os << "current simulated time: " << PhysiCell_globals.current_time << " " << 
			PhysiCell_settings.time_units << " (max: " << 
			PhysiCell_settings.max_time << " " << 
			PhysiCell_settings.time_units << ")" << std::endl; 
	}
	
	int local_cells  = all_cells->size(); 
	int global_cells; 
	MPI_Reduce(&local_cells, &global_cells, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm); 
	
	double global_total, global_time_diff, global_time_mech, global_time_pheno, global_time_intracell;
	double global_time_mech_sum = 0.0;
	double global_time_maboss_update = 0.0;
	double local_ram = getProcessMemoryGB();
	double global_ram = 0.0;
	unsigned long long global_maboss_updates = 0;
	unsigned long long global_mechanics_cells = 0;
	unsigned long long global_mechanics_neighbor_candidates = 0;
	unsigned long long global_mechanics_neighbor_interactions = 0;
	
	MPI_Reduce(&time_diff, &global_time_diff, 1, MPI_DOUBLE, MPI_MAX, 0, cart_topo.mpi_cart_comm);
	MPI_Reduce(&time_diff, &global_time_diff, 1, MPI_DOUBLE, MPI_MAX, 0, cart_topo.mpi_cart_comm);
	MPI_Reduce(&time_mechs, &global_time_mech, 1, MPI_DOUBLE, MPI_MAX, 0, cart_topo.mpi_cart_comm);
	MPI_Reduce(&time_mechs, &global_time_mech_sum, 1, MPI_DOUBLE, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	MPI_Reduce(&time_pheno, &global_time_pheno, 1, MPI_DOUBLE, MPI_MAX, 0, cart_topo.mpi_cart_comm);
	MPI_Reduce(&time_intracell, &global_time_intracell, 1, MPI_DOUBLE, MPI_MAX, 0, cart_topo.mpi_cart_comm);
	MPI_Reduce(&time_maboss_update, &global_time_maboss_update, 1, MPI_DOUBLE, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	MPI_Reduce(&count_maboss_updates, &global_maboss_updates, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	MPI_Reduce(&count_mechanics_cells, &global_mechanics_cells, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	MPI_Reduce(&count_mechanics_neighbor_candidates, &global_mechanics_neighbor_candidates, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	MPI_Reduce(&count_mechanics_neighbor_interactions, &global_mechanics_neighbor_interactions, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	MPI_Reduce(&local_ram, &global_ram, 1, MPI_DOUBLE, MPI_SUM, 0, cart_topo.mpi_cart_comm);

	if(IOProcessor(world))
	{	

		os << "total agents: " << global_cells << std::endl; 
	
		os << "interval wall time: ";
		BioFVM::TOC();
		BioFVM::display_stopwatch_value( os , BioFVM::stopwatch_value() ); 
		os <<  " | "  << global_ram << "GB "  << std::endl; 
		BioFVM::TIC(); 

			double global_total= global_time_diff + global_time_mech + global_time_pheno + global_time_intracell;
			os << "time scales: ";
			os << " diff " << (global_time_diff/ global_total) * 100 << "% | ";
			os << " mechs " << (global_time_mech/ global_total) * 100 << "% | ";
			os << " pheno " << (global_time_pheno/ global_total) * 100 << "% | ";
			os << " intracell " << (global_time_intracell/ global_total) * 100 << "%  " << std::endl;
			
			os << "MaBoSS updates: " << global_maboss_updates << " | avg wall time/update: ";
			if( global_maboss_updates > 0 )
			{
				os << ( global_time_maboss_update / global_maboss_updates ) / 1000.0 << " ms";
			}
			else
			{
				os << "n/a";
			}
			os << std::endl;

			os << "Mechanics cells: " << global_mechanics_cells << " | avg wall time/cell: ";
			if( global_mechanics_cells > 0 )
			{
				os << ( global_time_mech_sum / global_mechanics_cells ) / 1000.0 << " ms";
			}
			else
			{
				os << "n/a";
			}
			os << " | avg candidate neighbors/cell: ";
			if( global_mechanics_cells > 0 )
			{
				os << static_cast<double>( global_mechanics_neighbor_candidates ) / static_cast<double>( global_mechanics_cells );
			}
			else
			{
				os << "n/a";
			}
			os << " | avg interacting neighbors/cell: ";
			if( global_mechanics_cells > 0 )
			{
				os << static_cast<double>( global_mechanics_neighbor_interactions ) / static_cast<double>( global_mechanics_cells );
			}
			else
			{
				os << "n/a";
			}
			os << std::endl;

#ifdef MECHS_TIME
			if( count_mechanics_timing_steps > 0 )
			{
				const double avg_scale = 1.0 / static_cast<double>( count_mechanics_timing_steps ) / 1000.0;
				os << "MECHS_TIME avg/step: "
					<< "gradients " << time_mechs_gradient_total * avg_scale << " ms | "
					<< "halos " << time_mechs_halo_total * avg_scale << " ms | "
					<< "potentials " << time_mechs_potential_total * avg_scale << " ms" << std::endl;
				os << "                    "
					<< "interactions " << time_mechs_interactions_total * avg_scale << " ms | "
					<< "clear dummy " << time_mechs_clear_dummy_total * avg_scale << " ms | "
					<< "positions " << time_mechs_update_position_total * avg_scale << " ms" << std::endl;
				os << "                    "
					<< "pack " << time_mechs_pack_total * avg_scale << " ms | "
					<< "transfer " << time_mechs_transfer_total * avg_scale << " ms | "
					<< "unpack " << time_mechs_unpack_total * avg_scale << " ms | "
					<< "voxel update " << time_mechs_voxels_update_total * avg_scale << " ms" << std::endl;
			}
#endif
			
			os << "total wall time: "; 
		BioFVM::RUNTIME_TOC();
		BioFVM::display_stopwatch_value( os , BioFVM::runtime_stopwatch_value() ); 
		os << std::endl << std::endl;
	} 

	time_diff = 0.0;
	time_mechs = 0.0;
	time_pheno = 0.0;
	time_intracell = 0.0;
	time_maboss_update = 0.0;
	count_maboss_updates = 0;
	count_mechanics_cells = 0;
	count_mechanics_neighbor_candidates = 0;
	count_mechanics_neighbor_interactions = 0;
#ifdef MECHS_TIME
	time_mechs_gradient_total = 0.0;
	time_mechs_halo_total = 0.0;
	time_mechs_potential_total = 0.0;
	time_mechs_interactions_total = 0.0;
	time_mechs_clear_dummy_total = 0.0;
	time_mechs_update_position_total = 0.0;
	time_mechs_pack_total = 0.0;
	time_mechs_transfer_total = 0.0;
	time_mechs_unpack_total = 0.0;
	time_mechs_voxels_update_total = 0.0;
	count_mechanics_timing_steps = 0;
#endif
		
	return;
}

void log_output(double t, int output_index, Microenvironment microenvironment, std::ofstream& report_file)
{
	double scale=1000;
	int num_new_cells= 0;
	int num_deaths=0;
//	std::cout << "current simulated time: " << t   << " minutes " << std::endl; 
//	std::cout << "interval wall time: ";
//	BioFVM::TOC();
//	BioFVM::display_stopwatch_value( std::cout , BioFVM::stopwatch_value() ); 
//	std::cout << std::endl; 
//	std::cout << "total wall time: "; 
//	BioFVM::RUNTIME_TOC();
//	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 
//	std::cout << std::endl;
	
	std::cout << "time: "<<t<<std::endl;
	num_new_cells=t==0?all_basic_agents.size():((Cell_Container *)microenvironment.agent_container)->num_divisions_in_current_step;
	num_deaths=((Cell_Container *)microenvironment.agent_container)->num_deaths_in_current_step;
	std::cout<<"total number of agents (newly born, deaths): " << (*all_cells).size()<<"("<<num_new_cells<<", "<<num_deaths<<")" << std::endl; 
	report_file<<t<<"\t"<<(*all_cells).size()<<"\t"<<num_new_cells<<"\t"<<num_deaths<<"\t"<<BioFVM::stopwatch_value()<< std::endl; 
//	BioFVM::TIC();
	
	((Cell_Container *)microenvironment.agent_container)->num_divisions_in_current_step=0;
	((Cell_Container *)microenvironment.agent_container)->num_deaths_in_current_step=0;
	writePov(*all_cells, t, scale);
	writeCellReport(*all_cells, t);
	std::string filename; 
	filename.resize( 1024 , '\0' ); 
	sprintf( (char*) filename.c_str() , "output%08d.mat" , output_index ); 
	filename.resize( strlen( filename.c_str() ) ); 
	// std::cout << "\tWriting to file " << filename << " ... " << std::endl; 
	// microenvironment.write_to_matlab( filename ); 
	
	return;
}

void log_output(double t, int output_index, Microenvironment microenvironment, std::ofstream& report_file, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	double scale=1000;

	int global_num_new_cells= 0;
	int global_num_deaths=0;
	int global_cells = 0;
	int local_cells  = all_cells->size();

	MPI_Reduce(&((Cell_Container *)microenvironment.agent_container)->num_divisions_in_current_step, &global_num_new_cells, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	MPI_Reduce(&((Cell_Container *)microenvironment.agent_container)->num_deaths_in_current_step, &global_num_deaths, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	MPI_Reduce(&local_cells, &global_cells, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm);
	std::cout << "time: "<<t<<std::endl;
	
	std::cout<<"total number of agents (newly born, deaths): " << global_cells <<"("<<global_num_new_cells<<", "<<global_num_deaths<<")" << std::endl; 
	report_file<<t<<"\t"<<(*all_cells).size()<<"\t"<<global_num_new_cells<<"\t"<<global_num_deaths<<"\t"<<BioFVM::stopwatch_value()<< std::endl; 

	
	((Cell_Container *)microenvironment.agent_container)->num_divisions_in_current_step=0;
	((Cell_Container *)microenvironment.agent_container)->num_deaths_in_current_step=0;
	writePov(*all_cells, t, scale, world, cart_topo);
	writeCellReport(*all_cells, t, world, cart_topo);
	std::string filename; 
	filename.resize( 1024 , '\0' ); 
	sprintf( (char*) filename.c_str() , "output%08d.mat" , output_index ); 
	filename.resize( strlen( filename.c_str() ) ); 

	
	return;
}

	
};
