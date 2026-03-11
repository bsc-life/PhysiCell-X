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

#ifndef __PhysiCell_cell_container_h__
#define __PhysiCell_cell_container_h__

#include <vector>
#include <unordered_map>
#include "PhysiCell_cell.h"
#include "../BioFVM/BioFVM_agent_container.h"
#include "../BioFVM/BioFVM_mesh.h"
#include "../BioFVM/BioFVM_microenvironment.h"
#include "../DistPhy/DistPhy_Environment.h"
#include "../DistPhy/DistPhy_Cartesian.h"
#include "../DistPhy/DistPhy_Utils.h"

using namespace DistPhy::mpi; 

namespace PhysiCell{
extern double time_diff, time_mechs, time_pheno, time_intracell;
class Cell;

/*==========================================================================================*/
/* Gaurav Saxena added these 2 classes before class Cell_Container as objects of these 			*/
/* 2 classes are used inside Cell_Container class i.e. Moore_Cell_Info and Moore_Voxel_Info */
/* No need for forward declaration now as we have fully defined classes now 								*/
/*==========================================================================================*/

class Moore_Cell_Info
{
public:
	int ID;
	std::vector<double> position;
	double radius;
	double nuclear_radius;
	double cell_cell_repulsion_strength;
	double relative_maximum_adhesion_distance;
	double cell_cell_adhesion_strength;
	int type;

	Moore_Cell_Info() = default;
	
	Moore_Cell_Info(std::vector<char>& buffer, int size, int &pos);


};

class Moore_Voxel_Info 
{
public:
	int global_mesh_index;
	int no_of_cells_in_vxl;
	std::vector<double> center;
	double max_cell_interactive_distance_in_voxel;
	std::vector<Moore_Cell_Info> moore_cells;

	Moore_Voxel_Info() = default;

	Moore_Voxel_Info(std::vector<char>& buffer, int size, int &pos);
};

class Interacting_Cell_Info
{
public:
	int ID;
	int type;
	bool dead;
	std::vector<double> position;
	//Attack info
	bool attacked = false;
	double damage_suffered = 0;
	double time_attacked = 0;
	//Fuse info
	//IN
	int number_of_nuclei;
	double phenotype_volume;
	double cytoplasmic_fluid;
	double nuclear_fluid;
	double cytoplasmic_solid;
	double nuclear_solid;
	std::vector<double> internalized_substrates;
	double target_solid_cytoplasmic;
	double target_solid_nuclear;
	//OUT
	bool fused = false;

	//Ingest
	//IN
	std::vector<double> fraction_transferred_when_ingested;
	//OUT
	bool ingested = false;

};

class Interacting_Voxel
{
public:
	int global_mesh_index;  //used to id in u_map
	int no_of_cells_in_vxl; //necessary to unpack info
	//std::vector<double> center;
	//double max_cell_interactive_distance_in_voxel;
	std::vector<Interacting_Cell_Info> cells; 
};


class Cell_Container : public BioFVM::Agent_Container
{
 private:	
	std::vector<Cell*> cells_ready_to_divide; 			// the index of agents ready to divide
	std::vector<Cell*> cells_ready_to_die;
	int boundary_condition_for_pushed_out_agents; 	// what to do with pushed out cells
	bool initialzed = false;
	
 public:
	BioFVM::Cartesian_Mesh underlying_mesh;
	std::vector<double> max_cell_interactive_distance_in_voxel;
	int num_divisions_in_current_step = 0;	//In-class initializer, c++11
	int num_deaths_in_current_step = 0; //In-class initializer, c++11
	
	/*-----------------------------------------------------------*/
	/* Added by Gaurav Saxena: left/right refer to the 'process' */
	/* snd_buf_left - send to left process											 */
	/* snd_buf_right - send to right process										 */
	/* rcv_buf_left - rcv from left process											 */
	/* rcv_buf_right - rcv from right process									   */
	/*-----------------------------------------------------------*/
	
	std::vector<char> snd_buf_left;
	std::vector<char> snd_buf_right;
	std::vector<char> rcv_buf_left;
	std::vector<char> rcv_buf_right;
	
	/*--------------------------------------------------------------------------------------------------*/
	/* Added by Gaurav Saxena																																						*/
	/* position_left gives where the position in the snd_buf_left where the element being packed goes		*/
	/* position_right gives where the position in the snd_buf_right where the element being packed goes	*/
	/* Ultimately position_left/right give the length of the snd_buf's in bytes directly								*/
	/* They are reset to 0 in the pack() function.																											*/
	/*--------------------------------------------------------------------------------------------------*/
	
	int position_left;
	int position_right;
	
	/*--------------------------------------------------------------------------------------------------*/
	/* Added by Gaurav Saxena																																						*/
	/* no_cell_cross_left/right give the number of cells crossing to left or right. 										*/
	/* They are reset to a value 0 in the pack() function.																							*/
	/*--------------------------------------------------------------------------------------------------*/
	
	int no_cells_cross_left ;
	int no_cells_cross_right ; 
	
	/*--------------------------------------------------------------------------------------------------*/
	/* Added by Gaurav Saxena																																						*/
	/* no_of_cells_from_left/right give the number of cells coming from left/right neighbour 						*/
	/*--------------------------------------------------------------------------------------------------*/
	
	int no_of_cells_from_right;
	int no_of_cells_from_left; 

	double last_diffusion_time  = 0.0; 
	double last_cell_cycle_time = 0.0;
	double last_mechanics_time  = 0.0;
	
	Cell_Container();
 	void initialize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double voxel_size);
    
    /*-------------------------------------*/
    /* Parallel prototype of function above*/
    /*-------------------------------------*/
    void initialize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double voxel_size, mpi_Environment &world, mpi_Cartesian &cart_topo);
    
	void initialize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx, double dy, double dz);
    
    /*-------------------------------------*/
    /* Parallel prototype of function above*/
    /*-------------------------------------*/
    void initialize(double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx, double dy, double dz, mpi_Environment &world, mpi_Cartesian &cart_topo);
    
	std::vector<std::vector<Cell*> > agent_grid;
	std::vector<std::vector<Cell*> > agents_in_outer_voxels;
	
	void update_all_cells(double t);
	
		/*-------------------------------------*/
    /* Parallel prototype of function above*/
    /*-------------------------------------*/
    
    void update_all_cells(double t, mpi_Environment &world, mpi_Cartesian &cart_topo); 
	
	void update_all_cells(double t, double dt);
	void update_all_cells(double t, double phenotype_dt, double mechanics_dt);
	void update_all_cells(double t, double phenotype_dt, double mechanics_dt, double diffusion_dt ); 
	
		/*-------------------------------------*/
    /* Parallel prototype of function above*/
    /*-------------------------------------*/
  
    void update_all_cells(double t, double phenotype_dt, double mechanics_dt, double diffusion_dt, mpi_Environment &world, mpi_Cartesian &cart_topo );
	void update_cell_potentials(double time_since_last_mechanics,  mpi_Environment &world, mpi_Cartesian &cart_topo);

	void register_agent( Cell* agent );
	void add_agent_to_outer_voxel(Cell* agent);
	void remove_agent(Cell* agent );
	void remove_agent_from_voxel(Cell* agent, int voxel_index);
	void add_agent_to_voxel(Cell* agent, int voxel_index);
	
	void flag_cell_for_division( Cell* pCell ); 
	void flag_cell_for_removal( Cell* pCell ); 
	bool contain_any_cell(int voxel_index);

	Cell* find_cell(int local_voxel, int cell_id);
	
	/*-----------------------------------------------------------------*/
	/* Added by Gaurav Saxena, a new function which would 'byte' pack	 */
	/* the Cell data structure using MPI_Pack() and send it to the left*/
	/* or right process.																							 */
	/*-----------------------------------------------------------------*/
	
	void pack(std::vector<Cell*> *all_cells, mpi_Environment &world, mpi_Cartesian &cart_topo);
	void unpack(mpi_Environment &world, mpi_Cartesian &cart_topo); //As of now this is a dummy function
	
	/*------------------------------------------------------------------------*/
	/* Added by Gaurav Saxena 																								*/
	/* velocity of a cell is updated using (1) other cells in its voxel	 			*/
	/* (2) cells in 26 neighbouring voxels (a physically neighbouring 	 			*/
	/* voxel MAY NOT be a logical neighbour) 														 			*/
	/* Functions below help to pack, exchange, unpack some required fields		*/
	/* of a cell and voxels in neighbouring processes which are needed to 		*/
	/* update the velocity. 																									*/
	/* Later change name of function to exchange_moore_info 									*/
	/*------------------------------------------------------------------------*/
	void pack_moore_voxel(uint voxel_index, std::vector<char>& snd_buffer, int& len_buffer, int& position);
	void pack_moore_info(mpi_Environment &world, mpi_Cartesian &cart_topo);
	
	/*-----------------------------------------------------------------------------------------------------*/
	/* Added by Gaurav Saxena - mbfr is a vector that contains info of boundary voxels in adjacent process */
	/* and specific data fields of all the cells in these voxels, similarly mbfl 													 */
	/*-----------------------------------------------------------------------------------------------------*/
	
	std::vector<Moore_Voxel_Info> mbfr;		//Moore Boundary From Right (mbfr from right process)
	std::vector<Moore_Voxel_Info> mbfl; 	//Moore Boundary From Left 	(mbfl from left process)
	
	std::unordered_map<int, Moore_Voxel_Info> um_mbfl; //um_mbfl[globa_mesh_index]=Moore_Voxel_Info[index]
	std::unordered_map<int, Moore_Voxel_Info> um_mbfr; //um_mbfr[globa_mesh_index]=Moore_Voxel_Info[index]

	/*-----------------------------------------------------------------------------------------------------*/
	/* Internal cell information communication for cell cell interactions								   */
	/* It is required read information and return result to neighbours processes 						   */													 
	/*-----------------------------------------------------------------------------------------------------*/

	void pack_moore_info_parallel(mpi_Environment &world, mpi_Cartesian &cart_topo);

	void pack_cell_interact_info(mpi_Environment &world, mpi_Cartesian &cart_topo);
	void unpack_cell_interact_info(mpi_Environment &world, mpi_Cartesian &cart_topo);
	void cell_cell_interaction_with_border(std::vector<Interacting_Voxel> *iv);
	void exchange_mechanics_halos(mpi_Environment &world, mpi_Cartesian &cart_topo);
	void evaluate_cell_cell_interactions(double time_since_last_mechanics,  mpi_Environment &world, mpi_Cartesian &cart_topo);
	void evaluate_cell_elastic_interactions( PhysiCell::Cell *pCell , double dt_mec,  mpi_Environment &world, mpi_Cartesian &cart_topo);

	
	std::vector<Interacting_Voxel> ivfr;		//Voxels From Right (mbfr from right process)
	std::vector<Interacting_Voxel> ivfl; 	//Voxels Boundary From Left 	(mbfl from left process)
	
	std::unordered_map<int, Interacting_Voxel> um_ivfl; //um_mbfl[globa_mesh_index]=Moore_Voxel_Info[index]
	std::unordered_map<int, Interacting_Voxel> um_ivfr; //um_mbfr[globa_mesh_index]=Moore_Voxel_Info[index]

	std::vector<Interacting_Cell_Info> get_neighbour_interacting_cells(Cell* pCell, mpi_Environment &world, mpi_Cartesian &cart_topo);
	std::vector<Moore_Cell_Info> get_moore_neighbour_cells(Cell* pCell, mpi_Environment &world, mpi_Cartesian &cart_topo);
};

int find_escaping_face_index(Cell* agent);
extern std::vector<Cell*> *all_cells; 

Cell_Container* create_cell_container_for_microenvironment( BioFVM::Microenvironment& m , double mechanics_voxel_size );
Cell_Container* create_cell_container_for_microenvironment( BioFVM::Microenvironment& m , double mechanics_voxel_size, mpi_Environment &world, mpi_Cartesian &cart_topo ); 

/*-----------------------------------------------------------------------------*/
/* Added by Gaurav Saxena: To update the velocity a cell requires all the cells*/
/* in the neighbouring 26 voxels. For voxels at the sub-domain boundary, cells */
/* in the adjoining process boundary must all be packed and sent across.			 */
/* They need to be unpacked in a specially made class called Moore_Voxel_Info	 */
/* which contains Voxel specific info AND info of ALL the cells in that Voxel	 */
/* The latter info is a class in itself called Moore_Cell_Info 								 */
/*-----------------------------------------------------------------------------*/


};
#endif
