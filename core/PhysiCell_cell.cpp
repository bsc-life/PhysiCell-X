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

#include "PhysiCell_cell.h"
#include "PhysiCell_cell_container.h"
#include "PhysiCell_utilities.h"
#include "PhysiCell_constants.h"
#include "../BioFVM/BioFVM_vector.h" 
#include<limits.h>

namespace PhysiCell
{

Cell_Parameters::Cell_Parameters()
{
	o2_hypoxic_threshold = 15.0; // HIF-1alpha at half-max around 1.5-2%, and tumors often are below 2%
	o2_hypoxic_response = 8.0; // genomic / proteomic changes observed at 7-8 mmHg 
	o2_hypoxic_saturation = 4.0; // maximum HIF-1alpha at 0.5% o2 (McKeown)
	
	o2_necrosis_threshold = 5.0; 
	o2_necrosis_max = 2.5; 
	
	o2_proliferation_threshold = 5.0; // assume no proliferation at same level as starting necrosis 
	o2_proliferation_saturation = 160.0; // 5% = 38, 21% = 160 mmHg 
	o2_reference = 160.0; // assume all was measured in normoxic 21% o2 
	
	pReference_live_phenotype = NULL; // reference live (usually physioxic) phenotype 
	
	// necrosis parameters 
	
	max_necrosis_rate = 1.0 / (6.0 * 60.0); // assume cells survive 6 hours in very low oxygen 
	necrosis_type = PhysiCell_constants::deterministic_necrosis;;

	return; 
}

Cell_Definition::Cell_Definition()
{
	// set the microenvironment pointer 
	pMicroenvironment = NULL;
	if( BioFVM::get_default_microenvironment() != NULL )
	{ pMicroenvironment = BioFVM::get_default_microenvironment(); }

	// set up the default parameters 
		// the default Cell_Parameters constructor should take care of this
		
	type = 0; 
	name = "unnamed"; 

	parameters.pReference_live_phenotype = &phenotype; 
		
	// set up the default custom data 
		// the default Custom_Cell_Data constructor should take care of this
		
	// set up the default functions 
	functions.volume_update_function = NULL; // standard_volume_update_function;
	functions.update_migration_bias = NULL; 
	
	functions.update_phenotype = NULL; 
	functions.custom_cell_rule = NULL; 
	
	functions.update_velocity = NULL; // standard_update_cell_velocity;
	functions.add_cell_basement_membrane_interactions = NULL; 
	functions.calculate_distance_to_membrane = NULL; 
	
	functions.set_orientation = NULL;
	
	return; 
}

Cell_Definition::Cell_Definition( Cell_Definition& cd )
{
	// set the microenvironment pointer 
	pMicroenvironment = cd.pMicroenvironment;

	// set up the default parameters 
		// the default Cell_Parameters constructor should take care of this
		
	type = cd.type; 
	name = cd.name; 
	 
	parameters = cd.parameters;
	custom_data = cd.custom_data; 
	functions = cd.functions; 
	phenotype = cd.phenotype; 
	
	// this is the whole reason we need ot make a copy constructor 
	parameters.pReference_live_phenotype = &phenotype; 
	
	return; 
}

Cell_Definition& Cell_Definition::operator=( const Cell_Definition& cd )
{
	// set the microenvironment pointer 
	pMicroenvironment = cd.pMicroenvironment;

	// set up the default parameters 
		// the default Cell_Parameters constructor should take care of this
		
	type = cd.type; 
	name = cd.name; 
	 
	parameters = cd.parameters;
	custom_data = cd.custom_data; 
	functions = cd.functions; 
	phenotype = cd.phenotype; 
	
	// this is the whole reason we need ot make a copy constructor 
	parameters.pReference_live_phenotype = &phenotype; 
	
	return *this; 
}


Cell_Definition cell_defaults; 

Cell_State::Cell_State()
{
	neighbors.resize(0); 
	orientation.resize( 3 , 0.0 ); 
	
	simple_pressure = 0.0; 
	
	return; 
}

void Cell::update_motility_vector( double dt_ )
{
	if( phenotype.motility.is_motile == false )
	{
		phenotype.motility.motility_vector.assign( 3, 0.0 ); 
		return; 
	}
	
	if( UniformRandom() < dt_ / phenotype.motility.persistence_time || phenotype.motility.persistence_time < dt_ )
	{
		// choose a uniformly random unit vector 
		double temp_angle = 6.28318530717959*UniformRandom();
		double temp_phi = 3.1415926535897932384626433832795*UniformRandom();
		
		double sin_phi = sin(temp_phi);
		double cos_phi = cos(temp_phi);
		
		if( phenotype.motility.restrict_to_2D == true )
		{ 
			sin_phi = 1.0; 
			cos_phi = 0.0;
		}
		
		std::vector<double> randvec; 
		randvec.resize(3,sin_phi); 
		
		randvec[0] *= cos( temp_angle ); // cos(theta)*sin(phi)
		randvec[1] *= sin( temp_angle ); // sin(theta)*sin(phi)
		randvec[2] = cos_phi; //  cos(phi)
		
		// if the update_bias_vector function is set, use it  
		if( functions.update_migration_bias )
		{
			functions.update_migration_bias( this,phenotype,dt_ ); 
		}
		
		phenotype.motility.motility_vector = phenotype.motility.migration_bias_direction; // motiltiy = bias_vector
		phenotype.motility.motility_vector *= phenotype.motility.migration_bias; // motility = bias*bias_vector 
		
		double one_minus_bias = 1.0 - phenotype.motility.migration_bias; 
		
		axpy( &(phenotype.motility.motility_vector), one_minus_bias, randvec ); // motility = (1-bias)*randvec + bias*bias_vector
		
		normalize( &(phenotype.motility.motility_vector) ); 
		
		phenotype.motility.motility_vector *= phenotype.motility.migration_speed; 
	}	
	return; 
} 

void Cell::advance_bundled_phenotype_functions( double dt_ )
{
	// call the custom code to update the phenotype 
	if( functions.update_phenotype )
	{	functions.update_phenotype( this , phenotype , dt_ ); }
	
	// update volume 
	if( functions.volume_update_function )
	{
		functions.volume_update_function(this,phenotype,dt_); 
		
		// The following line is needed in every volume 
		// regulation method (it sets BioFVM total_volume)
		
		set_total_volume( phenotype.volume.total ); 
	}
	
	// update geometry
	phenotype.geometry.update( this, phenotype, dt_ );
	
	// check for new death events 
	if( phenotype.death.check_for_death( dt_ ) == true )
	{
		// if so, change the cycle model to the current death model 
		phenotype.cycle.sync_to_cycle_model( phenotype.death.current_model() ); 
		
		// also, turn off motility.
		
		phenotype.motility.is_motile = false; 
		phenotype.motility.motility_vector.assign( 3, 0.0 ); 
		functions.update_migration_bias = NULL;
		
		// turn off secretion, and reduce uptake by a factor of 10 
		phenotype.secretion.set_all_secretion_to_zero();
		phenotype.secretion.scale_all_uptake_by_factor( 0.10 );
		
		// make sure to run the death entry function 
		if( phenotype.cycle.current_phase().entry_function )
		{
			phenotype.cycle.current_phase().entry_function( this, phenotype, dt_ ); 
		}
	}
	
	// advance cycle model (for both cell cycle and death cycle models)
	phenotype.cycle.advance_cycle( this, phenotype, dt_ ); 
	if( phenotype.flagged_for_removal )
	{
		flag_for_removal(); 
		phenotype.flagged_for_removal = false; 
	}
	if( phenotype.flagged_for_division )
	{
		flag_for_division(); 
		phenotype.flagged_for_division = false; 
	}
	
	return; 
}

Cell::Cell()
{
	// use the cell defaults; 
	
	type = cell_defaults.type; 
	type_name = cell_defaults.name; 
	
	custom_data = cell_defaults.custom_data; 
	parameters = cell_defaults.parameters; 
	functions = cell_defaults.functions; 
	
	phenotype = cell_defaults.phenotype; 
	
	phenotype.molecular.sync_to_cell( this ); 
	
	// cell state should be fine by the default constructor 
	
	current_mechanics_voxel_index=-1;
	
	updated_current_mechanics_voxel_index = 0;
	
	is_movable = true;
	is_out_of_domain = false;
	
	/*-----------------------------------*/
	/* Added by Gaurav Saxena to say that*/
	/* when cell created, both = 0			 */
	/*-----------------------------------*/
	
	crossed_to_left_subdomain = false;
	crossed_to_right_subdomain = false; 

	
	displacement.resize(3,0.0); // state? 
	
	assign_orientation();
	container = NULL;
	
	
	return; 
}

/*--------------------------------------------------------------------------------------------------*/
/* Parallel version of Cell class Constructor which takes as input the cell id in parallel settings */
/* It must call explicitly a new constructor of Basic_Agent which takes as input this cell id				*/
/*--------------------------------------------------------------------------------------------------*/

Cell::Cell(int p_ID):Basic_Agent(p_ID) //----> Correct syntax for calling parametrized base class constructor
{
	// use the cell defaults; 
	
	//Basic_Agent(p_ID);---> this syntax is incorrect to call Base class parametrized constructor
	type = cell_defaults.type; 
	type_name = cell_defaults.name; 
	
	custom_data = cell_defaults.custom_data; 
	parameters 	= cell_defaults.parameters; 
	functions 	= cell_defaults.functions; 
	
	phenotype 	= cell_defaults.phenotype; 
	
	phenotype.molecular.sync_to_cell( this ); 
	
	// cell state should be fine by the default constructor 
	
	current_mechanics_voxel_index					=-1;	
	updated_current_mechanics_voxel_index = 0;
	
	is_movable 				= true;
	is_out_of_domain 	= false;
	displacement.resize(3,0.0); // state? 
	
	assign_orientation();				//Just assigns a random unit vector to the cell. 
	container = NULL;
	
	
	return; 
}


void Cell::flag_for_division( void )
{
	get_container()->flag_cell_for_division( this );
	return; 
}

void Cell::flag_for_removal( void )
{
	get_container()->flag_cell_for_removal( this );
	return;
}

void Cell::start_death( int death_model_index )
{
	// set the death data struture to the indicated death model 
	phenotype.death.trigger_death( death_model_index ); 
	// change the cycle model to the current death model 
	phenotype.cycle.sync_to_cycle_model( phenotype.death.current_model() ); 
		
	// turn off secretion, and reduce uptake by a factor of 10 
	phenotype.secretion.set_all_secretion_to_zero();
	phenotype.secretion.scale_all_uptake_by_factor( 0.10 );
		
	// turn off motility.
	phenotype.motility.is_motile = false; 
	phenotype.motility.motility_vector.assign( 3, 0.0 ); 
	functions.update_migration_bias = NULL;
		
	// make sure to run the death entry function 
	if( phenotype.cycle.current_phase().entry_function )
	{
		phenotype.cycle.current_phase().entry_function( this, phenotype, 0.0 ); 
	}

	return; 
}

void Cell::assign_orientation()
{
	state.orientation.resize(3,0.0);
	if( functions.set_orientation != NULL )
	{
		functions.set_orientation(this, phenotype, 0.0 );
	}
	else
	{
		//assign a random unit vector
		double theta= UniformRandom()*6.28318530717959; //rand*2*pi
		double z= 2* UniformRandom()-1;
		double temp= sqrt(1-z*z);
		state.orientation[0]= temp * cos(theta);
		state.orientation[1]= temp * sin(theta);
		state.orientation[2]= z;
	}
	
	return; 
}

Cell* Cell::divide( )
{
	// phenotype.flagged_for_division = false; 
	// phenotype.flagged_for_removal = false; 
	
	Cell* child = create_cell();
	child->copy_data( this );	
	child->copy_function_pointers(this);
	child->parameters = parameters;
	
	// evenly divide internalized substrates 
	// if these are not actively tracked, they are zero anyway 
	*internalized_substrates *= 0.5; 
	*(child->internalized_substrates) = *internalized_substrates ; 
	
	// The following is already performed by create_cell(). JULY 2017 ***
	// child->register_microenvironment( get_microenvironment() );
	
	// randomly place the new agent close to me, accounting for orientation and 
	// polarity (if assigned)
		
	double temp_angle = 6.28318530717959*UniformRandom();
	double temp_phi = 3.1415926535897932384626433832795*UniformRandom();
	
	double radius= phenotype.geometry.radius;
	std::vector<double> rand_vec (3, 0.0);
	
	rand_vec[0]= cos( temp_angle ) * sin( temp_phi );
	rand_vec[1]= sin( temp_angle ) * sin( temp_phi );
	rand_vec[2]= cos( temp_phi );
	rand_vec = rand_vec- phenotype.geometry.polarity*(rand_vec[0]*state.orientation[0]+ 
		rand_vec[1]*state.orientation[1]+rand_vec[2]*state.orientation[2])*state.orientation;
	
	if( norm(rand_vec) < 1e-16 )
	{
		std::cout<<"************ERROR********************"<<std::endl;
	}
	normalize( &rand_vec ); 
	// rand_vec/= norm(rand_vec);
	child->assign_position(position[0] + 0.5 * radius*rand_vec[0],
						 position[1] + 0.5 * radius*rand_vec[1],
						 position[2] + 0.5 * radius*rand_vec[2]);
	//change my position to keep the center of mass intact and then see if I need to update my voxel index
	static double negative_one_half = -0.5; 
	naxpy( &position, negative_one_half , rand_vec );// position = position - 0.5*rand_vec; 
	// position[0] -= 0.5*radius*rand_vec[0];
	// position[1] -= 0.5*radius*rand_vec[1]; 
	// position[2] -= 0.5*radius*rand_vec[2]; 
	
	//If this cell has been moved outside of the boundaries, mark it as such.
	//(If the child cell is outside of the boundaries, that has been taken care of in the assign_position function.)
	if( !get_container()->underlying_mesh.is_position_valid(position[0], position[1], position[2])){
		is_out_of_domain = true;
		is_active = false;
		is_movable = false;
	}	
	 
	update_voxel_in_container();
	phenotype.volume.divide(); 
	child->phenotype.volume.divide();
	child->set_total_volume(child->phenotype.volume.total);
	set_total_volume(phenotype.volume.total);
	
	// child->set_phenotype( phenotype ); 
	child->phenotype = phenotype; 
	
	return child;
}

/*-------------------------------------------------------------------------*/
/* Parallel version of the divide() function above												 */
/*-------------------------------------------------------------------------*/

Cell* Cell::divide(int p_ID, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	// phenotype.flagged_for_division = false; 
	// phenotype.flagged_for_removal = false; 
	
	Cell* child = create_cell(p_ID);					//Calls new version of create_cell function.
	
	child->copy_data( this );	
	child->copy_function_pointers(this);
	child->parameters = parameters;
	
	// evenly divide internalized substrates 
	// if these are not actively tracked, they are zero anyway 
	*internalized_substrates *= 0.5; 
	*(child->internalized_substrates) = *internalized_substrates ; 
	
	// The following is already performed by create_cell(). JULY 2017 ***
	// child->register_microenvironment( get_microenvironment() );
	
	// randomly place the new agent close to me, accounting for orientation and 
	// polarity (if assigned)
		
	double temp_angle = 6.28318530717959*UniformRandom();
	double temp_phi = 3.1415926535897932384626433832795*UniformRandom();
	
	double radius= phenotype.geometry.radius;
	std::vector<double> rand_vec (3, 0.0);
	
	rand_vec[0]= cos( temp_angle ) * sin( temp_phi );
	rand_vec[1]= sin( temp_angle ) * sin( temp_phi );
	rand_vec[2]= cos( temp_phi );
	rand_vec = rand_vec- phenotype.geometry.polarity*(rand_vec[0]*state.orientation[0]+ 
		rand_vec[1]*state.orientation[1]+rand_vec[2]*state.orientation[2])*state.orientation;
	
	if( norm(rand_vec) < 1e-16 )
	{
		std::cout<<"************ERROR********************"<<std::endl;
	}
	normalize( &rand_vec ); 
	// rand_vec/= norm(rand_vec);
	child->assign_position(position[0] + 0.5 * radius*rand_vec[0],
						 						 position[1] + 0.5 * radius*rand_vec[1],
						 						 position[2] + 0.5 * radius*rand_vec[2], 
						 						 world, cart_topo);
	//change my position to keep the center of mass intact and then see if I need to update my voxel index
	
	
	/*------------------------------------------------------------------------------------*/
	/* SIDENOTE: In naxpy() below we are passing &position. position is declared as 			*/
	/* std::vector<double> position. Now, "position" does not mean the address of the 		*/
	/* vector i.e. "position" is NOT a pointer. When we pass &position, it doesn't mean		*/
	/* that we are passing a pointer to a pointer. It JUST means we are passing the 			*/
	/* address of the vector - which is caught in std::vector<double> *y	i.e. a pointer	*/
	/* a type "vector". This is a big difference between arrays and vectors. 							*/
	/*------------------------------------------------------------------------------------*/
	
	/*------------------------------------------------------------------------------------*/
	/* TICKET: depending on the solution of the ticket, I may have to change the value of	*/
	/* negative_one_half = +0.5 and then call naxpy(). 																		*/
	/*------------------------------------------------------------------------------------*/
	
	static double negative_one_half = -0.5; 
	naxpy( &position, negative_one_half , rand_vec );// position = position - 0.5*rand_vec;
	 
	// position[0] -= 0.5*radius*rand_vec[0];
	// position[1] -= 0.5*radius*rand_vec[1]; 
	// position[2] -= 0.5*radius*rand_vec[2]; 
	
	//If this cell has been moved outside of the boundaries, mark it as such.
	//(If the child cell is outside of the boundaries, that has been taken care of in the assign_position function.)
	
	if( !get_container()->underlying_mesh.is_position_valid(position[0], position[1], position[2])){
		is_out_of_domain = true;
		is_active = false;
		is_movable = false;
	}	
	
	/*-----------------------------------------------------------------------------------*/
	/* Now the parent cell has also moved after division, so there is a possibility that */
	/* it may not be respecting sub-domain boundaries - therefore we need to call 			 */
	/* correct_position_within_subdomain() for the parent cell.													 */
	/*-----------------------------------------------------------------------------------*/
	
	get_container()->underlying_mesh.correct_position_within_subdomain(position, world, cart_topo); 
	 	
	update_voxel_in_container(world, cart_topo);			//Though world, cart_topo are not needed in this, just to distinguish
	
	phenotype.volume.divide(); 												//Don't think anything to parallelize
	child->phenotype.volume.divide();
	
	child->set_total_volume(child->phenotype.volume.total);
	set_total_volume(phenotype.volume.total);
	
	// child->set_phenotype( phenotype ); 
	child->phenotype = phenotype; 
	
	return child;
}


bool Cell::assign_position(std::vector<double> new_position)
{
	return assign_position(new_position[0], new_position[1], new_position[2]);
}

void Cell::set_previous_velocity(double xV, double yV, double zV)
{
	previous_velocity[0] = xV;
	previous_velocity[1] = yV;
	previous_velocity[2] = zV;

	return; 
}

bool Cell::assign_position(double x, double y, double z)
{
	position[0]=x;
	position[1]=y;
	position[2]=z;
	
	// update microenvironment current voxel index
	update_voxel_index();
	// update current_mechanics_voxel_index
	current_mechanics_voxel_index= get_container()->underlying_mesh.nearest_voxel_index( position );
	get_container()->register_agent(this);
	
	if( !get_container()->underlying_mesh.is_position_valid(x,y,z) )
	{	
		is_out_of_domain = true; 
		is_active = false; 
		is_movable = false; 
		
		return false;
	}
	
	return true;
}

/*-------------------------------------------------------*/
/* Parallel version of assign_position() function				 */
/*-------------------------------------------------------*/

bool Cell::assign_position(double x, double y, double z, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	position[0]=x;
	position[1]=y;
	position[2]=z;
	
	// update microenvironment current voxel index - this becomes the local voxel index in parallel settings
	
	/*------------------------------------------------------------------------------------*/
	/* update_voxel_index(world, cart_topo) is in BioFVM_basic_agent.cpp									*/
	/* First it calls is_position_valid() for the "Diffusion" mesh. 											*/
	/* If position invalid then sets out_of_domain=true, current_voxel_index=-1 & returns	*/
	/* Otherwise calls nearest_voxel_local_index()<---new function (GS).									*/
	/*------------------------------------------------------------------------------------*/

	update_voxel_index(world, cart_topo);
	
	//update current_mechanics_voxel_index
	//current_mechanics_voxel_index= get_container()->underlying_mesh.nearest_voxel_index( position);
	//Changed the call above to a parallel call.
	
	current_mechanics_voxel_index= get_container()->underlying_mesh.nearest_voxel_local_index( position, world, cart_topo);

	get_container()->register_agent(this);
	
	/*--------------------------------------------------------------------------------------------*/
	/* What is the need to call this from underlying_mesh ? the boundaries of the physical domain */
	/* remain the same.  																																					*/
	/* IMPORTANT: update_voxel_index() above calls a "position" modifying function created by 		*/
	/* Gaurav Saxena which adjusts cell positions if it doesn't respect sub-domain boundaries.		*/
	/* Hence, when calling !get_container()->underlying_mesh.is_position_valid(x,y,z), it is 			*/
	/* BETTER to replace (x,y,z) with (position[0], position[1], position[2]) as after calling		*/
	/* the function correct_position_within_subdomain(), MAYBE position[0] != x, position[1] != y,*/ 
	/* position[2] !=z. BUT for NOW let it remain like this, THINK then replace										*/
	/*--------------------------------------------------------------------------------------------*/

	if( !get_container()->underlying_mesh.is_position_valid(x,y,z) )
	{	
		is_out_of_domain = true; 
		is_active = false; 
		is_movable = false; 
		
		return false;
	}
	
	return true;
}

void Cell::set_total_volume(double volume)
{
	Basic_Agent::set_total_volume(volume);
	
	// If the new volume is significantly different than the 
	// current total volume, adjust all the sub-volumes 
	// proportionally. 
	
	// if( fabs( phenotype.volume.total - volume ) < 1e-16 )
	if( fabs( phenotype.volume.total - volume ) > 1e-16 )
	{
		double ratio= volume/ phenotype.volume.total;
		phenotype.volume.multiply_by_ratio(ratio);
	}
	
	phenotype.geometry.update( this, phenotype, 0.0 ); 
	// phenotype.update_radius();
	//if( get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] < 
	//	phenotype.geometry.radius * parameters.max_interaction_distance_factor )
	if( get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] < 
		phenotype.geometry.radius * phenotype.mechanics.relative_maximum_adhesion_distance )
	{
		// get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()]= phenotype.geometry.radius*parameters.max_interaction_distance_factor;
		get_container()->max_cell_interactive_distance_in_voxel[get_current_mechanics_voxel_index()] = phenotype.geometry.radius
			* phenotype.mechanics.relative_maximum_adhesion_distance;
	}
	
	return; 
}

double& Cell::get_total_volume(void)
{
	static bool I_warned_you = false; 
	if( I_warned_you == false )
	{
		std::cout << "Warning! Do not use " << __FUNCTION__ << "!" << std::endl 
			<< "Use (some_cell).phenotype.volume.total instead!" << std::endl; 
		I_warned_you = true; 
	}
	return phenotype.volume.total; 
}

void Cell::turn_off_reactions(double dt)
{	
	is_active = false;  
	
	for(int i=0;i< phenotype.secretion.uptake_rates.size();i++)
	{
		phenotype.secretion.uptake_rates[i] = 0.0;  
		phenotype.secretion.secretion_rates[i] = 0.0; 
	}
	set_internal_uptake_constants(dt);
	
	return; 
}

Cell_Container * Cell::get_container()
{
	if(container == NULL)
	{
		container = (Cell_Container *)get_microenvironment()->agent_container;
	}
	
	return container;
}

void Cell::die()
{
	delete_cell(this);
	return; 
}

/*--------------------------------------------------------------------------------------*/
/* A function needed in the distributed parallel scenario which removes a cell that has */ 
/* crossed to the neighbouring sub-domain. 																							*/
/*--------------------------------------------------------------------------------------*/

void Cell::remove_crossed_cell()
{
	 vaporize_teleported_cell(this->index); 
}

void Cell::update_position( double dt )
{
	// BioFVM Basic_Agent::update_position(dt) returns without doing anything. 
	// So we remove this to avoid any future surprises. 
	// 
	// Basic_Agent::update_position(dt);
		
	// use Adams-Bashforth 
	static double d1; 
	static double d2; 
	static bool constants_defined = false; 
	if( constants_defined == false )
	{
		d1 = dt; 
		d1 *= 1.5; 
		d2 = dt; 
		d2 *= -0.5; 
		constants_defined = true; 
	}
	
	// new AUgust 2017
	if( default_microenvironment_options.simulate_2D == true )
	{ velocity[2] = 0.0; }
	
	std::vector<double> old_position(position); 
	axpy( &position , d1 , velocity );  
	axpy( &position , d2 , previous_velocity );  
	// overwrite previous_velocity for future use 
	// if(sqrt(dist(old_position, position))>3* phenotype.geometry.radius)
		// std::cout<<sqrt(dist(old_position, position))<<"old_position: "<<old_position<<", new position: "<< position<<", velocity: "<<velocity<<", previous_velocity: "<< previous_velocity<<std::endl;
	
	previous_velocity = velocity; 
	
	velocity[0]=0; velocity[1]=0; velocity[2]=0;
	// #pragma omp critical
	//{update_voxel_in_container();}
	if(get_container()->underlying_mesh.is_position_valid(position[0],position[1],position[2]))
	{
		updated_current_mechanics_voxel_index=get_container()->underlying_mesh.nearest_voxel_index( position );
	}
	else
	{
		updated_current_mechanics_voxel_index=-1;
		
		is_out_of_domain = true; 
		is_active = false; 
		is_movable = false; 
	}
	return; 
}

/*----------------------------------------------*/
/* Parallel Implementation of update_position() */
/*----------------------------------------------*/

void Cell::update_position( double dt, mpi_Environment &world )
{
	/* use Adams-Bashforth */
	
	static double d1; 
	static double d2; 
	static bool constants_defined = false; 
	
	if( constants_defined == false )
	{
		d1  = dt; 
		d1 *= 1.5; 
		d2  = dt; 
		d2 *= -0.5; 
		constants_defined = true; 
	}
	

	if(default_microenvironment_options.simulate_2D == true ) 
		 velocity[2] = 0.0; 
	
	
	std::vector<double> old_position(position); 
	axpy( &position , d1 , velocity );  
	axpy( &position , d2 , previous_velocity );
		 
	/* overwrite previous_velocity for future use */
		
	previous_velocity = velocity; 
	
	velocity[0]=0; 
	velocity[1]=0; 
	velocity[2]=0;
	
	crossed_to_left_subdomain  = get_container()->underlying_mesh.has_it_crossed_to_left_subdomain(position[0], world);
	crossed_to_right_subdomain = get_container()->underlying_mesh.has_it_crossed_to_right_subdomain(position[0], world); 
		
	bool crossed = crossed_to_left_subdomain || crossed_to_right_subdomain;
	if(crossed == true)
	{
			std::cout<<std::endl; 
			std::cout<<"Rank="<<world.rank<<" Cell ID="<<ID<<" Old Position:("<<old_position[0]<<","<<old_position[1]<<","<<old_position[2]<<")"<<std::endl;  
			std::cout<<"Cell ID="<<ID<<" New Position:("<<position[0]<<","<<position[1]<<","<<position[2]<<")"<<std::endl;  
	}
	
	/*------------------------------------------------------------------------------*/
	/* If it has crossed to left OR right subdomain (from a non-boundary processes) */
	/* then we don't need to check if it is out-of-domain as it CANNOT be. 					*/
	/* If BOTH _left and _right above are false, then there is a possibility 				*/
	/* that process is boundary process and cell has gone out of physical boundary	*/
	/* In this case, the code below will execute 																		*/
	/* In case the cell is still within sub-domain, a new voxel index IS recomputed */
	/*------------------------------------------------------------------------------*/
	 
	if(crossed == false)
	{
		if(get_container()->underlying_mesh.is_position_valid(position[0],position[1],position[2]))
		{
			updated_current_mechanics_voxel_index=get_container()->underlying_mesh.nearest_voxel_index( position );
		}
		else
		{
			updated_current_mechanics_voxel_index=-1;
			is_out_of_domain = true; 
			is_active = false; 
			is_movable = false; 
		}
	}
	return; 
}

int Cell::get_current_mechanics_voxel_index()
{
	return current_mechanics_voxel_index;
}

void Cell::update_voxel_in_container()
{
	// call the method from BioFVM_basic_agent to update microenvironment's voxel index
	update_voxel_index();
	// int temp_current_voxel_index;
	// Check to see if we need to remove agents that are pushed out of boundary
	// if(!get_container()->underlying_mesh.is_position_valid(position[0],position[1],position[2]))	
		
	if(updated_current_mechanics_voxel_index==-1)// updated_current_mechanics_voxel_index is updated in update_position
	{
		// check if this agent has a valid voxel index, if so, remove it from previous voxel
		if( get_current_mechanics_voxel_index() >= 0)
		{
			// #pragma omp critical
			{get_container()->remove_agent_from_voxel(this, get_current_mechanics_voxel_index());}
		}
		// #pragma omp critical
		{get_container()->add_agent_to_outer_voxel(this);}
		// std::cout<<"cell out of boundary..."<< __LINE__<<" "<<ID<<std::endl;
		current_mechanics_voxel_index=-1;
		is_out_of_domain=true;
		is_active=false;
		return;
	}
	
	// temp_current_voxel_index= get_current_mechanics_voxel_index();
	// updated_current_mechanics_voxel_index=get_container()->underlying_mesh.nearest_voxel_index( position );
	
	// update mesh indices (if needed)
	if(updated_current_mechanics_voxel_index!= get_current_mechanics_voxel_index())
	{
		// #pragma omp critical
		{
			container->remove_agent_from_voxel(this, get_current_mechanics_voxel_index());
			container->add_agent_to_voxel(this, updated_current_mechanics_voxel_index);
		}
		current_mechanics_voxel_index=updated_current_mechanics_voxel_index;
	}
	
	return; 
}

/*------------------------------------------------------------------*/
/* Parallel version of function: update_voxel_in_container() above  */
/*------------------------------------------------------------------*/

void Cell::update_voxel_in_container(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	
	update_voxel_index(world, cart_topo);							 //This is in BioFVM_basic_agent.cpp
	
	// Check to see if we need to remove agents that are pushed out of boundary
			
	if(updated_current_mechanics_voxel_index==-1)			//This is updated in update_position() and -1 indicates outside physical boundary 
	{
		if( get_current_mechanics_voxel_index() >= 0)		//If agent had valid voxel index, need to remove agent from that voxel
			  get_container()->remove_agent_from_voxel(this, get_current_mechanics_voxel_index());			
		
		
		get_container()->add_agent_to_outer_voxel(this); //std::cout<<"cell out of boundary..."<< __LINE__<<" "<<ID<<std::endl;		
		current_mechanics_voxel_index=-1;
		is_out_of_domain=true;
		is_active=false;
		return;
	}
	
	// update mesh indices (if needed)
	if(updated_current_mechanics_voxel_index!= get_current_mechanics_voxel_index())
	{
			container->remove_agent_from_voxel(this, get_current_mechanics_voxel_index());
			container->add_agent_to_voxel(this, updated_current_mechanics_voxel_index);
			current_mechanics_voxel_index=updated_current_mechanics_voxel_index;
	}	
	return; 
}

void Cell::copy_data(Cell* copy_me)
{
	// phenotype=copyMe-> phenotype; //it is taken care in set_phenotype
	type = copy_me->type; 
	type_name = copy_me->type_name; 
	
	custom_data = copy_me->custom_data; 
	parameters = copy_me->parameters; 
	
	velocity = copy_me->velocity; 
	// expected_phenotype = copy_me-> expected_phenotype; //it is taken care in set_phenotype
	cell_source_sink_solver_temp1 = std::vector<double>(copy_me->cell_source_sink_solver_temp1);
	cell_source_sink_solver_temp2 = std::vector<double>(copy_me->cell_source_sink_solver_temp2);
	
	return; 
}

void Cell::copy_function_pointers(Cell* copy_me)
{
	functions = copy_me->functions; 
	return; 
}

void Cell::add_potentials(Cell* other_agent)
{
	if( this->ID == other_agent->ID )
	{ return; }

	// 12 uniform neighbors at a close packing distance, after dividing out all constants
	static double simple_pressure_scale = 0.027288820670331; // 12 * (1 - sqrt(pi/(2*sqrt(3))))^2 
	// 9.820170012151277; // 12 * ( 1 - sqrt(2*pi/sqrt(3)))^2

	double distance = 0; 
	for( int i = 0 ; i < 3 ; i++ ) 
	{ 
		displacement[i] = position[i] - (*other_agent).position[i]; 
		distance += displacement[i] * displacement[i]; 
	}
	// Make sure that the distance is not zero
	
	distance = std::max(sqrt(distance), 0.00001); 
	
	//Repulsive
	double R = phenotype.geometry.radius+ (*other_agent).phenotype.geometry.radius; 
	
	double RN = phenotype.geometry.nuclear_radius + (*other_agent).phenotype.geometry.nuclear_radius;	
	double temp_r, c;
	if( distance > R ) 
	{
		temp_r=0;
	}
	// else if( distance < RN ) 
	// {
		// double M = 1.0; 
		// c = 1.0 - RN/R; 
		// c *= c; 
		// c -= M; 
		// temp_r = ( c*distance/RN  + M  ); 
	// }
	else
	{
		// temp_r = 1 - distance/R;
		temp_r = -distance; // -d
		temp_r /= R; // -d/R
		temp_r += 1.0; // 1-d/R
		temp_r *= temp_r; // (1-d/R)^2 
		
		// add the relative pressure contribution 
		state.simple_pressure += ( temp_r / simple_pressure_scale ); // New July 2017 
	}
	
	// August 2017 - back to the original if both have same coefficient 
	double effective_repulsion = sqrt( phenotype.mechanics.cell_cell_repulsion_strength * other_agent->phenotype.mechanics.cell_cell_repulsion_strength ); 
	temp_r *= effective_repulsion; 
	
	// temp_r *= phenotype.mechanics.cell_cell_repulsion_strength; // original 
	//////////////////////////////////////////////////////////////////
	
	// Adhesive
	//double max_interactive_distance = parameters.max_interaction_distance_factor * phenotype.geometry.radius + 
	//	(*other_agent).parameters.max_interaction_distance_factor * (*other_agent).phenotype.geometry.radius;
		
	double max_interactive_distance = phenotype.mechanics.relative_maximum_adhesion_distance * phenotype.geometry.radius + 
		(*other_agent).phenotype.mechanics.relative_maximum_adhesion_distance * (*other_agent).phenotype.geometry.radius;
		
	if(distance < max_interactive_distance ) 
	{	
		// double temp_a = 1 - distance/max_interactive_distance; 
		double temp_a = -distance; // -d
		temp_a /= max_interactive_distance; // -d/S
		temp_a += 1.0; // 1 - d/S 
		temp_a *= temp_a; // (1-d/S)^2 
		// temp_a *= phenotype.mechanics.cell_cell_adhesion_strength; // original 
		
		// August 2017 - back to the original if both have same coefficient 
		double effective_adhesion = sqrt( phenotype.mechanics.cell_cell_adhesion_strength * other_agent->phenotype.mechanics.cell_cell_adhesion_strength ); 
		temp_a *= effective_adhesion; 
		
		temp_r -= temp_a;
	}
	/////////////////////////////////////////////////////////////////
	if( fabs(temp_r) < 1e-16 )
	{ return; }
	temp_r /= distance;
	// for( int i = 0 ; i < 3 ; i++ ) 
	// {
	//	velocity[i] += displacement[i] * temp_r; 
	// }
	axpy( &velocity , temp_r , displacement ); 
	
	return;
}

Cell* create_cell( void )
{
	Cell* pNew; 
	pNew = new Cell;		
	(*all_cells).push_back( pNew ); 
	pNew->index=(*all_cells).size()-1;
	
	// new usability enhancements in May 2017 
	
	if( BioFVM::get_default_microenvironment() )
	{
		pNew->register_microenvironment( BioFVM::get_default_microenvironment() );
	}

	// All the phenotype and other data structures are already set 
	// by virtue of the default Cell constructor. 
	
	return pNew; 
}

// in the future, I might swap this with create_cell(): 
// In that "create_cell()" uses "create_cell( cell_defaults )" 
Cell* create_cell( Cell_Definition& cd )
{
	Cell* pNew = create_cell(); 
	
	// use the cell defaults; 
	pNew->type = cd.type; 
	pNew->type_name = cd.name; 
	
	pNew->custom_data = cd.custom_data; 
	pNew->parameters = cd.parameters; 
	pNew->functions = cd.functions; 
	
	pNew->phenotype = cd.phenotype; 
	pNew->is_movable = true;
	pNew->is_out_of_domain = false;
	pNew->displacement.resize(3,0.0); // state? 
	
	pNew->assign_orientation();
	
	return pNew; 
}

/*---------------------------------------------------------------------*/
/* New function for parallel environment: Cell* create_cell(int p_ID)  */
/* Very similar body to the serial Cell* create_cell( void ) above 	 	 */
/* It will call a new Constructor pNew = new Cell(p_ID);							 */
/* This new Constructor will explicitly call a new Constructor of		 	 */
/* the Basic_Agent class : Basic_Agent(int p_ID)											 */ 
/*---------------------------------------------------------------------*/

Cell* create_cell( int p_ID )
{
	Cell* pNew; 
	pNew = new Cell(p_ID);		
	(*all_cells).push_back( pNew ); 
	pNew->index=(*all_cells).size()-1;
	
	
	if( BioFVM::get_default_microenvironment() )
	{
		pNew->register_microenvironment( BioFVM::get_default_microenvironment() );
	} 
	
	return pNew; 
}

void Cell::convert_to_cell_definition( Cell_Definition& cd )
{
	
	// use the cell defaults; 
	type = cd.type; 
	type_name = cd.name; 
	
	custom_data = cd.custom_data; 
	parameters = cd.parameters; 
	functions = cd.functions; 
	
	phenotype = cd.phenotype; 
	// is_movable = true;
	// is_out_of_domain = false;
	
	// displacement.resize(3,0.0); // state? 
	
	assign_orientation();	
	
	return; 
}

void delete_cell( int index )
{
	// released internalized substrates (as of 1.5.x releases)
	(*all_cells)[index]->release_internalized_substrates(); 
	
	// deregister agent in from the agent container
	(*all_cells)[index]->get_container()->remove_agent((*all_cells)[index]);
	// de-allocate (delete) the cell; 
	delete (*all_cells)[index]; 

	// performance goal: don't delete in the middle -- very expensive reallocation
	// alternative: copy last element to index position, then shrink vector by 1 at the end O(constant)

	// move last item to index location  
	(*all_cells)[ (*all_cells).size()-1 ]->index=index;
	(*all_cells)[index] = (*all_cells)[ (*all_cells).size()-1 ];
	// shrink the vector
	(*all_cells).pop_back();	
	return; 
}

void delete_cell( Cell* pDelete )
{
	delete_cell(pDelete->index);
	return; 
}

/*--------------------------------------------------------------------------------------------*/
/* The function vaporize_teleported_cell(int index) is only for distributed parallel scenario */
/* It is exactly like delete(int) function above except for it doesn't call 								  */
/* 'release_internalized_substrates()' as the cell does not die but only leaves the current		*/
/* sub-domain and enters a neighbouring sub-domain 																						*/
/*--------------------------------------------------------------------------------------------*/

void vaporize_teleported_cell(int index)
{
	
	// deregister agent in from the agent container
	/* get the container of the cell and remove the basic agent (basically cell) from it */
	(*all_cells)[index]->get_container()->remove_agent((*all_cells)[index]);
	
	// de-allocate (delete) the cell;
	/* Will it not create a memory leak as there are dynamic members of cell class ?*/ 
	delete (*all_cells)[index]; 

	// performance goal: don't delete in the middle -- very expensive reallocation
	// alternative: copy last element to index position, then shrink vector by 1 at the end O(constant)
	// move last item to index location 
	 
	(*all_cells)[ (*all_cells).size()-1 ]->index=index;
	(*all_cells)[index] = (*all_cells)[ (*all_cells).size()-1 ];
	
	// shrink the vector
	(*all_cells).pop_back();	
	
}

bool is_neighbor_voxel(Cell* pCell, std::vector<double> my_voxel_center, std::vector<double> other_voxel_center, int other_voxel_index)
{
	double max_interactive_distance = pCell->phenotype.mechanics.relative_maximum_adhesion_distance * pCell->phenotype.geometry.radius 
		+ pCell->get_container()->max_cell_interactive_distance_in_voxel[other_voxel_index];
	
	int comparing_dimension = -1, comparing_dimension2 = -1;
	if(my_voxel_center[0] == other_voxel_center[0] && my_voxel_center[1] == other_voxel_center[1])
	{
		comparing_dimension = 2;
	}
	else if(my_voxel_center[0] == other_voxel_center[0] && my_voxel_center[2] == other_voxel_center[2])
	{
		comparing_dimension = 1;
	}
	else if(my_voxel_center[1] == other_voxel_center[1] && my_voxel_center[2] == other_voxel_center[2])
	{
		comparing_dimension = 0;
	}
	
	if(comparing_dimension != -1) 
	{ //then it is an immediate neighbor (through side faces)
		double surface_coord= 0.5*(my_voxel_center[comparing_dimension] + other_voxel_center[comparing_dimension]);
		if(std::fabs(pCell->position[comparing_dimension] - surface_coord) > max_interactive_distance)
		{ return false; }
		return true;
	}
	comparing_dimension=-1;
	
	if(my_voxel_center[0] == other_voxel_center[0])
	{
		comparing_dimension = 1; comparing_dimension2 = 2;
	}
	else if(my_voxel_center[1] == other_voxel_center[1])
	{
		comparing_dimension=0; comparing_dimension2 = 2;
	}
	else if(my_voxel_center[2] == other_voxel_center[2])
	{
		comparing_dimension = 0; comparing_dimension2=1;
	}
	if(comparing_dimension != -1)
	{
		double line_coord1= 0.5*(my_voxel_center[comparing_dimension] + other_voxel_center[comparing_dimension]);
		double line_coord2= 0.5*(my_voxel_center[comparing_dimension2] + other_voxel_center[comparing_dimension2]);
		double distance_squared= std::pow( pCell->position[comparing_dimension] - line_coord1,2)+ std::pow( pCell->position[comparing_dimension2] - line_coord2,2);
		if(distance_squared > max_interactive_distance * max_interactive_distance)
		{ return false; }
		return true;
	}
	std::vector<double> corner_point= 0.5*(my_voxel_center+other_voxel_center);
	double distance_squared= (corner_point[0]-pCell->position[0])*(corner_point[0]-pCell->position[0])
		+(corner_point[1]-pCell->position[1])*(corner_point[1]-pCell->position[1]) 
		+(corner_point[2]-pCell->position[2]) * (corner_point[2]-pCell->position[2]);
	if(distance_squared > max_interactive_distance * max_interactive_distance)
	{ return false; }
	return true;
}

std::vector<Cell*>& Cell::cells_in_my_container( void )
{
	return get_container()->agent_grid[get_current_mechanics_voxel_index()];
}

void Cell::ingest_cell( Cell* pCell_to_eat )
{
	// absorb all the volume(s)

	// absorb fluid volume (all into the cytoplasm) 
	phenotype.volume.cytoplasmic_fluid += pCell_to_eat->phenotype.volume.fluid; 
	pCell_to_eat->phenotype.volume.cytoplasmic_fluid = 0.0; 
	
	// absorb nuclear and cyto solid volume (into the cytoplasm) 
	phenotype.volume.cytoplasmic_solid += pCell_to_eat->phenotype.volume.cytoplasmic_solid; 
	pCell_to_eat->phenotype.volume.cytoplasmic_solid = 0.0; 
	
	phenotype.volume.cytoplasmic_solid += pCell_to_eat->phenotype.volume.nuclear_solid; 
	pCell_to_eat->phenotype.volume.nuclear_solid = 0.0; 
	
	// consistency calculations 
	
	phenotype.volume.fluid = phenotype.volume.nuclear_fluid + 
		phenotype.volume.cytoplasmic_fluid; 
	pCell_to_eat->phenotype.volume.fluid = 0.0; 
	
	phenotype.volume.solid = phenotype.volume.cytoplasmic_solid + 
		phenotype.volume.nuclear_solid; 
	pCell_to_eat->phenotype.volume.solid = 0.0; 
	
	// no change to nuclear volume (initially) 
	pCell_to_eat->phenotype.volume.nuclear = 0.0; 
	pCell_to_eat->phenotype.volume.nuclear_fluid = 0.0; 
	
	phenotype.volume.cytoplasmic = phenotype.volume.cytoplasmic_solid + 
		phenotype.volume.cytoplasmic_fluid; 
	pCell_to_eat->phenotype.volume.cytoplasmic = 0.0; 
	
	phenotype.volume.total = phenotype.volume.nuclear + 
		phenotype.volume.cytoplasmic; 
	pCell_to_eat->phenotype.volume.total = 0.0; 

	phenotype.volume.fluid_fraction = phenotype.volume.fluid / 
		(  phenotype.volume.total ); 
	pCell_to_eat->phenotype.volume.fluid_fraction = 0.0; 

	phenotype.volume.cytoplasmic_to_nuclear_ratio = phenotype.volume.cytoplasmic_solid / 
		( phenotype.volume.nuclear_solid + 1e-16 );
		
	// update corresponding BioFVM parameters (self-consistency) 
	set_total_volume( phenotype.volume.total ); 
	pCell_to_eat->set_total_volume( 0.0 ); 
	
	// absorb the internalized substrates 
	
	// multiply by the fraction that is supposed to be ingested (for each substrate) 
	*(pCell_to_eat->internalized_substrates) *= 
		*(pCell_to_eat->fraction_transferred_when_ingested); // 
	
	*internalized_substrates += *(pCell_to_eat->internalized_substrates); 
	static int n_substrates = internalized_substrates->size(); 
	pCell_to_eat->internalized_substrates->assign( n_substrates , 0.0 ); 	
	
	// trigger removal from the simulation 
	// pCell_to_eat->die(); // I don't think this is safe if it's in an OpenMP loop 
	// flag it for removal 
	pCell_to_eat->flag_for_removal(); 
	// mark it as dead 
	pCell_to_eat->phenotype.death.dead = true; 
	// set secretion and uptake to zero 
	pCell_to_eat->phenotype.secretion.set_all_secretion_to_zero( );  
	pCell_to_eat->phenotype.secretion.set_all_uptake_to_zero( ); 

	
	// deactivate all custom function 
	pCell_to_eat->functions.custom_cell_rule = NULL; 
	pCell_to_eat->functions.update_phenotype = NULL; 
	pCell_to_eat->functions.contact_function = NULL; 
	
	// set it to zero mechanics 
	pCell_to_eat->functions.custom_cell_rule = NULL; 
	
	return; 
}

void Cell::lyse_cell( void )
{
	flag_for_removal(); 
	// mark it as dead 
	phenotype.death.dead = true; 
	// set secretion and uptake to zero 
	phenotype.secretion.set_all_secretion_to_zero( );  
	phenotype.secretion.set_all_uptake_to_zero( ); 
	
	// deactivate all custom function 
	functions.custom_cell_rule = NULL; 
	functions.update_phenotype = NULL; 
	functions.contact_function = NULL; 
	
	// set it to zero mechanics 
	functions.custom_cell_rule = NULL; 

	return; 
}

/*---------------------------------------------*/
/* Added by Gaurav Saxena to print 'most' of	 */
/* the fields of the Cell data structure,			 */
/* useful when comparing cell fields after 		 */
/* pack_and_send() function. 									 */
/*---------------------------------------------*/
void Cell::print_cell(int id)
{
	/*-------------------------------------------------------------------------------------*/
	/* Trying to print all the data fields of the Cell class 															 */
	/* Since Cell is publicly derived from the Basic_Agent class, the protected and	public */
	/* member of Basic_Agent can be accessed in member functions of the Cell class but the */
	/* private members of the Basic_Agent class cannot be accessed by member functions of  */
	/* Derived class. For example: "volume" is not directly accessible but cell_source_sink*/
	/* solver_temp1[i] is directly accessible. 																						 */
	/*-------------------------------------------------------------------------------------*/
	
	static int prnt_counter = 0; 
	if(prnt_counter % 1000 == 0)
	{
	
	std::cout<<std::endl; 
	std::cout<<"=================CELL DATA : "<<prnt_counter<<"======================"; 
	std::cout<<std::endl; 
	
	/*================CELL===================== */

	std::cout<<"type_name:"<<type_name<<std::endl;
	std::cout<<"is_out_of_domain:"<<is_out_of_domain<<std::endl; 
	std::cout<<"is_movable:"<<is_movable<<std::endl;
	std::cout<<"Disp[0]:"<<displacement[0]<<" Disp[1]:"<<displacement[1]<<" Disp[2]:"<<displacement[2]<<std::endl; 
	
	/*===============BASIC_AGENT================ */
	
	std::cout<<endl; 
	std::cout<<"Volume:"<<get_total_volume()<<std::endl;
	std::cout<<"volume_is_changed:"<<get_is_volume_changed()<<std::endl; 
	
	std::cout<<"cell_source_sink_solver_temp1"<<std::endl;
	for(int i=0; i<cell_source_sink_solver_temp1.size();i++)
	 	std::cout<<i<<":"<<cell_source_sink_solver_temp1[i]<<std::endl; 	
	 
	
	std::cout<<"cell_source_sink_solver_temp2"<<std::endl;
	for(int i=0; i<cell_source_sink_solver_temp2.size();i++)
	 	std::cout<<i<<":"<<cell_source_sink_solver_temp2[i]<<std::endl; 	
	
	 
	std::cout<<"Prev_Vel[0]:"<<previous_velocity[0]<<" Prev_Vel[1]:"<<previous_velocity[1]<<" Prev_Vel[2]:"<<previous_velocity[2]<<std::endl;
	std::cout<<"is_active:"<<is_active<<std::endl; 
	
	std::cout<<"total_extracellular_substrate_change"<<std::endl;
	for(int i=0; i<total_extracellular_substrate_change.size();i++)
	 	std::cout<<i<<":"<<total_extracellular_substrate_change[i]<<std::endl; 	
	
	
	std::cout<<"Secretion Rates"<<std::endl;
	for(int i=0; i<secretion_rates->size();i++)
	 	std::cout<<i<<":"<<secretion_rates[i]<<std::endl; 	
	
	
	std::cout<<"Saturation Densities"<<std::endl;
	for(int i=0; i<saturation_densities->size();i++)
	 	std::cout<<i<<":"<<saturation_densities[i]<<std::endl; 	
	
	
	std::cout<<"Uptake Rates"<<std::endl;
	for(int i=0; i<uptake_rates->size();i++)
	 	std::cout<<i<<":"<<uptake_rates[i]<<std::endl; 	
	
	
	std::cout<<"Internalized Substrates"<<std::endl;
	for(int i=0; i<internalized_substrates->size();i++)
	 	std::cout<<i<<":"<<internalized_substrates[i]<<std::endl; 	
	
	
	std::cout<<"Fraction Released At Death"<<std::endl;
	for(int i=0; i<fraction_released_at_death->size();i++)
	 	std::cout<<i<<":"<<fraction_released_at_death[i]<<std::endl; 	
	
	
	std::cout<<"Fraction Transferred When Ingested;"<<std::endl;
	for(int i=0; i<fraction_transferred_when_ingested->size();i++)
	 	std::cout<<i<<":"<<fraction_transferred_when_ingested[i]<<std::endl; 	
	
	
	std::cout<<"ID Number:"<<ID<<std::endl;
	std::cout<<"Type:"<<type<<std::endl; 
	
	std::cout<<"Position[0]: "<<position[0]<<" Position[1]:"<<position[1]<<" Position[2]:"<<position[2]<<std::endl;
	std::cout<<"Vel[0]:"<< velocity[0]<<" Vel[1]:"<<velocity[1]<<" Vel[2]:"<<velocity[2]<<std::endl;
	
	/*==============CUSTOM_CELL_DATA====================*/
	
	std::cout<<endl; 
	std::cout<<"Custom Cell Data"<<std::endl; 
	std::cout<<"Name to Index Map"<<std::endl;
	for(auto it = (this->custom_data.get_name_to_index_map()).cbegin(); it != (this->custom_data.get_name_to_index_map()).cend(); it++)
		std::cout<<"{"<<(*it).first<<":"<<(*it).second<<"}"<<std::endl; 
	
	/*===============CELL_PARAMETERS====================*/
	
	std::cout<<std::endl; 
	std::cout<<"Cell Parameters"<<std::endl;
	std::cout<<"o2_hypoxic_threshold:"<<this->parameters.o2_hypoxic_threshold<<std::endl;
	std::cout<<"o2_hypoxic_response:"<<this->parameters.o2_hypoxic_response<<std::endl; 
	std::cout<<"o2_hypoxic_saturation:"<<this->parameters.o2_hypoxic_saturation<<std::endl;
	std::cout<<"o2_proliferation_saturation:"<<this->parameters.o2_proliferation_saturation<<std::endl;
	std::cout<<"o2_proliferation_threshold"<<this->parameters.o2_proliferation_threshold<<std::endl;
	std::cout<<"o2_reference"<<this->parameters.o2_reference<<std::endl;
	std::cout<<"o2_necrosis_threshold"<<this->parameters.o2_necrosis_threshold<<std::endl;
	std::cout<<"o2_necrosis_max"<<this->parameters.o2_necrosis_max<<std::endl;
	std::cout<<"max_necrosis_rate"<<this->parameters.max_necrosis_rate<<std::endl;
	std::cout<<"necrosis_type"<<this->parameters.necrosis_type<<std::endl;
	if(this->parameters.pReference_live_phenotype != NULL)
		std::cout<<"pReference_live_phenotype pointer in Cell_Parameters is NOT NULL"<<std::endl;
	if(this->parameters.pReference_live_phenotype == &phenotype)
		std::cout<<"pReference_live_phenotype pointer in Cell_Parameters is pointing to object phenotype of class Phenotype"<<std::endl; 
	
	/*=========CELL_STATE============*/
	
	std::cout<<std::endl; 
	std::cout<<"Orientation"<<std::endl;
	std::cout<<"Orient[0]:"<< (this->state).orientation[0] << " Orient[1]:" << (this->state).orientation[1] << " Orient[2]:" << (this->state).orientation[2] << std::endl; 
	std::cout<<"Simple Pressure:"<<this->state.simple_pressure<<std::endl; 
	
	/*=========PHENOTYPE===============*/
	
	std::cout<<std::endl; 
	std::cout<<"Phenotype Data"<<std::endl;
	std::cout<<"flagged_for_division:"<<this->phenotype.flagged_for_division<<std::endl;
	std::cout<<"flagged_for_removal:"<<this->phenotype.flagged_for_removal<<std::endl; 
	
	
	/*===========VARIABLE==============*/
	
	std::cout<<std::endl; 
	std::cout<<"Custom Cell Data variables"<<std::endl;
	for(int i=0; i<this->custom_data.variables.size(); i++)
		std::cout<<"i:"<<i<<" name:"<<this->custom_data.variables[i].name<<" value:"<<this->custom_data.variables[i].value << " units:"<<this->custom_data.variables[i].units<<std::endl;  
	
	/*=========VECTOR_VARIABLE==========*/
	
	std::cout<<std::endl; 
	std::cout<<"Custom Cell Data Variables Vectors"<<std::endl; 
	for(int i=0; i<this->custom_data.vector_variables.size();i++)
	{
		std::cout<<"i:"<<i<<" name:"<<this->custom_data.vector_variables[i].name<<" units:"<<this->custom_data.variables[i].units<<std::endl;
		for(int j=0; j<this->custom_data.vector_variables[i].value.size();j++)
			std::cout<<" value "<<j<<":"<<this->custom_data.vector_variables[i].value[j]<<std::endl;
	}
	
	/*===========CYCLE_MODEL=============*/
	
	std::cout<<std::endl; 
	std::cout<<"Cycle Model data"<<std::endl;
	std::cout<<"name:"<<this->functions.cycle_model.name<<" code:"<<this->functions.cycle_model.code<<std::endl; 
	std::cout<<"default_phase_index:"<<this->functions.cycle_model.default_phase_index<<std::endl; 
	std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_model = this->functions.cycle_model.get_inverse_index_maps();
	std::cout<<"Vector of unordered maps"<<std::endl;  
	for(int i=0; i<inverse_index_maps_cycle_model.size(); i++)
	{
		for(auto it = inverse_index_maps_cycle_model[i].cbegin(); it != inverse_index_maps_cycle_model[i].cend(); it++)
			std::cout<<"{"<<(*it).first<<":"<<(*it).second<<"}"<<std::endl;
		std::cout<<endl; 
	}
	
	/*==============CYCLE=================*/
	
	std::cout<<std::endl; 
	std::cout<<"Cycle data - only pointer checking"<<std::endl;
	if(this->phenotype.cycle.pCycle_Model != NULL)
		std::cout<<"Pointer pCycle_Model in Cycle class is NOT NULL"<<std::endl; 
	
	
	/*===============DEATH=================*/
	
	std::cout<<std::endl; 
	std::cout<<"Death Data"<<std::endl;
	std::cout<<"dead:"<<this->phenotype.death.dead<<std::endl; 
	std::cout<<"current_death_model_index:"<<this->phenotype.death.current_death_model_index<<std::endl; 
	for(int i=0; i < this->phenotype.death.rates.size(); i++)
		std::cout<<"rates "<<i<<":"<<this->phenotype.death.rates[i]<<std::endl;
	std::cout<<"Checking if vector of Cycle_Model pointers contains NULL entries"<<std::endl;
	for(int i=0; i< this->phenotype.death.models.size(); i++)
		if(this->phenotype.death.models[i] != NULL)
			cout<<"models["<<i<<"]"<<" is NOT NULL"<<std::endl; 
		
	/*===============VOLUME================*/
	
	std::cout<<std::endl; 
	std::cout<<"Volume Data"<<std::endl;
	std::cout<<"total:"<<this->phenotype.volume.total<<std::endl; 
	std::cout<<"solid:"<<this->phenotype.volume.solid<<std::endl;
	std::cout<<"fluid:"<<this->phenotype.volume.fluid<<std::endl;
	std::cout<<"fluid_fraction:"<<this->phenotype.volume.fluid_fraction<<std::endl;
	std::cout<<"nuclear:"<<this->phenotype.volume.nuclear<<std::endl;
	std::cout<<"nuclear_fluid:"<<this->phenotype.volume.nuclear_fluid<<std::endl; 
	std::cout<<"nuclear_solid:"<<this->phenotype.volume.nuclear_solid<<std::endl;
	std::cout<<"cytoplasmic:"<<this->phenotype.volume.cytoplasmic<<std::endl;
	std::cout<<"cytoplasmic_fluid:"<<this->phenotype.volume.cytoplasmic_fluid<<std::endl;
	std::cout<<"cytoplasmic_solid:"<<this->phenotype.volume.cytoplasmic_solid<<std::endl;
	std::cout<<"calcified_fraction:"<<this->phenotype.volume.calcified_fraction<<std::endl; 
	std::cout<<"cytoplasmic_to_nuclear_ratio:"<<this->phenotype.volume.cytoplasmic_to_nuclear_ratio<<std::endl;
	std::cout<<"rupture_volume:"<<this->phenotype.volume.rupture_volume<<std::endl;
	std::cout<<"cytoplasmic_biomass_change_rate:"<<this->phenotype.volume.cytoplasmic_biomass_change_rate<<std::endl;
	std::cout<<"nuclear_biomass_change_rate:"<<this->phenotype.volume.nuclear_biomass_change_rate<<std::endl;
	std::cout<<"fluid_change_rate:"<<this->phenotype.volume.fluid_change_rate<<std::endl; 
	std::cout<<"calcification_rate:"<<this->phenotype.volume.calcification_rate<<std::endl;
	std::cout<<"target_solid_cytoplasmic:"<<this->phenotype.volume.target_solid_cytoplasmic<<std::endl;
	std::cout<<"target_solid_nuclear:"<<this->phenotype.volume.target_solid_nuclear<<std::endl;
	std::cout<<"target_fluid_fraction:"<<this->phenotype.volume.target_fluid_fraction<<std::endl;
	std::cout<<"target_cytoplasmic_to_nuclear_ratio:"<<this->phenotype.volume.target_cytoplasmic_to_nuclear_ratio<<std::endl;
	std::cout<<"relative_rupture_volume:"<<this->phenotype.volume.relative_rupture_volume<<std::endl;

	/*===============GEOMETRY============*/
	
	std::cout<<std::endl; 
	std::cout<<"Geometry Data"<<std::endl;
	std::cout<<"radius:"<<this->phenotype.geometry.radius<<std::endl;
	std::cout<<"nuclear_radius:"<<this->phenotype.geometry.nuclear_radius<<std::endl;
	std::cout<<"surface_area:"<<this->phenotype.geometry.surface_area<<std::endl;
	std::cout<<"polarity:"<<this->phenotype.geometry.polarity<<std::endl;
	
	/*===============MECHANICS============*/
	
	std::cout<<std::endl; 
	std::cout<<"Mechanics Data"<<std::endl;
	std::cout<<"cell_adhesion_strength:"<<this->phenotype.mechanics.cell_cell_adhesion_strength<<std::endl;
	std::cout<<"cell_BM_adhesion_strength:"<<this->phenotype.mechanics.cell_BM_adhesion_strength<<std::endl;
	std::cout<<"cell_cell_repulsion_strength:"<<this->phenotype.mechanics.cell_cell_repulsion_strength<<std::endl;
	std::cout<<"cell_BM_repulsion_strength:"<<this->phenotype.mechanics.cell_BM_repulsion_strength<<std::endl;	 
	std::cout<<"relative_maximum_adhesion_distance:"<<this->phenotype.mechanics.relative_maximum_adhesion_distance<<std::endl;
	
	/*=================MOTILITY===========*/
	
	std::cout<<std::endl; 
	std::cout<<"Motility Data"<<std::endl;
	std::cout<<"is_motile:"<<this->phenotype.motility.is_motile<<std::endl;
	std::cout<<"persistence_time:"<<this->phenotype.motility.persistence_time<<std::endl;
	std::cout<<"migration_speed:"<<this->phenotype.motility.migration_speed<<std::endl;
	std::cout<<"migration_bias:"<<this->phenotype.motility.migration_bias<<std::endl;
	std::cout<<"restrict_to_2D:"<<this->phenotype.motility.restrict_to_2D<<std::endl;
	for(int i=0; i<this->phenotype.motility.migration_bias_direction.size(); i++)
		std::cout<<"bias direction "<<i<<":"<<this->phenotype.motility.migration_bias_direction[i]<<std::endl;
	for(int i=0; i<this->phenotype.motility.motility_vector.size(); i++)
		std::cout<<"motility vector "<<i<<":"<<this->phenotype.motility.motility_vector[i]<<std::endl; 
	
	/*=================SECRETION============*/
	
	std::cout<<std::endl; 
	std::cout<<"Secretion Data"<<std::endl;
	for(int i=0; i<this->phenotype.secretion.secretion_rates.size();i++)
			std::cout<<"secretion_rates "<<i<<":"<<this->phenotype.secretion.secretion_rates[i]<<std::endl; 
	for(int i=0; i<this->phenotype.secretion.uptake_rates.size();i++)
			std::cout<<"uptake_rates "<<i<<":"<<this->phenotype.secretion.uptake_rates[i]<<std::endl;
	for(int i=0; i<this->phenotype.secretion.saturation_densities.size();i++)
			std::cout<<"saturation densities "<<i<<":"<<this->phenotype.secretion.saturation_densities[i]<<std::endl;
			
	/*=================MOLECULAR============*/
		
	std::cout<<std::endl; 
	std::cout<<"Molecular Data"<<std::endl;
	for(int i=0; i<this->phenotype.molecular.internalized_total_substrates.size();i++)
			std::cout<<"internalized_total_substrates "<<i<<":"<<this->phenotype.molecular.internalized_total_substrates[i]<<std::endl; 
	for(int i=0; i<this->phenotype.molecular.fraction_released_at_death.size();i++)
			std::cout<<"fraction_released_at_death "<<i<<":"<<this->phenotype.molecular.fraction_released_at_death[i]<<std::endl;
	for(int i=0; i<this->phenotype.molecular.fraction_transferred_when_ingested.size();i++)
			std::cout<<"fraction_transferred_when_ingested "<<i<<":"<<this->phenotype.molecular.fraction_transferred_when_ingested[i]<<std::endl;	
		
	/*===================PHASE===============*/
	
	std::cout<<std::endl; 
	std::cout<<"Phase Data"<<std::endl;
	for(int i=0; i<this->functions.cycle_model.phases.size();i++)
	{
		std::cout<<i<<":"<<std::endl;
		std::cout<<"index:"<<this->functions.cycle_model.phases[i].index<<std::endl;
		std::cout<<"code:"<<this->functions.cycle_model.phases[i].code<<std::endl; 
		std::cout<<"name"<<this->functions.cycle_model.phases[i].name<<std::endl;
		std::cout<<"division_at_phase_exit"<<this->functions.cycle_model.phases[i].division_at_phase_exit<<std::endl;
		std::cout<<"removal_at_phase_exit"<<this->functions.cycle_model.phases[i].removal_at_phase_exit<<std::endl; 
	}
	
	/*================PHASE_LINKS==============*/
	
	std::cout<<std::endl; 
	std::cout<<"Phase Links Data"<<std::endl;
	for(int i=0; i<this->functions.cycle_model.phase_links.size();i++)
		{
			std::cout<<"Vector:"<<i<<std::endl;
			for(int j=0; j<this->functions.cycle_model.phase_links[i].size();j++)
			{
				std::cout<<"Vector:"<<j<<std::endl;
				std::cout<<"start_phase_index"<<this->functions.cycle_model.phase_links[i][j].start_phase_index<<std::endl;
				std::cout<<"end_phase_index"<<this->functions.cycle_model.phase_links[i][j].end_phase_index<<std::endl;
				std::cout<<"fixed_duration"<<this->functions.cycle_model.phase_links[i][j].fixed_duration<<std::endl;
			}
		}
	
	/*==================CYCLE_DATA==============*/
	
	std::cout<<std::endl; 
	std::cout<<"Cycle_Data"<<std::endl;
	std::cout<<"time_units:"<<this->phenotype.cycle.data.time_units<<std::endl;
	std::cout<<"current_phase_index:"<<this->phenotype.cycle.data.current_phase_index<<std::endl;
	std::cout<<"elapsed_time_in_phase:"<<this->phenotype.cycle.data.elapsed_time_in_phase<<std::endl;
	for(int i=0; i<this->phenotype.cycle.data.transition_rates.size();i++)
		{
			std::cout<<"Transition Rates Vector:"<<i<<std::endl;
			for(int j=0; j<this->phenotype.cycle.data.transition_rates[i].size();j++)
				std::cout<<"rate "<<j<<":"<<this->phenotype.cycle.data.transition_rates[i][j]<<std::endl; 
		}
	std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_data = this->phenotype.cycle.data.get_inverse_index_maps();
	std::cout<<"Vector of unordered maps";  
	for(int i=0; i<inverse_index_maps_cycle_data.size(); i++)
	{
		for(auto it = inverse_index_maps_cycle_data[i].cbegin(); it != inverse_index_maps_cycle_data[i].cend(); it++)
			std::cout<<"{"<<(*it).first<<":"<<(*it).second<<"}"<<std::endl;
		std::cout<<endl; 
	}
	std::cout<<"Checking if pCycle_Model pointer is NULL"<<std::endl;
	if(this->phenotype.cycle.data.pCycle_Model != NULL)
		std::cout<<"pCycle_Model is NOT NULL"; 
		
	/*=============DEATH_PARAMETERS===============*/
		
	std::cout<<std::endl; 
	std::cout<<"Death Parameters"<<std::endl;
	for(int i=0; i<this->phenotype.death.parameters.size(); i++)
		{	
			std::cout<<i<<":"<<std::endl;
			std::cout<<"time_units:"<<this->phenotype.death.parameters[i].time_units<<std::endl;
			std::cout<<"unlysed_fluid_change_rate:"<<this->phenotype.death.parameters[i].unlysed_fluid_change_rate<<std::endl;
			std::cout<<"lysed_fluid_change_rate:"<<this->phenotype.death.parameters[i].lysed_fluid_change_rate<<std::endl;
			std::cout<<"cytoplasmic_biomass_change_rate:"<<this->phenotype.death.parameters[i].cytoplasmic_biomass_change_rate<<std::endl;
			std::cout<<"nuclear_biomass_change_rate:"<<this->phenotype.death.parameters[i].nuclear_biomass_change_rate<<std::endl;
			std::cout<<"calcification rate:"<<this->phenotype.death.parameters[i].calcification_rate<<std::endl;
			std::cout<<"relative_rupture_volume:"<<this->phenotype.death.parameters[i].relative_rupture_volume<<std::endl;
		}	
	}
	prnt_counter++; 
}


/*---------------------------------------------*/
/* Added by Gaurav Saxena to pack almost all 	 */
/* fields of Cell data structure when cells		 */
/* cross sub-domain boundaries.								 */
/*---------------------------------------------*/

void Cell_Container::pack(std::vector<Cell*> *all_cells, mpi_Environment &world, mpi_Cartesian &cart_topo)
{

/*-----------------------------------------------------------------------------*/	
/* Right now test for a single cell. Then later (1) Loop through all cells		 */
/* (2) Test if crossed_to_left_subdomain == true or crossed_to_right_subdomain */
/* == true. Only one can be true. (3) If anyone is true pack cells contents 	 */
/* in snd_buf_to_left or snd_buf_to_right data members, respectively. 				 */
/*-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/* Temporary help variables to contain lengths and temporarily variables */
/* Need to maintain two positions: one for left buffer and one for right */
/* buffer. 
/*-----------------------------------------------------------------------*/	

	
	
	int len_snd_buf_left  	 = 0; 
	int len_snd_buf_right 	 = 0; 
	
	 
	
	int 	 																	temp_int;
	double 																	temp_double;
	std::string 														temp_str;
	std::unordered_map<std::string, int> :: iterator it;
	Cell 																		*pCell; 
	  
	
	int len_int 						 = 0; 
	int len_double 					 = 0; 
	int len_str 						 = 0; 
	int len_vector					 = 0; 
	
	/*---------------------------------------------------------*/
	/* All 6 members below are now in the Cell_Container class */
	/*---------------------------------------------------------*/
	
	position_left 			 		= 0;					//Must be initialized to 0
	position_right 			 		= 0;					//Must be initialized to 0
	no_cells_cross_left  		= 0;					//Must be initialized to 0
	no_cells_cross_right		= 0;					//Must be initialized to 0
	no_of_cells_from_right 	= 0;
	no_of_cells_from_left 	= 0; 
	
	snd_buf_left.resize(0); 					//When we enter function again, this is reset
	snd_buf_right.resize(0);					//Reset the others also
	rcv_buf_left.resize(0);
	rcv_buf_right.resize(0); 

	
	/*--------------------------------------------------------------------------------*/
	/* First count the number of cells crossing to left subdomain and right subdomain */
	/*--------------------------------------------------------------------------------*/
	
	for(int i=0; i<(*all_cells).size();i++)
	{
		pCell = (*all_cells)[i]; 		
		if(pCell->crossed_to_left_subdomain == true)
			no_cells_cross_left++; 
		if(pCell->crossed_to_right_subdomain == true)
			no_cells_cross_right++; 		
	}

	/*-----------------------------------------------------------------------------------*/	
	/* There is no need to pack crossed_to_left_subdomain and crossed_to_right_subdomain */
	/* as we send these two fields as the First Communication between MPI processes. 		 */
	/* Also try to send the length of the buffer to be allocated when sending.					 */
	/*-----------------------------------------------------------------------------------*/	
	
  //std::cout<<"Total cells crossing to left in Rank "<<world.rank<<":"<<no_cells_cross_left<<std::endl;
	//std::cout<<"Total cells crossing to right in Rank "<<world.rank<<":"<<no_cells_cross_right<<std::endl;

	
	/* IMPORTANT: CANNOT USE #pragma omp for HERE AS ALL THREADS WILL WRITE TO THE SAME SHARED BUFFER */
	for(int i=0; i<(*all_cells).size();i++)
	{
		pCell = (*all_cells)[i]; 
					
		if(pCell->crossed_to_left_subdomain == true)
		{			
			//pCell->crossed_to_left_subdomain = false; //RESET IT BUT LATER REMOVE THIS LINE
			/* Cell ID first - needed to create a cell */
						
			len_snd_buf_left = position_left + sizeof(pCell->ID);
			snd_buf_left.resize(len_snd_buf_left);
			MPI_Pack(&(pCell->ID), 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
			
			/* Pack double position[0], position[1], position[2] - needed to assign position 							 */
			/* This is the only array where I am not packing the length of the array with the array values */
			/* as we're only dealing with 3-dimensions here, hence position has x, y, z values						 */
			
			len_snd_buf_left = position_left + 3 * sizeof(pCell->position[0]); 
			snd_buf_left.resize(len_snd_buf_left); 
			MPI_Pack(&(pCell->position[0]), 3, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			
			/* Pack length of string (MPI_INT) and then string pCell->type_name (MPI_CHAR) */
			
			len_str = pCell->type_name.length(); 
			len_snd_buf_left = position_left + sizeof(len_str) + len_str;
			snd_buf_left.resize(len_snd_buf_left); 
			MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
			MPI_Pack(&(pCell->type_name[0]), len_str , MPI_CHAR, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
			
			/* Packing unordered_map<std::string, int> in Custom_Data class */
			/* First pack the number of entries in the unordered_map 				*/
			
			len_int = pCell->custom_data.get_name_to_index_map().size();
			len_snd_buf_left = position_left + sizeof(len_int);
			snd_buf_left.resize(len_snd_buf_left);  
			MPI_Pack(&len_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);  
			
			/* Now pack the <string, int> pairs but store length(string) before string */
					
			for(auto it = (pCell->custom_data.get_name_to_index_map()).cbegin(); it != (pCell->custom_data.get_name_to_index_map()).cend(); it++)
			{
				/* resize buffer to contain length of string, string and integer */
				
				temp_str = it->first; 
				temp_int = it->second; 
				len_str = temp_str.length();
				len_snd_buf_left = position_left + sizeof(len_str) + len_str + sizeof(temp_int); 
				snd_buf_left.resize(len_snd_buf_left);
				 
				MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
				MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
				MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD ); 			
			}	
			
		 /* Packing members of class 'Variable' (as we are going depth first i.e. Inorder traversal of data members) */
		 
		 	for(int i=0; i<pCell->custom_data.variables.size(); i++)
		 	{
		    temp_str = pCell->custom_data.variables[i].name; 
		    len_str = temp_str.length();
		    len_snd_buf_left = position_left + len_str + sizeof(len_str); 
		    snd_buf_left.resize(len_snd_buf_left);
		    MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		    MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		    
		    /* There was a mistake in MPI_Pack here, instead of &snd_buf_left[0], I had &snd_buf_left */
		    temp_double = pCell->custom_data.variables[i].value;
		    len_snd_buf_left = position_left + sizeof(temp_double); 
		    MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		    
		    temp_str = pCell->custom_data.variables[i].units; 
		    len_str = temp_str.length();
		    len_snd_buf_left = position_left + len_str + sizeof(len_str); 
		    snd_buf_left.resize(len_snd_buf_left);
		    MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		    MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);		     
			}
			
			/* Packing members of class 'Vector_Variable' i.e. string name, a vector<double> value and string units */
			
			for(int i=0; i<pCell->custom_data.vector_variables.size(); i++)
			{
				temp_str = pCell->custom_data.vector_variables[i].name; 
		    len_str = temp_str.length();
		    len_snd_buf_left = position_left + len_str + sizeof(len_str); 
		    snd_buf_left.resize(len_snd_buf_left);
		    MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		    MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		    
		    len_vector = pCell->custom_data.vector_variables[i].value.size(); 
		    len_snd_buf_left = position_left + len_vector * sizeof(double);
		    snd_buf_left.resize(len_snd_buf_left); 
		    MPI_Pack(&(pCell->custom_data.vector_variables[i].value[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		    
		    temp_str = pCell->custom_data.vector_variables[i].units; 
		    len_str = temp_str.length();
		    len_snd_buf_left = position_left + len_str + sizeof(len_str); 
		    snd_buf_left.resize(len_snd_buf_left);
		    MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		    MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);  	    
			}
			
			/* Packing members of class 'Cell_Parameters' except for pReference_live_phenotype pointer */
			/* There are 9 doubles and 1 int in this class except for the pointer */
			
			/*=====================================================================*/
			/* pReference_live_phenotype * IS NOT PACKED 													 */
			/*=====================================================================*/
			
			len_snd_buf_left = position_left + 9 * sizeof(double) + 1 * sizeof(int); 
			snd_buf_left.resize(len_snd_buf_left); 
			
			MPI_Pack(&(pCell->parameters.o2_hypoxic_threshold), 				1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.o2_hypoxic_response), 					1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.o2_hypoxic_saturation), 				1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.o2_proliferation_saturation), 	1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.o2_proliferation_threshold), 	1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.o2_reference), 								1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.o2_necrosis_threshold), 				1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.o2_necrosis_max), 							1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.max_necrosis_rate), 						1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.necrosis_type), 								1, MPI_INT, 	 &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		
		
		
		/* Packing data members of Cell_Functions in depth first manner */
		/* It only has Cycle_Model data member, so pack Cycle_Model data members depth first */
		
		std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_model = pCell->functions.cycle_model.get_inverse_index_maps();
		len_vector = inverse_index_maps_cycle_model.size(); 
		len_snd_buf_left = position_left + sizeof(len_vector); 
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		
		for(int i=0; i<len_vector; i++)
		{
			 len_int = inverse_index_maps_cycle_model[i].size();
			 len_snd_buf_left = position_left + sizeof(len_int) + len_int * 2 * sizeof(int); 
			 snd_buf_left.resize(len_snd_buf_left);
			 MPI_Pack(&len_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			 
			for(auto it = inverse_index_maps_cycle_model[i].cbegin(); it != inverse_index_maps_cycle_model[i].cend(); it++)
			{	 
				MPI_Pack(&(it->first),  1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
				MPI_Pack(&(it->second), 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD ); 			
			}				 
		}
		
		temp_str = pCell->functions.cycle_model.name; 
		len_str  = temp_str.length();
		len_snd_buf_left = position_left + sizeof(len_str) + len_str; 
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		
		len_snd_buf_left = position_left + sizeof(int);
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&(pCell->functions.cycle_model.code), 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		
		/* Packing of data members of Phases as its object is a member of Cycle_Model */
		/* std::vector<Phase> phases; is a member of it - its a vector */		
		
		len_vector = pCell->functions.cycle_model.phases.size(); 
		len_snd_buf_left = position_left + sizeof(len_vector);
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		
		for(int i=0; i<len_vector; i++)
		{
			len_snd_buf_left = position_left + 2 * sizeof(int);
			snd_buf_left.resize(len_snd_buf_left); 
			 
			temp_int = pCell->functions.cycle_model.phases[i].index; 
			MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
			temp_int = pCell->functions.cycle_model.phases[i].code;
			MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
			
			temp_str = pCell->functions.cycle_model.phases[i].name; 
			len_str  = temp_str.length();
			len_snd_buf_left = position_left + sizeof(len_str) + len_str; 
			snd_buf_left.resize(len_snd_buf_left); 
			MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			
			/* There are 2 remaining bool values which I will pack in 2 ints */
			len_snd_buf_left = position_left + 2 * sizeof(int);
			snd_buf_left.resize(len_snd_buf_left);
			
			if(pCell->functions.cycle_model.phases[i].division_at_phase_exit == true)
				temp_int = 1; 
			else
				temp_int = 0;
			MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
			
			if(pCell->functions.cycle_model.phases[i].removal_at_phase_exit == true)
				temp_int = 1; 
			else
				temp_int = 0;
			MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);		
		}
		
		/* Packing members of Phase_Link - std::vector<std::vector<Phase_Link>> phase_links is a member of */
		/* class Cycle_Model which is a member of Functions - we're packing in DFS i.e. Inorder traversal	 */
		
		len_vector = pCell->functions.cycle_model.phase_links.size();
		len_snd_buf_left = position_left + sizeof(len_vector);
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		for(int i=0; i<len_vector; i++)
		{
			len_int = pCell->functions.cycle_model.phase_links[i].size();
			len_snd_buf_left = position_left + sizeof(len_int) + len_int * 3 * sizeof(int); 
			snd_buf_left.resize(len_snd_buf_left); 
			MPI_Pack(&len_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			for(int j=0; j<len_int; j++)
			{
				temp_int = pCell->functions.cycle_model.phase_links[i][j].start_phase_index; 
				MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
				temp_int = pCell->functions.cycle_model.phase_links[i][j].end_phase_index;
				MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
				if(pCell->functions.cycle_model.phase_links[i][j].fixed_duration == true)
					temp_int = 1;
				else
					temp_int = 0; 
				MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
			}
		}
		
		/* Going back to packing the remaning data fields of Cycle_Model */
		/* Next field is is 'default_phase_index' of Cycle_Model class */
		
		temp_int = pCell->functions.cycle_model.default_phase_index; 
		len_snd_buf_left = position_left + sizeof(temp_int);
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		
		 /* Packing data fields of Cycle_Data as it is an object of Cycle_Model class above */
		 /* First field is std::vector<std::unordered_map<int,int>> inverse_index_maps;			*/
		 
		std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_data = pCell->functions.cycle_model.data.get_inverse_index_maps();	 
		len_vector = inverse_index_maps_cycle_data.size(); 
		len_snd_buf_left = position_left + sizeof(len_vector); 
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		
		for(int i=0; i<len_vector; i++)
		{
			 len_int = inverse_index_maps_cycle_data[i].size();
			 len_snd_buf_left = position_left + sizeof(len_int) + len_int * 2 * sizeof(int); 
			 snd_buf_left.resize(len_snd_buf_left);
			 MPI_Pack(&len_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			 
			for(auto it = inverse_index_maps_cycle_data[i].cbegin(); it != inverse_index_maps_cycle_data[i].cend(); it++)
			{	 
				MPI_Pack(&(it->first),  1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
				MPI_Pack(&(it->second), 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD ); 			
			}				 
		}
		
			/*=====================================================================*/
			/* Cycle_Model* pCycle_Model;  IS NOT PACKED 													 */
			/*=====================================================================*/
			
			temp_str = pCell->functions.cycle_model.data.time_units;
			len_str = temp_str.length();
			len_snd_buf_left = position_left + sizeof(len_str) + len_str; 
			snd_buf_left.resize(len_snd_buf_left);
			MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		 
		 len_vector = pCell->functions.cycle_model.data.transition_rates.size(); 
		 len_snd_buf_left = position_left + sizeof(len_vector);
		 snd_buf_left.resize(len_snd_buf_left);
		 MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		 
		 for(int i=0; i<len_vector; i++)
		 {
		 	len_int = pCell->functions.cycle_model.data.transition_rates[i].size();
		 	len_snd_buf_left = position_left + sizeof(len_int) + len_int * sizeof(double);
		 	snd_buf_left.resize(len_snd_buf_left);
		 	MPI_Pack(&len_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		 	
		 	for(int j=0; j<len_int; j++)
		 	{
		 		temp_double = pCell->functions.cycle_model.data.transition_rates[i][j];
		 		MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		 	}
		 }
		 
		 /* 1 int and 1 double are next to be packed in class Cycle_Data */
		 
		 len_snd_buf_left = position_left + sizeof(int) + sizeof(double);
		 snd_buf_left.resize(len_snd_buf_left);
		 MPI_Pack(&(pCell->functions.cycle_model.data.current_phase_index), 	1, MPI_INT, 	 &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		 MPI_Pack(&(pCell->functions.cycle_model.data.elapsed_time_in_phase), 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 

		/* Now packing members of class Cell_State as this is the next member in class Cell */
		/* has std::vector<double> orientation and double simple_pressure 									*/
		
		len_vector = pCell->state.orientation.size();
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double) + sizeof(double); 
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		MPI_Pack(&(pCell->state.orientation[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);		
		MPI_Pack(&(pCell->state.simple_pressure), 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		
		/* Packing data members of class Phenotype next */
		
		/* First it has 2 boolean variables - allocate sizeof(int) space for both of them */
		len_snd_buf_left = position_left + 2 * sizeof(int);
		snd_buf_left.resize(len_snd_buf_left);
		if(pCell->phenotype.flagged_for_division == true)
			temp_int = 1;
		else
			temp_int = 0;
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		if(pCell->phenotype.flagged_for_removal == true)
			temp_int = 1;
		else
			temp_int = 0;
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		/* Packing of members of class Cycle as the class Phenotype contains a data member of this class */
		/* Not packing pCell->phenotype.cycle.pCycle_Model pointer data member *
		
			/*=====================================================================*/
			/* Cycle_Model* pCycle_Model;  IS NOT PACKED 													 */
			/*=====================================================================*/
		
		 /* Packing data fields of Cycle_Data as it is an object of Cycle class which is    */
		 /* a data member of the Phenotype class which is a member of class Cell						*/
		 /* First field is std::vector<std::unordered_map<int,int>> inverse_index_maps;			*/
		 
		std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_data_1 = pCell->phenotype.cycle.data.get_inverse_index_maps();	 
		len_vector = inverse_index_maps_cycle_data_1.size(); 
		len_snd_buf_left = position_left + sizeof(len_vector);
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		
		for(int i=0; i<len_vector; i++)
		{
			 len_int = inverse_index_maps_cycle_data_1[i].size();
			 len_snd_buf_left = position_left + sizeof(len_int) + len_int * 2 * sizeof(int); 
			 snd_buf_left.resize(len_snd_buf_left);
			 MPI_Pack(&len_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			 
			for(auto it = inverse_index_maps_cycle_data_1[i].cbegin(); it != inverse_index_maps_cycle_data_1[i].cend(); it++)
			{	 
				MPI_Pack(&(it->first),  1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
				MPI_Pack(&(it->second), 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD ); 			
			}				 
		}
		
			/*=====================================================================*/
			/* Cycle_Model* pCycle_Model;  IS NOT PACKED  												 */
			/*=====================================================================*/
			
			temp_str = pCell->phenotype.cycle.data.time_units;
			len_str = temp_str.length();
			len_snd_buf_left = position_left + sizeof(len_str) + len_str; 
			snd_buf_left.resize(len_snd_buf_left);
			MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		 
		 len_vector = pCell->phenotype.cycle.data.transition_rates.size(); 
		 len_snd_buf_left = position_left + sizeof(len_vector);
		 snd_buf_left.resize(len_snd_buf_left);
		 MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		 
		 for(int i=0; i<len_vector; i++)
		 {
		 	len_int = pCell->phenotype.cycle.data.transition_rates[i].size();
		 	len_snd_buf_left = position_left + sizeof(len_int) + len_int * sizeof(double);
		 	snd_buf_left.resize(len_snd_buf_left);
		 	MPI_Pack(&len_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		 	
		 	for(int j=0; j<len_int; j++)
		 	{
		 		temp_double = pCell->phenotype.cycle.data.transition_rates[i][j];
		 		MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		 	}
		 }
		 
		 /* 1 int and 1 double are next to be packed in class Cycle_Data */
		 
		 len_snd_buf_left = position_left + sizeof(int) + sizeof(double);
		 snd_buf_left.resize(len_snd_buf_left);
		 MPI_Pack(&(pCell->phenotype.cycle.data.current_phase_index), 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		 MPI_Pack(&(pCell->phenotype.cycle.data.elapsed_time_in_phase), 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		/* Now packing data members of class Death as there is an object of class Death in Phenotype class */
		
		len_vector = pCell->phenotype.death.rates.size();
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double); 
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.death.rates[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		
		/* Next field is std::vector<Cycle_Model*> models ;  this is a vector of Cycle_Model * pointers */
		
		/*=====================================================================*/
		/* std::vector<Cycle_Model*> models ;  IS NOT PACKED  								 */
		/*=====================================================================*/
		
		/* Now packing members of class Death_Parameters as class Death contains */
		/* std::vector<Death_Parameters> parameters 														 */
		
		len_vector = pCell->phenotype.death.parameters.size();
		len_snd_buf_left = position_left + sizeof(len_vector);
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		 
		for(int i=0; i<len_vector; i++)
		{
			temp_str = pCell->phenotype.death.parameters[i].time_units;
			len_str  = temp_str.length();
			len_snd_buf_left = position_left + sizeof(len_str) + len_str; 
			snd_buf_left.resize(len_snd_buf_left); 
			MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
			MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
			
			/* 6 double data types are to be packed now */
			len_snd_buf_left = position_left + 6 * sizeof(double);
			snd_buf_left.resize(len_snd_buf_left);
			
			temp_double = pCell->phenotype.death.parameters[i].unlysed_fluid_change_rate;
			MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
			temp_double = pCell->phenotype.death.parameters[i].lysed_fluid_change_rate;
			MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);			
			temp_double = pCell->phenotype.death.parameters[i].cytoplasmic_biomass_change_rate;
			MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
			temp_double = pCell->phenotype.death.parameters[i].nuclear_biomass_change_rate;
			MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
			temp_double = pCell->phenotype.death.parameters[i].calcification_rate;
			MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
			temp_double = pCell->phenotype.death.parameters[i].relative_rupture_volume;
			MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		}
		 
		 /* Next 1 bool value (allocate sizeof(int)) and 1 int value are to be packed */
		 
		 len_snd_buf_left = position_left + 2 * sizeof(int);
		 snd_buf_left.resize(len_snd_buf_left);
		 
		 if(pCell->phenotype.death.dead == true)
		 		temp_int = 1;
		 else
		 		temp_int = 0; 
		 MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		 
		 temp_int = pCell->phenotype.death.current_death_model_index;
		 MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		/* Next data members of class Volume are packed as class Phenotype contains its object */
		/* 22 doubles are to be packed */
		
		len_snd_buf_left = position_left + 22 * sizeof(double); 
		snd_buf_left.resize(len_snd_buf_left);
		
		MPI_Pack(&(pCell->phenotype.volume.total), 						 							 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.solid), 						 							 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.fluid), 						 							 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.fluid_fraction), 	 							 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.nuclear), 					 							 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.nuclear_fluid), 		 							 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.nuclear_solid), 		 							 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.cytoplasmic), 			 							 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.cytoplasmic_fluid), 							 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.cytoplasmic_solid), 							 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.calcified_fraction), 						 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.cytoplasmic_to_nuclear_ratio), 	 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.rupture_volume), 								 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.cytoplasmic_biomass_change_rate), 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.nuclear_biomass_change_rate), 		 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.fluid_change_rate), 							 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.calcification_rate), 						 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.target_solid_cytoplasmic), 			 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.target_solid_nuclear), 					 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.target_fluid_fraction), 					 		 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.target_cytoplasmic_to_nuclear_ratio), 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.relative_rupture_volume), 						 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		/* Next packing data members of class Geometry as this class has an object in class Phenotype */
		/* 4 doubles are to be packed */
		
		len_snd_buf_left = position_left + 4 * sizeof(double);
		snd_buf_left.resize(len_snd_buf_left);
		
		MPI_Pack(&(pCell->phenotype.geometry.radius), 				1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.geometry.nuclear_radius), 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.geometry.surface_area), 	1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.geometry.polarity), 			1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		/* Next packing object of class Mechanics as class Phenotype contains object of this class */
		/* 5 doubles are to be packed */
		
		len_snd_buf_left = position_left + 5 * sizeof(double);
		snd_buf_left.resize(len_snd_buf_left);

		MPI_Pack(&(pCell->phenotype.mechanics.cell_cell_adhesion_strength), 			 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.mechanics.cell_BM_adhesion_strength), 				 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.mechanics.cell_cell_repulsion_strength), 			 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.mechanics.cell_BM_repulsion_strength), 				 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.mechanics.relative_maximum_adhesion_distance), 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);

		/* Next Packing data members of class Motility */
		/* First bool value is to be packed - allocate sizeof(int) */
		
		len_snd_buf_left = position_left + sizeof(int);
		snd_buf_left.resize(len_snd_buf_left);
		if(pCell->phenotype.motility.is_motile == true)
			temp_int = 1;
		else
			temp_int = 0;
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		/* Now pack 2 doubles */
		
		len_snd_buf_left = position_left + 2 * sizeof(double);		
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&(pCell->phenotype.motility.persistence_time), 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.motility.migration_speed),  1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		len_vector = pCell->phenotype.motility.migration_bias_direction.size();
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.motility.migration_bias_direction[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		len_snd_buf_left = position_left + sizeof(double);		
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&(pCell->phenotype.motility.migration_bias), 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		len_snd_buf_left = position_left + sizeof(int);		
		snd_buf_left.resize(len_snd_buf_left);
		if(pCell->phenotype.motility.restrict_to_2D == true)
			temp_int = 1; 
		else
			temp_int = 0; 
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
	
		len_vector = pCell->phenotype.motility.motility_vector.size();
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.motility.motility_vector[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
	
		/* Now packing data members of class Secretion as its object is a data member of class Phenotype */
		
		len_vector = pCell->phenotype.secretion.secretion_rates.size();
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.secretion.secretion_rates[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		len_vector = pCell->phenotype.secretion.uptake_rates.size();
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.secretion.uptake_rates[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		len_vector = pCell->phenotype.secretion.saturation_densities.size();
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.secretion.saturation_densities[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		/* Now packing data members of class Molcular as its object is a data member of class Phenotype */
		
		len_vector = pCell->phenotype.molecular.internalized_total_substrates.size();
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.molecular.internalized_total_substrates[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		len_vector = pCell->phenotype.molecular.fraction_released_at_death.size();
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.molecular.fraction_released_at_death[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		len_vector = pCell->phenotype.molecular.fraction_transferred_when_ingested.size();
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.molecular.fraction_transferred_when_ingested[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		/* Returning to class Cell to pack the remaining data members : 2 bools + 1 std::vector<double> */
		
		len_snd_buf_left = position_left + 2 * sizeof(int);
		snd_buf_left.resize(len_snd_buf_left);
		
		if(pCell->is_out_of_domain == true)
			temp_int = 1;
		else
			temp_int = 0;
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		
		if(pCell->is_movable == true)
			temp_int = 1;
		else
			temp_int = 0;
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		len_vector = pCell->displacement.size();
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->displacement[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		/* Now Packing data members of class Basic_Agent - the parent class of class Cell */
		/* First 2 private data members are to be packed: double and bool 								*/

		temp_double = pCell->get_total_volume();
		len_snd_buf_left = position_left + sizeof(temp_double);
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		if(pCell->get_is_volume_changed() == true)
			temp_int = 1;
		else
			temp_int = 0; 
		
		len_snd_buf_left = position_left + sizeof(temp_int);
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		
		/* Now 3 std::vector<double> are to be packed - but are protected, hence write helper functions */
		
		
		std::vector<double> & temp_cell_source_sink_solver_temp1 = pCell->get_cell_source_sink_solver_temp1();
		len_vector = temp_cell_source_sink_solver_temp1.size(); 
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&temp_cell_source_sink_solver_temp1[0], len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD); 
		
		std::vector<double> & temp_cell_source_sink_solver_temp2 = pCell->get_cell_source_sink_solver_temp2();
		len_vector = temp_cell_source_sink_solver_temp2.size(); 
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&temp_cell_source_sink_solver_temp2[0], len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);

		std::vector<double> & temp_previous_velocity = pCell->get_previous_velocity();
		len_vector = temp_previous_velocity.size(); 
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&temp_previous_velocity[0], len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);

		if(pCell->get_is_active() == true)
			temp_int = 1;
		else
			temp_int = 0;
			
		len_snd_buf_left = position_left + sizeof(int);
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		std::vector<double> & temp_total_extracellular_substrate_change = pCell->get_total_extracellular_substrate_change();
		len_vector = temp_total_extracellular_substrate_change.size(); 
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&temp_total_extracellular_substrate_change[0], len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		/* Now need to pack multiple std::vector<double> *ptr ; type pointers to vectors */
		
		len_vector = pCell->secretion_rates->size(); 
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double); 
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->secretion_rates[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		len_vector = pCell->saturation_densities->size(); 
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double); 
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->saturation_densities[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		len_vector = pCell->uptake_rates->size(); 
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double); 
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->uptake_rates[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);

		len_vector = pCell->internalized_substrates->size(); 
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double); 
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->internalized_substrates[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);

		len_vector = pCell->fraction_released_at_death->size(); 
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double); 
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->fraction_released_at_death[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);		

		len_vector = pCell->fraction_transferred_when_ingested->size(); 
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double); 
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->fraction_transferred_when_ingested[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);	
		
		/* 1 single int type; to be packed next, remember ID and position are already packed */
		
		len_snd_buf_left = position_left + sizeof(int);
		snd_buf_left.resize(len_snd_buf_left);
		MPI_Pack(&(pCell->type), 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		
		len_vector = pCell->velocity.size(); 
		len_snd_buf_left = position_left + sizeof(len_vector) + len_vector * sizeof(double); 
		snd_buf_left.resize(len_snd_buf_left); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->velocity[0]), len_vector, MPI_DOUBLE, &snd_buf_left[0], len_snd_buf_left, &position_left, MPI_COMM_WORLD);
					
	}
	
		/*----------------------------------------------------------------------*/
		/* PACKING OF BUFFER TO BE SENT TO LEFT PROCESS ENDS 										*/
		/* THE RIGHT BUFFER IS IDENTICAL EXCEPT FOR _LEFT IS REPLACED BY _RIGHT */
		/*----------------------------------------------------------------------*/
		 
		if(pCell->crossed_to_right_subdomain == true)
		{
		
			//pCell->crossed_to_right_subdomain = false;
			/* Cell ID first - needed to create a cell */
						
			len_snd_buf_right = position_right + sizeof(pCell->ID);
			snd_buf_right.resize(len_snd_buf_right);
			MPI_Pack(&(pCell->ID), 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
			
			/* Pack double position[0], position[1], position[2] - needed to assign position 							 */
			/* This is the only array where I am not packing the length of the array with the array values */
			/* as we're only dealing with 3-dimensions here, hence position has x, y, z values						 */
			
			len_snd_buf_right = position_right + 3 * sizeof(pCell->position[0]); 
			snd_buf_right.resize(len_snd_buf_right); 
			MPI_Pack(&(pCell->position[0]), 3, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			
			/* Pack length of string (MPI_INT) and then string pCell->type_name (MPI_CHAR) */
			
			len_str = pCell->type_name.length(); 
			len_snd_buf_right = position_right + sizeof(len_str) + len_str;
			snd_buf_right.resize(len_snd_buf_right); 
			MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
			MPI_Pack(&(pCell->type_name[0]), len_str , MPI_CHAR, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
			
			/* Packing unordered_map<std::string, int> in Custom_Data class */
			/* First pack the number of entries in the unordered_map 				*/
			
			len_int = pCell->custom_data.get_name_to_index_map().size();
			len_snd_buf_right = position_right + sizeof(len_int);
			snd_buf_right.resize(len_snd_buf_right);  
			MPI_Pack(&len_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);  
			
			/* Now pack the <string, int> pairs but store length(string) before string */
					
			for(auto it = (pCell->custom_data.get_name_to_index_map()).cbegin(); it != (pCell->custom_data.get_name_to_index_map()).cend(); it++)
			{
				/* resize buffer to contain length of string, string and integer */
				
				temp_str = it->first; 
				temp_int = it->second; 
				len_str = temp_str.length();
				len_snd_buf_right = position_right + sizeof(len_str) + len_str + sizeof(temp_int); 
				snd_buf_right.resize(len_snd_buf_right);
				 
				MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
				MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
				MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD ); 			
			}	
			
		 /* Packing members of class 'Variable' (as we are going depth first i.e. Inorder traversal of data members) */
		 
		 	for(int i=0; i<pCell->custom_data.variables.size(); i++)
		 	{
		    temp_str = pCell->custom_data.variables[i].name; 
		    len_str = temp_str.length();
		    len_snd_buf_right = position_right + len_str + sizeof(len_str); 
		    snd_buf_right.resize(len_snd_buf_right);
		    MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		    MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		    
		    /* There was a mistake in MPI_Pack here, instead of &snd_buf_right[0], I had &snd_buf_right */
		    temp_double = pCell->custom_data.variables[i].value;
		    len_snd_buf_right = position_right + sizeof(temp_double); 
		    MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		    
		    temp_str = pCell->custom_data.variables[i].units; 
		    len_str = temp_str.length();
		    len_snd_buf_right = position_right + len_str + sizeof(len_str); 
		    snd_buf_right.resize(len_snd_buf_right);
		    MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		    MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);		     
			}
			
			/* Packing members of class 'Vector_Variable' i.e. string name, a vector<double> value and string units */
			
			for(int i=0; i<pCell->custom_data.vector_variables.size(); i++)
			{
				temp_str = pCell->custom_data.vector_variables[i].name; 
		    len_str = temp_str.length();
		    len_snd_buf_right = position_right + len_str + sizeof(len_str); 
		    snd_buf_right.resize(len_snd_buf_right);
		    MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		    MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		    
		    len_vector = pCell->custom_data.vector_variables[i].value.size(); 
		    len_snd_buf_right = position_right + len_vector * sizeof(double);
		    snd_buf_right.resize(len_snd_buf_right); 
		    MPI_Pack(&(pCell->custom_data.vector_variables[i].value[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		    
		    temp_str = pCell->custom_data.vector_variables[i].units; 
		    len_str = temp_str.length();
		    len_snd_buf_right = position_right + len_str + sizeof(len_str); 
		    snd_buf_right.resize(len_snd_buf_right);
		    MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		    MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);  	    
			}
			
			/* Packing members of class 'Cell_Parameters' except for pReference_live_phenotype pointer */
			/* There are 9 doubles and 1 int in this class except for the pointer */
			
			/*=====================================================================*/
			/* pReference_live_phenotype * IS NOT PACKED 													 */
			/*=====================================================================*/
			
			len_snd_buf_right = position_right + 9 * sizeof(double) + 1 * sizeof(int); 
			snd_buf_right.resize(len_snd_buf_right); 
			
			MPI_Pack(&(pCell->parameters.o2_hypoxic_threshold), 				1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.o2_hypoxic_response), 					1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.o2_hypoxic_saturation), 				1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.o2_proliferation_saturation), 	1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.o2_proliferation_threshold), 	1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.o2_reference), 								1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.o2_necrosis_threshold), 				1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.o2_necrosis_max), 							1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.max_necrosis_rate), 						1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			MPI_Pack(&(pCell->parameters.necrosis_type), 								1, MPI_INT, 	 &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		
		
		
		/* Packing data members of Cell_Functions in depth first manner */
		/* It only has Cycle_Model data member, so pack Cycle_Model data members depth first */
		
		std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_model = pCell->functions.cycle_model.get_inverse_index_maps();
		len_vector = inverse_index_maps_cycle_model.size(); 
		len_snd_buf_right = position_right + sizeof(len_vector); 
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		
		for(int i=0; i<len_vector; i++)
		{
			 len_int = inverse_index_maps_cycle_model[i].size();
			 len_snd_buf_right = position_right + sizeof(len_int) + len_int * 2 * sizeof(int); 
			 snd_buf_right.resize(len_snd_buf_right);
			 MPI_Pack(&len_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			 
			for(auto it = inverse_index_maps_cycle_model[i].cbegin(); it != inverse_index_maps_cycle_model[i].cend(); it++)
			{	 
				MPI_Pack(&(it->first),  1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
				MPI_Pack(&(it->second), 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD ); 			
			}				 
		}
		
		temp_str = pCell->functions.cycle_model.name; 
		len_str  = temp_str.length();
		len_snd_buf_right = position_right + sizeof(len_str) + len_str; 
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		
		len_snd_buf_right = position_right + sizeof(int);
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&(pCell->functions.cycle_model.code), 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		
		/* Packing of data members of Phases as its object is a member of Cycle_Model */
		/* std::vector<Phase> phases; is a member of it - its a vector */		
		
		len_vector = pCell->functions.cycle_model.phases.size(); 
		len_snd_buf_right = position_right + sizeof(len_vector);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		
		for(int i=0; i<len_vector; i++)
		{
			len_snd_buf_right = position_right + 2 * sizeof(int);
			snd_buf_right.resize(len_snd_buf_right); 
			 
			temp_int = pCell->functions.cycle_model.phases[i].index; 
			MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
			temp_int = pCell->functions.cycle_model.phases[i].code;
			MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
			
			temp_str = pCell->functions.cycle_model.phases[i].name; 
			len_str  = temp_str.length();
			len_snd_buf_right = position_right + sizeof(len_str) + len_str; 
			snd_buf_right.resize(len_snd_buf_right); 
			MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			
			/* There are 2 remaining bool values which I will pack in 2 ints */
			len_snd_buf_right = position_right + 2 * sizeof(int);
			snd_buf_right.resize(len_snd_buf_right);
			
			if(pCell->functions.cycle_model.phases[i].division_at_phase_exit == true)
				temp_int = 1; 
			else
				temp_int = 0;
			MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
			
			if(pCell->functions.cycle_model.phases[i].removal_at_phase_exit == true)
				temp_int = 1; 
			else
				temp_int = 0;
			MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);		
		}
		
		/* Packing members of Phase_Link - std::vector<std::vector<Phase_Link>> phase_links is a member of */
		/* class Cycle_Model which is a member of Functions - we're packing in DFS i.e. Inorder traversal	 */
		
		len_vector = pCell->functions.cycle_model.phase_links.size();
		len_snd_buf_right = position_right + sizeof(len_vector);
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		for(int i=0; i<len_vector; i++)
		{
			len_int = pCell->functions.cycle_model.phase_links[i].size();
			len_snd_buf_right = position_right + sizeof(len_int) + len_int * 3 * sizeof(int); 
			snd_buf_right.resize(len_snd_buf_right); 
			MPI_Pack(&len_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			for(int j=0; j<len_int; j++)
			{
				temp_int = pCell->functions.cycle_model.phase_links[i][j].start_phase_index; 
				MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
				temp_int = pCell->functions.cycle_model.phase_links[i][j].end_phase_index;
				MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
				if(pCell->functions.cycle_model.phase_links[i][j].fixed_duration == true)
					temp_int = 1;
				else
					temp_int = 0; 
				MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
			}
		}
		
		/* Going back to packing the remaning data fields of Cycle_Model */
		/* Next field is is 'default_phase_index' of Cycle_Model class */
		
		temp_int = pCell->functions.cycle_model.default_phase_index; 
		len_snd_buf_right = position_right + sizeof(temp_int);
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		
		 /* Packing data fields of Cycle_Data as it is an object of Cycle_Model class above */
		 /* First field is std::vector<std::unordered_map<int,int>> inverse_index_maps;			*/
		 
		std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_data = pCell->functions.cycle_model.data.get_inverse_index_maps();	 
		len_vector = inverse_index_maps_cycle_data.size(); 
		len_snd_buf_right = position_right + sizeof(len_vector); 
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		
		for(int i=0; i<len_vector; i++)
		{
			 len_int = inverse_index_maps_cycle_data[i].size();
			 len_snd_buf_right = position_right + sizeof(len_int) + len_int * 2 * sizeof(int); 
			 snd_buf_right.resize(len_snd_buf_right);
			 MPI_Pack(&len_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			 
			for(auto it = inverse_index_maps_cycle_data[i].cbegin(); it != inverse_index_maps_cycle_data[i].cend(); it++)
			{	 
				MPI_Pack(&(it->first),  1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
				MPI_Pack(&(it->second), 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD ); 			
			}				 
		}
		
			/*=====================================================================*/
			/* Cycle_Model* pCycle_Model;  IS NOT PACKED 													 */
			/*=====================================================================*/
			
			temp_str = pCell->functions.cycle_model.data.time_units;
			len_str = temp_str.length();
			len_snd_buf_right = position_right + sizeof(len_str) + len_str; 
			snd_buf_right.resize(len_snd_buf_right);
			MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		 
		 len_vector = pCell->functions.cycle_model.data.transition_rates.size(); 
		 len_snd_buf_right = position_right + sizeof(len_vector);
		 snd_buf_right.resize(len_snd_buf_right);
		 MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		 
		 for(int i=0; i<len_vector; i++)
		 {
		 	len_int = pCell->functions.cycle_model.data.transition_rates[i].size();
		 	len_snd_buf_right = position_right + sizeof(len_int) + len_int * sizeof(double);
		 	snd_buf_right.resize(len_snd_buf_right);
		 	MPI_Pack(&len_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		 	
		 	for(int j=0; j<len_int; j++)
		 	{
		 		temp_double = pCell->functions.cycle_model.data.transition_rates[i][j];
		 		MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		 	}
		 }
		 
		 /* 1 int and 1 double are next to be packed in class Cycle_Data */
		 
		 len_snd_buf_right = position_right + sizeof(int) + sizeof(double);
		 snd_buf_right.resize(len_snd_buf_right);
		 MPI_Pack(&(pCell->functions.cycle_model.data.current_phase_index), 	1, MPI_INT, 	 &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		 MPI_Pack(&(pCell->functions.cycle_model.data.elapsed_time_in_phase), 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 

		/* Now packing members of class Cell_State as this is the next member in class Cell */
		/* has std::vector<double> orientation and double simple_pressure 									*/
		
		len_vector = pCell->state.orientation.size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double) + sizeof(double); 
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		MPI_Pack(&(pCell->state.orientation[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);		
		MPI_Pack(&(pCell->state.simple_pressure), 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		
		/* Packing data members of class Phenotype next */
		
		/* First it has 2 boolean variables - allocate sizeof(int) space for both of them */
		len_snd_buf_right = position_right + 2 * sizeof(int);
		snd_buf_right.resize(len_snd_buf_right);
		if(pCell->phenotype.flagged_for_division == true)
			temp_int = 1;
		else
			temp_int = 0;
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		if(pCell->phenotype.flagged_for_removal == true)
			temp_int = 1;
		else
			temp_int = 0;
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		/* Packing of members of class Cycle as the class Phenotype contains a data member of this class */
		/* Not packing pCell->phenotype.cycle.pCycle_Model pointer data member *
		
			/*=====================================================================*/
			/* Cycle_Model* pCycle_Model;  IS NOT PACKED 													 */
			/*=====================================================================*/
		
		 /* Packing data fields of Cycle_Data as it is an object of Cycle class which is    */
		 /* a data member of the Phenotype class which is a member of class Cell						*/
		 /* First field is std::vector<std::unordered_map<int,int>> inverse_index_maps;			*/
		 
		std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_data_1 = pCell->phenotype.cycle.data.get_inverse_index_maps();	 
		len_vector = inverse_index_maps_cycle_data_1.size(); 
		len_snd_buf_right = position_right + sizeof(len_vector);
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		
		for(int i=0; i<len_vector; i++)
		{
			 len_int = inverse_index_maps_cycle_data_1[i].size();
			 len_snd_buf_right = position_right + sizeof(len_int) + len_int * 2 * sizeof(int); 
			 snd_buf_right.resize(len_snd_buf_right);
			 MPI_Pack(&len_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			 
			for(auto it = inverse_index_maps_cycle_data_1[i].cbegin(); it != inverse_index_maps_cycle_data_1[i].cend(); it++)
			{	 
				MPI_Pack(&(it->first),  1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
				MPI_Pack(&(it->second), 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD ); 			
			}				 
		}
		
			/*=====================================================================*/
			/* Cycle_Model* pCycle_Model;  IS NOT PACKED  												 */
			/*=====================================================================*/
			
			temp_str = pCell->phenotype.cycle.data.time_units;
			len_str = temp_str.length();
			len_snd_buf_right = position_right + sizeof(len_str) + len_str; 
			snd_buf_right.resize(len_snd_buf_right);
			MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		 
		 len_vector = pCell->phenotype.cycle.data.transition_rates.size(); 
		 len_snd_buf_right = position_right + sizeof(len_vector);
		 snd_buf_right.resize(len_snd_buf_right);
		 MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		 
		 for(int i=0; i<len_vector; i++)
		 {
		 	len_int = pCell->phenotype.cycle.data.transition_rates[i].size();
		 	len_snd_buf_right = position_right + sizeof(len_int) + len_int * sizeof(double);
		 	snd_buf_right.resize(len_snd_buf_right);
		 	MPI_Pack(&len_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		 	
		 	for(int j=0; j<len_int; j++)
		 	{
		 		temp_double = pCell->phenotype.cycle.data.transition_rates[i][j];
		 		MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		 	}
		 }
		 
		 /* 1 int and 1 double are next to be packed in class Cycle_Data */
		 
		 len_snd_buf_right = position_right + sizeof(int) + sizeof(double);
		 snd_buf_right.resize(len_snd_buf_right);
		 MPI_Pack(&(pCell->phenotype.cycle.data.current_phase_index), 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		 MPI_Pack(&(pCell->phenotype.cycle.data.elapsed_time_in_phase), 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		/* Now packing data members of class Death as there is an object of class Death in Phenotype class */
		
		len_vector = pCell->phenotype.death.rates.size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double); 
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.death.rates[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		
		/* Next field is std::vector<Cycle_Model*> models ;  this is a vector of Cycle_Model * pointers */
		
		/*=====================================================================*/
		/* std::vector<Cycle_Model*> models ;  IS NOT PACKED  								 */
		/*=====================================================================*/
		
		/* Now packing members of class Death_Parameters as class Death contains */
		/* std::vector<Death_Parameters> parameters 														 */
		
		len_vector = pCell->phenotype.death.parameters.size();
		len_snd_buf_right = position_right + sizeof(len_vector);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		 
		for(int i=0; i<len_vector; i++)
		{
			temp_str = pCell->phenotype.death.parameters[i].time_units;
			len_str  = temp_str.length();
			len_snd_buf_right = position_right + sizeof(len_str) + len_str; 
			snd_buf_right.resize(len_snd_buf_right); 
			MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
			MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
			
			/* 6 double data types are to be packed now */
			len_snd_buf_right = position_right + 6 * sizeof(double);
			snd_buf_right.resize(len_snd_buf_right);
			
			temp_double = pCell->phenotype.death.parameters[i].unlysed_fluid_change_rate;
			MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
			temp_double = pCell->phenotype.death.parameters[i].lysed_fluid_change_rate;
			MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);			
			temp_double = pCell->phenotype.death.parameters[i].cytoplasmic_biomass_change_rate;
			MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
			temp_double = pCell->phenotype.death.parameters[i].nuclear_biomass_change_rate;
			MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
			temp_double = pCell->phenotype.death.parameters[i].calcification_rate;
			MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
			temp_double = pCell->phenotype.death.parameters[i].relative_rupture_volume;
			MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		}
		 
		 /* Next 1 bool value (allocate sizeof(int)) and 1 int value are to be packed */
		 
		 len_snd_buf_right = position_right + 2 * sizeof(int);
		 snd_buf_right.resize(len_snd_buf_right);
		 
		 if(pCell->phenotype.death.dead == true)
		 		temp_int = 1;
		 else
		 		temp_int = 0; 
		 MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		 
		 temp_int = pCell->phenotype.death.current_death_model_index;
		 MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		/* Next data members of class Volume are packed as class Phenotype contains its object */
		/* 22 doubles are to be packed */
		
		len_snd_buf_right = position_right + 22 * sizeof(double); 
		snd_buf_right.resize(len_snd_buf_right);
		
		MPI_Pack(&(pCell->phenotype.volume.total), 						 							 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.solid), 						 							 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.fluid), 						 							 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.fluid_fraction), 	 							 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.nuclear), 					 							 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.nuclear_fluid), 		 							 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.nuclear_solid), 		 							 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.cytoplasmic), 			 							 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.cytoplasmic_fluid), 							 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.cytoplasmic_solid), 							 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.calcified_fraction), 						 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.cytoplasmic_to_nuclear_ratio), 	 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.rupture_volume), 								 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.cytoplasmic_biomass_change_rate), 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.nuclear_biomass_change_rate), 		 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.fluid_change_rate), 							 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.calcification_rate), 						 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.target_solid_cytoplasmic), 			 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.target_solid_nuclear), 					 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.target_fluid_fraction), 					 		 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.target_cytoplasmic_to_nuclear_ratio), 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.volume.relative_rupture_volume), 						 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		/* Next packing data members of class Geometry as this class has an object in class Phenotype */
		/* 4 doubles are to be packed */
		
		len_snd_buf_right = position_right + 4 * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		
		MPI_Pack(&(pCell->phenotype.geometry.radius), 				1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.geometry.nuclear_radius), 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.geometry.surface_area), 	1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.geometry.polarity), 			1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		/* Next packing object of class Mechanics as class Phenotype contains object of this class */
		/* 5 doubles are to be packed */
		
		len_snd_buf_right = position_right + 5 * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);

		MPI_Pack(&(pCell->phenotype.mechanics.cell_cell_adhesion_strength), 			 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.mechanics.cell_BM_adhesion_strength), 				 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.mechanics.cell_cell_repulsion_strength), 			 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.mechanics.cell_BM_repulsion_strength), 				 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.mechanics.relative_maximum_adhesion_distance), 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		/* Next Packing data members of class Motility */
		/* First bool value is to be packed - allocate sizeof(int) */
		
		len_snd_buf_right = position_right + sizeof(int);
		snd_buf_right.resize(len_snd_buf_right);
		if(pCell->phenotype.motility.is_motile == true)
			temp_int = 1;
		else
			temp_int = 0;
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		/* Now pack 2 doubles */
		
		len_snd_buf_right = position_right + 2 * sizeof(double);		
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&(pCell->phenotype.motility.persistence_time), 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.motility.migration_speed),  1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		len_vector = pCell->phenotype.motility.migration_bias_direction.size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.motility.migration_bias_direction[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		len_snd_buf_right = position_right + sizeof(double);		
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&(pCell->phenotype.motility.migration_bias), 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		len_snd_buf_right = position_right + sizeof(int);		
		snd_buf_right.resize(len_snd_buf_right);
		if(pCell->phenotype.motility.restrict_to_2D == true)
			temp_int = 1; 
		else
			temp_int = 0; 
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
	
		len_vector = pCell->phenotype.motility.motility_vector.size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.motility.motility_vector[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
	
		/* Now packing data members of class Secretion as its object is a data member of class Phenotype */
		
		len_vector = pCell->phenotype.secretion.secretion_rates.size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.secretion.secretion_rates[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		len_vector = pCell->phenotype.secretion.uptake_rates.size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.secretion.uptake_rates[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		len_vector = pCell->phenotype.secretion.saturation_densities.size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.secretion.saturation_densities[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		/* Now packing data members of class Molcular as its object is a data member of class Phenotype */
		
		len_vector = pCell->phenotype.molecular.internalized_total_substrates.size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.molecular.internalized_total_substrates[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		len_vector = pCell->phenotype.molecular.fraction_released_at_death.size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.molecular.fraction_released_at_death[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		len_vector = pCell->phenotype.molecular.fraction_transferred_when_ingested.size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.molecular.fraction_transferred_when_ingested[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		/* Returning to class Cell to pack the remaining data members : 2 bools + 1 std::vector<double> */
		
		len_snd_buf_right = position_right + 2 * sizeof(int);
		snd_buf_right.resize(len_snd_buf_right);
		
		if(pCell->is_out_of_domain == true)
			temp_int = 1;
		else
			temp_int = 0;
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		
		if(pCell->is_movable == true)
			temp_int = 1;
		else
			temp_int = 0;
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		len_vector = pCell->displacement.size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->displacement[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		/* Now Packing data members of class Basic_Agent - the parent class of class Cell */
		/* First 2 private data members are to be packed: double and bool 								*/

		temp_double = pCell->get_total_volume();
		len_snd_buf_right = position_right + sizeof(temp_double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		if(pCell->get_is_volume_changed() == true)
			temp_int = 1;
		else
			temp_int = 0; 
		
		len_snd_buf_right = position_right + sizeof(temp_int);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		
		/* Now 3 std::vector<double> are to be packed - but are protected, hence write helper functions */
		
		
		std::vector<double> & temp_cell_source_sink_solver_temp1 = pCell->get_cell_source_sink_solver_temp1();
		len_vector = temp_cell_source_sink_solver_temp1.size(); 
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&temp_cell_source_sink_solver_temp1[0], len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		
		std::vector<double> & temp_cell_source_sink_solver_temp2 = pCell->get_cell_source_sink_solver_temp2();
		len_vector = temp_cell_source_sink_solver_temp2.size(); 
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&temp_cell_source_sink_solver_temp2[0], len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		std::vector<double> & temp_previous_velocity = pCell->get_previous_velocity();
		len_vector = temp_previous_velocity.size(); 
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&temp_previous_velocity[0], len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		if(pCell->get_is_active() == true)
			temp_int = 1;
		else
			temp_int = 0;
			
		len_snd_buf_right = position_right + sizeof(int);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&temp_int, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		std::vector<double> & temp_total_extracellular_substrate_change = pCell->get_total_extracellular_substrate_change();
		len_vector = temp_total_extracellular_substrate_change.size(); 
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&temp_total_extracellular_substrate_change[0], len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		/* Now need to pack multiple std::vector<double> *ptr ; type pointers to vectors */
		
		len_vector = pCell->secretion_rates->size(); 
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double); 
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->secretion_rates[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		len_vector = pCell->saturation_densities->size(); 
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double); 
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->saturation_densities[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		len_vector = pCell->uptake_rates->size(); 
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double); 
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->uptake_rates[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		len_vector = pCell->internalized_substrates->size(); 
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double); 
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->internalized_substrates[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		len_vector = pCell->fraction_released_at_death->size(); 
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double); 
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->fraction_released_at_death[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);		

		len_vector = pCell->fraction_transferred_when_ingested->size(); 
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double); 
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->fraction_transferred_when_ingested[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);	
		
		/* 1 single int type; to be packed next, remember ID and position are already packed */
		
		len_snd_buf_right = position_right + sizeof(int);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&(pCell->type), 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		
		len_vector = pCell->velocity.size(); 
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double); 
		snd_buf_right.resize(len_snd_buf_right); 
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->velocity[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
			
		}
	}
}

void Cell_Container::unpack(mpi_Environment &world)
{
		
		int position_left = 0, position_right = 0;
		int size;
		int cell_ID;  
		
		if(no_of_cells_from_right > 0)
		{
			size = rcv_buf_right.size();
			MPI_Unpack(&rcv_buf_right[0], size, &position_right, &cell_ID, 1, MPI_INT, MPI_COMM_WORLD);
			std::cout<<"Rank = "<<world.rank<< " First Cell ID received = "<<cell_ID<<std::endl;  
		}
		
		if(no_of_cells_from_left > 0)
		{
			size = rcv_buf_left.size();
			MPI_Unpack(&rcv_buf_left[0], size, &position_left, &cell_ID, 1, MPI_INT, MPI_COMM_WORLD);
			std::cout<<"Rank = "<<world.rank<< " First Cell ID received = "<<cell_ID<<std::endl;  
		}
		
		
}

};

