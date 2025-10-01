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
#include "../../../../modules/PhysiCell_settings.h"
#include "../../../../DistPhy/DistPhy_Utils.h"
#include "../../../../DistPhy/DistPhy_Collective.h"

using namespace DistPhy::mpi;

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

void create_cell_types( void )
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 

	Cell_Definition* pCD = find_cell_definition( "cancer cell"); 
	pCD->functions.update_phenotype = tumor_cell_phenotype_with_oncoprotein; 

	pCD->parameters.o2_proliferation_saturation = 38; 
	pCD->parameters.o2_reference = 38; 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	//display_cell_definitions( std::cout ); 
	
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

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
  double y_spacing= cell_radius*2;
  double z_spacing= cell_radius*sqrt(3);
	
	//Attempt to generate very small number of cells
	 //double x_spacing= 10;
	 //double y_spacing= 10;
	 //double z_spacing= 10;

	std::vector<double> tempPoint(3,0.0);
	// std::vector<double> cylinder_center(3,0.0);

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

/*------------------------------------------------------------------------*/
/* Parallel version of setup_tissue(), replacing this function completely */
/* by Miguel's version of setup_tissue and then parallelizing             */
/*------------------------------------------------------------------------*/

void setup_tissue(Microenvironment &m, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
  
	Cell* pC;

	// Place a cluster of tumor cells at the center
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double tumor_radius = parameters.doubles( "tumor_radius" ); // 250.0; now changed to 150 in PhysiCell_settings.xml file	
	int i = parameters.doubles.find_index( "tumor_radius" ); 
	
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
        generated_positions_at_root = create_cell_sphere_positions(cell_radius,tumor_radius);   //Generate the cell positions
        
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
	double x_outer = tumor_radius; 
	double y = 0.0; 

	Cell_Definition* pCD = find_cell_definition( "cancer cell"); 
	
	
	for( int i=0; i < mc.my_no_of_cell_IDs; i++ )
	{
		
		//pCell = create_cell(mc.my_cell_IDs[i]); // tumor cell --> This has to be replaced by create_cell(mc.my_cell_IDs[i])
		pCell = create_cell((*pCD), mc.my_cell_IDs[i]);	
		pCell->assign_position(mc.my_cell_coords[3*i],mc.my_cell_coords[3*i+1],mc.my_cell_coords[3*i+2],world, cart_topo); //pCell->assign_position( positions[i] );
				
	} 
	
	return; 
}

// Custom cell phenotype function to scale immunostimulatory factor with hypoxia 

void tumor_cell_phenotype_with_oncoprotein( Cell* pCell, Phenotype& phenotype, double dt )
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



/* https://www.karger.com/Article/Fulltext/494069 */ 

void macrophage_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	// find my cell definition 
	static Cell_Definition* pCD = find_cell_definition( pCell->type_name ); 

	// sample environment 

	static int nPIF = microenvironment.find_density_index( "pro-inflammatory" ); 
	static int nDebris = microenvironment.find_density_index( "debris"); 
	static int nQ = microenvironment.find_density_index("quorum");

	// if dead, release debris
	if( phenotype.death.dead == true )
	{
		phenotype.secretion.net_export_rates[nDebris] = phenotype.volume.total; 
		pCell->functions.update_phenotype = NULL; 
		return;
	}	

	double * samples = pCell->nearest_density_vector(); 
	double PIF = samples[nPIF];
	double debris = samples[nDebris]; 
	double Q = samples[nQ];

	// sample contacts 

	static int bacteria_type = find_cell_definition( "bacteria")->type; 

	int num_bacteria = 0; 
	int num_dead = 0; 
	for( int n=0; n < pCell->state.neighbors.size(); n++ )
	{
		Cell* pC = pCell->state.neighbors[n]; 
		if( pC->phenotype.death.dead == true )
		{ num_dead++; }
		else
		{ 
			if( pC->type == bacteria_type )
			{ num_bacteria++; }
		}
	}

	// contact with dead cells or bacteria, or debris 
	// increases secretion of pro-inflammatory 

	static double secretion_dead_sensitivity = 1; 
	static double secretion_bacteria_sensitivity = 1; 
	static double secretion_debris_sensitivity = 2; 
	static double secretion_quorum_sensitivity = 5; 

	double base_val = pCD->phenotype.secretion.secretion_rates[nPIF]; 
	double max_response = 10; // phenotype.volume.total; 
	double signal = 
		secretion_dead_sensitivity*num_dead + 
		secretion_bacteria_sensitivity*num_bacteria + 
		secretion_debris_sensitivity*debris + 
		secretion_quorum_sensitivity*Q; 
	double half_max = pCD->custom_data["secretion_halfmax"]; // 0.5; // 0.5; 
	double hill = Hill_response_function( signal , half_max , 1.5 ); 
	
	
	phenotype.secretion.secretion_rates[nPIF] = base_val + (max_response-base_val)*hill; 
	
/*	
	#pragma omp critical
	{
	std::cout << "secretion index: " << nPIF << " base: " << base_val << " max: " << max_response << " actual: " << phenotype.secretion.secretion_rates[nPIF] << std::endl; 
	std::cout << "\tsignal: " << signal << " vs halfmax: " << half_max << std::endl; 
	std::cout << "\t\tdead: " << num_dead << " bac: " << num_bacteria << " debris: " << debris << " Q: " << Q << std::endl; 
	std::cout << "\t\t\tsaturation: " << phenotype.secretion.saturation_densities[nPIF]<< std::endl; 
	}
*/	

	// chemotaxis bias increases with debris or quorum factor 

	static double bias_debris_sensitivity = 0.1; 
	static double bias_quorum_sensitivity = 1; 

	base_val = pCD->phenotype.motility.migration_bias; 
	max_response = 0.75; 
	signal = bias_debris_sensitivity*debris + 
		bias_quorum_sensitivity*Q ; // + 10 * PIF; 
	half_max = pCD->custom_data["migration_bias_halfmax"]; // 0.01 // 0.005 //0.1 // 0.05
	hill = Hill_response_function( signal , half_max , 1.5 ); 
	phenotype.motility.migration_bias = base_val + (max_response-base_val)*hill; 	

/*
	#pragma omp critical 
	{
	std::cout << "signal: " << signal << " halfmax: " << half_max 
	<< " hill: " << hill << std::endl; 
	
	std::cout << "\tbase: " << base_val 
	<< " max: " << max_response 
	<< " actual: " << phenotype.motility.migration_bias << std::endl; 
	}
*/

	// migration speed slows down in the presence of debris or quorum factor 

	base_val = pCD->phenotype.motility.migration_speed; 
	max_response = 0.1 * base_val; 
	signal = bias_debris_sensitivity*debris + 
		bias_quorum_sensitivity*Q ; // + 10 * PIF; 
	half_max = pCD->custom_data["migration_speed_halfmax"]; // 0.1 // 0.05 
	hill = Hill_response_function( signal , half_max , 1.5 ); 
	phenotype.motility.migration_speed = base_val + (max_response-base_val)*hill; 	

	return; 
}

void CD8Tcell_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	// find my cell definition 
	static Cell_Definition* pCD = find_cell_definition( pCell->type_name ); 

	// sample environment 

	static int nR = microenvironment.find_density_index( "resource");
	static int nTox = microenvironment.find_density_index( "toxin");
	static int nDebris = microenvironment.find_density_index( "debris" );
	static int nPIF = microenvironment.find_density_index( "pro-inflammatory"); 
	
	double * samples = pCell->nearest_density_vector(); 
	double PIF = samples[nPIF];	
	
	// if dead, release debris
	if( phenotype.death.dead == true )
	{
		phenotype.secretion.net_export_rates[nDebris] = phenotype.volume.total; 
		pCell->functions.update_phenotype = NULL; 
		return;
	}
	
	// migration bias increases with pro-inflammatory 

	double signal = PIF; 
	double base_val = pCD->phenotype.motility.migration_bias; 
	double max_val = 0.75; 
	double half_max = pCD->custom_data["migration_bias_halfmax"]; // 0.05 // 0.25 
	double hill = Hill_response_function( PIF , half_max , 1.5 ); 

	phenotype.motility.migration_bias = base_val + (max_val-base_val)*hill; 
	
/*	
	#pragma omp critical 
	{
		std::cout << "signal: " << signal << " halfmax: " << half_max 
		<< " hill: " << hill << std::endl; 
		
		std::cout << "\tbase: " << base_val 
		<< " max: " << max_val 
		<< " actual: " << phenotype.motility.migration_bias << std::endl; 
	}	
*/	

	return; 
}

void neutrophil_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	// find my cell definition 
	static Cell_Definition* pCD = find_cell_definition( pCell->type_name ); 

	// sample environment 

	static int nR = microenvironment.find_density_index( "resource");
	static int nTox = microenvironment.find_density_index( "toxin");
	static int nDebris = microenvironment.find_density_index( "debris" );
	static int nPIF = microenvironment.find_density_index( "pro-inflammatory"); 
	
	double * samples = pCell->nearest_density_vector(); 
	double PIF = samples[nPIF];	
	
	// if dead, release debris
	if( phenotype.death.dead == true )
	{
		phenotype.secretion.net_export_rates[nDebris] = phenotype.volume.total; 
		pCell->functions.update_phenotype = NULL; 
		return;
	}

	// migration bias increases with pro-inflammatory 

	double signal = PIF; 
	double base_val = pCD->phenotype.motility.migration_bias; 
	double max_val = 0.75; 
	double half_max = pCD->custom_data["migration_bias_halfmax"]; // 0.25 
	double hill = Hill_response_function( PIF , half_max , 1.5 ); 

	phenotype.motility.migration_bias = base_val + (max_val-base_val)*hill; 

	return; 
}

void stem_cell_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	// find my cell definition 
	static Cell_Definition* pCD = find_cell_definition( pCell->type_name ); 

	// sample environment 

	static int nR = microenvironment.find_density_index( "resource");
	static int nTox = microenvironment.find_density_index( "toxin");
	static int nDebris = microenvironment.find_density_index( "debris" ); 

	// if dead, release debris
	if( phenotype.death.dead == true )
	{
		phenotype.secretion.net_export_rates[nDebris] = phenotype.volume.total; 
		pCell->functions.update_phenotype = NULL; 
		return;
	}

	double * samples = pCell->nearest_density_vector(); 
	double R = samples[nR];
	double toxin = samples[nTox];

	// sample contacts 

	static int stem_type = find_cell_definition( "stem")->type; 
	static int diff_type = find_cell_definition( "differentiated")->type; 
	static int bacteria_type = find_cell_definition( "bacteria")->type; 

	int num_stem = 0; 
	int num_differentiated = 0; 
	int num_bacteria = 0; 
	int num_dead = 0; 
	for( int n=0; n < pCell->state.neighbors.size(); n++ )
	{
		Cell* pC = pCell->state.neighbors[n]; 
		if( pC->phenotype.death.dead == true )
		{ num_dead++; }
		else
		{ 
			if( pC->type == stem_type )
			{ num_stem++; }
			if( pC->type == num_differentiated )
			{ num_differentiated++; }
			if( pC->type == bacteria_type )
			{ num_bacteria++; }
		}
	}

	// contact with a stem cell increases differentiation 
	static double max_stem_diff = parameters.doubles("max_stem_differentiation"); // 0.0075 
	static double stem_diff_halfmax = pCD->custom_data["differentiation_contact_halfmax"]; // 0.1 

	double base_val = 0; // phenotype.cell_transformations.transformation_rates[diff_type]; 
	double max_val = max_stem_diff; // 0.0075;
	double signal = num_stem; 
	double half_max = stem_diff_halfmax; // 0.1; 
	double hill = Hill_response_function( signal, half_max , 1.5 ); 
	phenotype.cell_transformations.transformation_rates[diff_type] = base_val + (max_val-base_val)*hill; 

	// contact with a differentiated cell reduces proliferation 
	// high rate of proliferation unless in contact with a differentiated cell 

	static double stem_cycling_halfmax = pCD->custom_data["cycling_contact_halfmax"]; // 0.1; 

	base_val = pCD->phenotype.cycle.data.exit_rate(0); // 0.002; 
	max_val = 0.0; 
	signal = num_differentiated; 
	half_max = stem_cycling_halfmax; //  0.1; 
	hill = Hill_response_function( signal, half_max , 1.5 ); 
	phenotype.cycle.data.exit_rate(0) = base_val + (max_val-base_val)*hill; 

	// resource reduces necrotic death 

	max_val = 0.0028;  
	static int nNecrosis = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model );
	static double stem_saturation_necrosis = pCD->custom_data["necrosis_saturation_resource"];
	static double stem_threshold_necrosis = pCD->custom_data["necrosis_threshold_resource"];

	phenotype.death.rates[nNecrosis] = max_val * 
		decreasing_linear_response_function( R, stem_saturation_necrosis, stem_threshold_necrosis );

	// toxin increases apoptotic death 
	
	static int nApoptosis = phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );

	static double toxicity_halfmax = pCD->custom_data["toxicity_halfmax"]; // 0.4 
	static double relative_max_toxicity = pCD->custom_data["relative_max_toxicity"]; 

	signal = toxin; 
	base_val = pCD->phenotype.death.rates[nApoptosis]; 
	max_val = base_val * relative_max_toxicity; // 100*base_val;

	hill = Hill_response_function( signal , toxicity_halfmax , 1.5 ); 
	phenotype.death.rates[nApoptosis] = base_val + (max_val-base_val)*hill; 
	
	return; 
}

void differentiated_cell_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	// find my cell definition 
	static Cell_Definition* pCD = find_cell_definition( pCell->type_name ); 

	// sample environment 

	static int nR = microenvironment.find_density_index( "resource");
	static int nTox = microenvironment.find_density_index( "toxin");
	static int nDebris = microenvironment.find_density_index( "debris" );
	
	// if dead, release debris
	if( phenotype.death.dead == true )
	{
		phenotype.secretion.net_export_rates[nDebris] = phenotype.volume.total; 
		pCell->functions.update_phenotype = NULL; 
		return;
	}
	
	double * samples = pCell->nearest_density_vector(); 
	double R = samples[nR];
	double toxin = samples[nTox];


	double signal = 0.0; 
	double hill = 0.0; 

	// pressure reduces proliferation 
	signal = pCell->state.simple_pressure;  
	static double pressure_halfmax = pCD->custom_data["cycling_pressure_halfmax"]; // 0.5 
	hill = Hill_response_function( signal, pressure_halfmax , 1.5 );  
	double base_val = pCD->phenotype.cycle.data.exit_rate(0); 

	phenotype.cycle.data.exit_rate(0) = (1-hill)*base_val; 

	// resource reduces necrotic death 

	double max_val = 0.0028;  
	static int nNecrosis = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model );

	// get same from bacteria
	static double necrosis_saturation = pCD->custom_data["necrosis_saturation_resource"]; // 0.075 
	static double necrosis_threshold = pCD->custom_data["necrosis_threshold_resource"]; // 0.15 

	phenotype.death.rates[nNecrosis] = max_val * 
		decreasing_linear_response_function( R, necrosis_saturation, necrosis_threshold ); 

	// toxin increases apoptotic death 
	
	static int nApoptosis = phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model );

	static double toxicity_halfmax = pCD->custom_data["toxicity_halfmax"]; // 0.2 
	static double relative_max_tox_death = pCD->custom_data["relative_max_toxicity"]; // 100 

	signal = toxin; 
	base_val = pCD->phenotype.death.rates[nApoptosis]; 
	double max_response = base_val * relative_max_tox_death; 
	hill = Hill_response_function( signal , toxicity_halfmax , 1.5 ); 
	// std::cout << "tox: " << signal << " " << hill << std::endl; 
	phenotype.death.rates[nApoptosis] = base_val + (max_response-base_val)*hill; 

	return; 
}
