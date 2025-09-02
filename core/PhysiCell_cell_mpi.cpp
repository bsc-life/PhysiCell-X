#include "PhysiCell_cell.h"
#include "PhysiCell_cell_container.h"
#include "PhysiCell_utilities.h"
#include "PhysiCell_constants.h"
#include "MPI_helper.h"
#include "../BioFVM/BioFVM_vector.h"


#ifdef ADDON_PHYSIBOSS
#include "../addons/PhysiBoSS/src/maboss_intracellular.h"
#endif
#ifdef ADDON_ROADRUNNER
#include "../addons/libRoadrunner/src/librr_intracellular.h"
#endif
#ifdef ADDON_PHYSIDFBA
#include "../addons/dFBA/src/dfba_intracellular.h"
#endif

#include<limits.h>
#include <signal.h> //for segfault
#include <algorithm>
#include <iterator>

using namespace PhysiCell;

//To consider: USE MPI_Pack_size() to get the size of the packed buffer
void Cell::pack(std::vector<char>& snd_buffer, int& len_buffer, int &position){
	int len_int;
	std::string temp_str;
	int temp_int;		
	temp_int = this->ID;
	pack_buff(temp_int, snd_buffer, len_buffer, position);

    /* Pack double position[0], position[1], position[2] - needed to assign position 							 */
    /* This is the only array where I am not packing the length of the array with the array values */
    /* as we're only dealing with 3-dimensions here, hence position has x, y, z values						 */

	//vector double
	pack_buff(this->position, snd_buffer, len_buffer, position);

	/* Pack length of string (MPI_INT) and then string this->type_name (MPI_CHAR) */

	//std::cout<<"Type Name in Packing = "<<this->type_name<<std::endl;
	pack_buff(this->type_name, snd_buffer, len_buffer, position);

    /* Packing unordered_map<std::string, int> in Custom_Data class */
    /* First pack the number of entries in the unordered_map 				*/

    len_int = this->custom_data.get_name_to_index_map().size();
	pack_buff(len_int, snd_buffer, len_buffer, position);

    /* Now pack the <string, int> pairs but store length(string) before string */

    for(auto it = (this->custom_data.get_name_to_index_map()).cbegin(); it != (this->custom_data.get_name_to_index_map()).cend(); it++)
    {
        /* resize buffer to contain length of string, string and integer */
        pack_buff(it->first, snd_buffer, len_buffer, position);
		pack_buff(it->second, snd_buffer, len_buffer, position);
    }

    /* Packing members of class 'Variable' (as we are going depth first i.e. Inorder traversal of data members) */

	pack_buff(static_cast<int>(this->custom_data.variables.size()), snd_buffer, len_buffer, position);
    for(int i=0; i<this->custom_data.variables.size(); i++)
    {
        
		pack_buff(this->custom_data.variables[i].name, snd_buffer, len_buffer, position);
        /* There was a mistake in MPI_Pack here, instead of &snd_buffer[0], I had &snd_buffer */
        /* A second mistake is the 'missing buffer resize statement ! 														*/
        //Value
		pack_buff(this->custom_data.variables[i].value, snd_buffer, len_buffer, position);

        //units
		pack_buff(this->custom_data.variables[i].units, snd_buffer, len_buffer, position);

        //1.14 includes conserved quantity
		pack_buff(this->custom_data.variables[i].conserved_quantity, snd_buffer, len_buffer, position);
    }

	/* Packing members of class 'Vector_Variable' i.e. string name, a vector<double> value and string units */
	pack_buff(static_cast<int>(this->custom_data.vector_variables.size()), snd_buffer, len_buffer, position);

    for(int i=0; i<this->custom_data.vector_variables.size(); i++)
    {
		pack_buff(this->custom_data.vector_variables[i].name, snd_buffer, len_buffer, position);
	

		pack_buff(this->custom_data.vector_variables[i].value, snd_buffer, len_buffer, position);

		pack_buff(this->custom_data.vector_variables[i].units, snd_buffer, len_buffer, position);

		pack_buff(this->custom_data.vector_variables[i].conserved_quantity, snd_buffer, len_buffer, position);
		 
    }

    /* Packing members of class 'Cell_Parameters' except for pReference_live_phenotype pointer */
    /* There are 9 doubles and 1 int in this class except for the pointer */

    /*=====================================================================*/
    /* pReference_live_phenotype * IS NOT PACKED 													 */
    /*=====================================================================*/

   // Determine buffer size and pack using helper functions
	pack_buff(this->parameters.o2_hypoxic_threshold,      snd_buffer, len_buffer, position);
	pack_buff(this->parameters.o2_hypoxic_response,       snd_buffer, len_buffer, position);
	pack_buff(this->parameters.o2_hypoxic_saturation,     snd_buffer, len_buffer, position);
	pack_buff(this->parameters.o2_proliferation_saturation, snd_buffer, len_buffer, position);
	pack_buff(this->parameters.o2_proliferation_threshold,  snd_buffer, len_buffer, position);
	pack_buff(this->parameters.o2_reference,              snd_buffer, len_buffer, position);
	pack_buff(this->parameters.o2_necrosis_threshold,     snd_buffer, len_buffer, position);
	pack_buff(this->parameters.o2_necrosis_max,           snd_buffer, len_buffer, position);
	pack_buff(this->parameters.max_necrosis_rate,         snd_buffer, len_buffer, position);
	pack_buff(this->parameters.necrosis_type,             snd_buffer, len_buffer, position);

	/* Packing data members of Cell_Functions in depth first manner */
	/* It only has Cycle_Model data member, so pack Cycle_Model data members depth first */

	std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_model = this->functions.cycle_model.get_inverse_index_maps();
	

	pack_buff( static_cast<int>(inverse_index_maps_cycle_model.size()), snd_buffer, len_buffer, position);

    for(int i=0; i<inverse_index_maps_cycle_model.size(); i++)
    {
		pack_buff( static_cast<int>(inverse_index_maps_cycle_model[i].size()), snd_buffer, len_buffer, position);

        for(auto it = inverse_index_maps_cycle_model[i].cbegin(); it != inverse_index_maps_cycle_model[i].cend(); it++)
        {
            //MPI_Pack(&(it->first),  1, MPI_INT, &snd_buffer[0], len_buffer, &position, MPI_COMM_WORLD);
            //MPI_Pack(&(it->second), 1, MPI_INT, &snd_buffer[0], len_buffer, &position, MPI_COMM_WORLD );
            pack_buff( it->first, snd_buffer, len_buffer, position);
            pack_buff( it->second, snd_buffer, len_buffer, position);
        }
    }

	pack_buff( this->functions.cycle_model.name, snd_buffer, len_buffer, position);

	pack_buff( this->functions.cycle_model.code, snd_buffer, len_buffer, position);

	/* Packing of data members of Phases as its object is a member of Cycle_Model */
	/* std::vector<Phase> phases; is a member of it - its a vector */

	
	pack_buff( static_cast<int>(this->functions.cycle_model.phases.size()), snd_buffer, len_buffer, position);

    for(int i=0; i<this->functions.cycle_model.phases.size(); i++)
    {

		pack_buff( this->functions.cycle_model.phases[i].index, snd_buffer, len_buffer, position);
		pack_buff( this->functions.cycle_model.phases[i].code, snd_buffer, len_buffer, position);


		pack_buff( this->functions.cycle_model.phases[i].name, snd_buffer, len_buffer, position);

		pack_buff( this->functions.cycle_model.phases[i].division_at_phase_exit, snd_buffer, len_buffer, position);

		pack_buff( this->functions.cycle_model.phases[i].removal_at_phase_exit, snd_buffer, len_buffer, position);
    }

	/* Packing members of Phase_Link - std::vector<std::vector<Phase_Link>> phase_links is a member of */
	/* class Cycle_Model which is a member of Functions - we're packing in DFS i.e. Inorder traversal	 */

	
	pack_buff( static_cast<int>(this->functions.cycle_model.phase_links.size()), snd_buffer, len_buffer, position);

    for(int i=0; i<this->functions.cycle_model.phase_links.size(); i++)
    {
		pack_buff(static_cast<int>(this->functions.cycle_model.phase_links[i].size()), snd_buffer, len_buffer, position);
        for(int j=0; j<this->functions.cycle_model.phase_links[i].size(); j++)
        {
			pack_buff(this->functions.cycle_model.phase_links[i][j].start_phase_index, snd_buffer, len_buffer, position);
			pack_buff(this->functions.cycle_model.phase_links[i][j].end_phase_index, snd_buffer, len_buffer, position);
			pack_buff(this->functions.cycle_model.phase_links[i][j].fixed_duration, snd_buffer, len_buffer, position);
        }
    }

    /* Going back to packing the remaning data fields of Cycle_Model */
    /* Next field is is 'default_phase_index' of Cycle_Model class */

	pack_buff(this->functions.cycle_model.default_phase_index, snd_buffer, len_buffer, position);

	/* Packing data fields of Cycle_Data as it is an object of Cycle_Model class above */
	/* First field is std::vector<std::unordered_map<int,int>> inverse_index_maps;			*/

    std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_data = this->functions.cycle_model.data.get_inverse_index_maps();

	pack_buff( static_cast<int>(inverse_index_maps_cycle_data.size()), snd_buffer, len_buffer, position);
	
    for(int i=0; i<inverse_index_maps_cycle_data.size(); i++)
    {
		pack_buff( static_cast<int>(inverse_index_maps_cycle_data[i].size()), snd_buffer, len_buffer, position);

        for(auto it = inverse_index_maps_cycle_data[i].cbegin(); it != inverse_index_maps_cycle_data[i].cend(); it++)
        {
			pack_buff( it->first, snd_buffer, len_buffer, position);
			pack_buff( it->second, snd_buffer, len_buffer, position);
        }
    }

    /*=====================================================================*/
    /* Cycle_Model* pCycle_Model;  IS NOT PACKED 						   */
    /*=====================================================================*/

	pack_buff( this->functions.cycle_model.data.time_units, snd_buffer, len_buffer, position);

	pack_buff( static_cast<int>(this->functions.cycle_model.data.transition_rates.size()), snd_buffer, len_buffer, position);
    for(int i=0; i<this->functions.cycle_model.data.transition_rates.size(); i++)
    {
		pack_buff(  this->functions.cycle_model.data.transition_rates[i], snd_buffer, len_buffer, position);

    /* this->functions.cycle_data.data.transition_rates[i] - is INCORRECT, as MPI needs base address of a double array and NOT a vector address */
    //MPI_Pack(&(this->functions.cycle_model.data.transition_rates[i][0]), len_int, MPI_DOUBLE, &snd_buffer[0], len_buffer, &position, MPI_COMM_WORLD);
    }

    /* 1 int and 1 double are next to be packed in class Cycle_Data */
	pack_buff( this->functions.cycle_model.data.current_phase_index, snd_buffer, len_buffer, position);
	pack_buff( this->functions.cycle_model.data.elapsed_time_in_phase, snd_buffer, len_buffer, position);


	/* Now packing members of class Cell_State as this is the next member in class Cell */
	/* has std::vector<double> orientation and double simple_pressure 									*/

	pack_buff( this->state.orientation, snd_buffer, len_buffer, position);
	pack_buff( this->state.simple_pressure, snd_buffer, len_buffer, position);
	pack_buff( this->state.number_of_nuclei, snd_buffer, len_buffer, position);
	pack_buff( this->state.total_attack_time, snd_buffer, len_buffer, position);
	pack_buff( this->state.contact_with_basement_membrane, snd_buffer, len_buffer, position);
	/* Packing data members of class Phenotype next */

	/* First it has 2 boolean variables - allocate sizeof(int) space for both of them */
	
	pack_buff( this->phenotype.flagged_for_division, snd_buffer, len_buffer, position);
    
	pack_buff( this->phenotype.flagged_for_removal, snd_buffer, len_buffer, position);
	//Jose: primer batch es por aqui
    /* Packing of members of class Cycle as the class Phenotype contains a data member of this class */
    /* Not packing this->phenotype.cycle.pCycle_Model pointer data member *

    /*=====================================================================*/
    /* Cycle_Model* pCycle_Model;  IS NOT PACKED 													 */
    /*=====================================================================*/

    /* Packing data fields of Cycle_Data as it is an object of Cycle class which is    */
    /* a data member of the Phenotype class which is a member of class Cell						*/
    /* First field is std::vector<std::unordered_map<int,int>> inverse_index_maps;			*/

	std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_data_1 = this->phenotype.cycle.data.get_inverse_index_maps();
	
	pack_buff( static_cast<int>(inverse_index_maps_cycle_data_1.size()), snd_buffer, len_buffer, position);

    for(int i=0; i<inverse_index_maps_cycle_data_1.size(); i++)
    {

		pack_buff( static_cast<int>(inverse_index_maps_cycle_data_1[i].size()), snd_buffer, len_buffer, position);

        for(auto it = inverse_index_maps_cycle_data_1[i].cbegin(); it != inverse_index_maps_cycle_data_1[i].cend(); it++)
        {
			pack_buff( it->first, snd_buffer, len_buffer, position);
			pack_buff( it->second, snd_buffer, len_buffer, position);
        }
    }

    /*=====================================================================*/
    /* Cycle_Model* pCycle_Model;  IS NOT PACKED  						   */
    /*=====================================================================*/
    //pack string time_units
	pack_buff( this->phenotype.cycle.data.time_units, snd_buffer, len_buffer, position);
    //pack transition_rates
	pack_buff( static_cast<int>(this->phenotype.cycle.data.transition_rates.size()), snd_buffer, len_buffer, position);

    for(int i=0; i<this->phenotype.cycle.data.transition_rates.size(); i++)
    {
		pack_buff( this->phenotype.cycle.data.transition_rates[i], snd_buffer, len_buffer, position);
    }

    /* 1 int and 1 double are next to be packed in class Cycle_Data */
	pack_buff( this->phenotype.cycle.data.current_phase_index, snd_buffer, len_buffer, position);
	pack_buff( this->phenotype.cycle.data.elapsed_time_in_phase, snd_buffer, len_buffer, position);
	//new member of class Cycle -> asymmetric_division
	pack_buff( this->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities, snd_buffer, len_buffer, position);

    /* Now packing data members of class Death as there is an object of class Death in Phenotype class */

	pack_buff( this->phenotype.death.rates, snd_buffer, len_buffer, position);

    /* Next field is std::vector<Cycle_Model*> models ;  this is a vector of Cycle_Model * pointers */

    /*=====================================================================*/
    /* std::vector<Cycle_Model*> models ;  IS NOT PACKED  								 */
    /*=====================================================================*/

    /* Now packing members of class Death_Parameters as class Death contains */
    /* std::vector<Death_Parameters> parameters 														 */
	
	pack_buff( static_cast<int>(this->phenotype.death.parameters.size()), snd_buffer, len_buffer, position);

    for(int i=0; i<this->phenotype.death.parameters.size(); i++)
    {
		pack_buff(this->phenotype.death.parameters[i].time_units, snd_buffer, len_buffer, position);

        pack_buff(this->phenotype.death.parameters[i].unlysed_fluid_change_rate,      snd_buffer, len_buffer, position);
		pack_buff(this->phenotype.death.parameters[i].lysed_fluid_change_rate,        snd_buffer, len_buffer, position);
		pack_buff(this->phenotype.death.parameters[i].cytoplasmic_biomass_change_rate,snd_buffer, len_buffer, position);
		pack_buff(this->phenotype.death.parameters[i].nuclear_biomass_change_rate,    snd_buffer, len_buffer, position);
		pack_buff(this->phenotype.death.parameters[i].calcification_rate,             snd_buffer, len_buffer, position);
		pack_buff(this->phenotype.death.parameters[i].relative_rupture_volume,        snd_buffer, len_buffer, position);

	}

    /* Next 1 bool value (allocate sizeof(int)) and 1 int value are to be packed */

	pack_buff(this->phenotype.death.dead ,snd_buffer, len_buffer, position);

	pack_buff(this->phenotype.death.current_death_model_index ,snd_buffer, len_buffer, position);

    /* Next data members of class Volume are packed as class Phenotype contains its object */
    /* 22 doubles are to be packed */

	pack_buff(this->phenotype.volume.total ,snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.volume.solid ,snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.volume.fluid ,snd_buffer, len_buffer, position);

	pack_buff(this->phenotype.volume.fluid_fraction ,snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.volume.nuclear,snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.volume.nuclear_fluid ,snd_buffer, len_buffer, position);

    pack_buff(this->phenotype.volume.nuclear_solid,                            snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.volume.cytoplasmic,                              snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.volume.cytoplasmic_fluid,                        snd_buffer, len_buffer, position);

	pack_buff(this->phenotype.volume.cytoplasmic_solid,                        snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.volume.calcified_fraction,                       snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.volume.cytoplasmic_to_nuclear_ratio,            snd_buffer, len_buffer, position);

	pack_buff(this->phenotype.volume.rupture_volume,                           snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.volume.cytoplasmic_biomass_change_rate,         snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.volume.nuclear_biomass_change_rate,             snd_buffer, len_buffer, position);

	pack_buff(this->phenotype.volume.fluid_change_rate,                        snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.volume.calcification_rate,                      snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.volume.target_solid_cytoplasmic,                snd_buffer, len_buffer, position);

	pack_buff(this->phenotype.volume.target_solid_nuclear,                    snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.volume.target_fluid_fraction,                   snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.volume.target_cytoplasmic_to_nuclear_ratio,     snd_buffer, len_buffer, position);

	pack_buff(this->phenotype.volume.relative_rupture_volume,                 snd_buffer, len_buffer, position);


    /* Next packing data members of class Geometry as this class has an object in class Phenotype */
    /* 4 doubles are to be packed */

    pack_buff(this->phenotype.geometry.radius,         snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.geometry.nuclear_radius, snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.geometry.surface_area,   snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.geometry.polarity,       snd_buffer, len_buffer, position);
    /* Next packing object of class Mechanics as class Phenotype contains object of this class */
    /* 5 doubles are to be packed, now v1.8 has added 4 more doubles and 1 int */
	pack_buff(this->phenotype.mechanics.cell_cell_adhesion_strength,              snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.mechanics.cell_BM_adhesion_strength,               snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.mechanics.cell_cell_repulsion_strength,            snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.mechanics.cell_BM_repulsion_strength,              snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.mechanics.relative_maximum_adhesion_distance,      snd_buffer, len_buffer, position);
	
	pack_buff(this->phenotype.mechanics.attachment_rate,                         snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.mechanics.detachment_rate,                        snd_buffer, len_buffer, position);
	
	pack_buff(this->phenotype.mechanics.relative_maximum_attachment_distance,    snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.mechanics.relative_detachment_distance,            snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.mechanics.attachment_elastic_constant,             snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.mechanics.maximum_attachment_rate,                 snd_buffer, len_buffer, position);
	
	pack_buff(this->phenotype.mechanics.maximum_number_of_attachments,           snd_buffer, len_buffer, position);
	

    //pack Phenotype::mechanics::cell_adhesion_affinities
    pack_buff(this->phenotype.mechanics.cell_adhesion_affinities, snd_buffer, len_buffer, position);
    /* Next Packing data members of class Motility */
    /* First bool value is to be packed - allocate sizeof(int) */

	pack_buff(this->phenotype.motility.is_motile, snd_buffer, len_buffer, position);

    /* Now pack 2 doubles */

    pack_buff(this->phenotype.motility.persistence_time, snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.motility.migration_speed,  snd_buffer, len_buffer, position);

    pack_buff(this->phenotype.motility.migration_bias_direction, snd_buffer, len_buffer, position);

	pack_buff(this->phenotype.motility.restrict_to_2D, snd_buffer, len_buffer, position);

	pack_buff(this->phenotype.motility.motility_vector, snd_buffer, len_buffer, position);
    /*================================================================*/
    /* Version 1.7 introduces 2 new 'int' variables in class Motility */
    /* Now these will be packed 																			*/
    /*================================================================*/

    pack_buff(this->phenotype.motility.chemotaxis_index, snd_buffer, len_buffer, position);
	pack_buff(this->phenotype.motility.chemotaxis_direction, snd_buffer, len_buffer, position);
	
	//Pack vector<double> chemotactic_sensitivities
	pack_buff(this->phenotype.motility.chemotactic_sensitivities, snd_buffer, len_buffer, position);

    /* Now packing data members of class Secretion as its object is a data member of class Phenotype */
	pack_buff(this->phenotype.secretion.secretion_rates, snd_buffer, len_buffer, position);

	pack_buff(this->phenotype.secretion.uptake_rates, snd_buffer, len_buffer, position);

	pack_buff(this->phenotype.secretion.saturation_densities, snd_buffer, len_buffer, position);

	pack_buff(this->phenotype.secretion.net_export_rates, snd_buffer, len_buffer, position);

    /* Now packing data members of class Molcular as its object is a data member of class Phenotype */

	pack_buff(this->phenotype.molecular.internalized_total_substrates, snd_buffer, len_buffer, position);

	pack_buff(this->phenotype.molecular.fraction_released_at_death, snd_buffer, len_buffer, position);

	pack_buff(this->phenotype.molecular.fraction_transferred_when_ingested, snd_buffer, len_buffer, position);
    /* Now packing data members of class Intracellular as its object is a data member of class Phenotype 		  */
    /* In PhysiCell v1.9.0, 'type' has changed to 'intracellular_type' => different from physiboss-dev branch */
    
    /* Gaurav Saxena is removing the following code because there might be some problems in 			*/
    /* comparing the string type to NULL (maybe nullptr should be used). The intracellular_type 	*/
    /* need NOT be packed because in the constructor for MaBossIntracellular we can set the 			*/
    /* intracellular_type="maboss" (it is set in all constructors except the one where unpacking) */
    /* is taking place. Thus, I am shifting setting this value while unpacking when constructor		*/
    /* of this class is called 																																		*/

    
    if (this->phenotype.intracellular != NULL) 											//this is NULL and NOT nullptr
    {
        temp_str = this->phenotype.intracellular->intracellular_type; 				//This was the original code
    }
    else
        temp_str = "not-maboss";
                
	pack_buff(temp_str, snd_buffer, len_buffer, position);
			
		
#ifdef ADDON_PHYSIBOSS

		if(this->phenotype.intracellular != NULL)			//Added by Gaurav Saxena
			if (this->phenotype.intracellular->intracellular_type.compare("maboss") == 0) 
			{
				MaBoSSIntracellular* t_intracellular = static_cast<MaBoSSIntracellular*>(this->phenotype.intracellular); 
				t_intracellular->pack(snd_buffer, len_buffer, position);	
			}
				
#endif
		
		/* Returning to class Cell to pack the remaining data members : 2 bools + 1 std::vector<double> */

		pack_buff(this->is_out_of_domain, snd_buffer, len_buffer, position);
		
		pack_buff(this->is_movable, snd_buffer, len_buffer, position);

		pack_buff(this->displacement, snd_buffer, len_buffer, position);

		/* Now Packing data members of class Basic_Agent - the parent class of class Cell */
		/* First 2 private data members are to be packed: double and bool 								*/

		pack_buff(this->get_total_volume(), snd_buffer, len_buffer, position);

		pack_buff(this->get_is_volume_changed(), snd_buffer, len_buffer, position);

		/* Now 3 std::vector<double> are to be packed - but are protected, hence write helper functions */


		std::vector<double> & temp_cell_source_sink_solver_temp1 = this->get_cell_source_sink_solver_temp1();
		pack_buff(temp_cell_source_sink_solver_temp1, snd_buffer, len_buffer, position);

		std::vector<double> & temp_cell_source_sink_solver_temp2 = this->get_cell_source_sink_solver_temp2();
		pack_buff(temp_cell_source_sink_solver_temp2, snd_buffer, len_buffer, position);

		/*=============================================================================================*/
		/* v1.7 introduces 2 more new protected vector of doubles - helper functions have been written */
		/* Now they will be packed - similar to the way the other 3 protected vector doubles 					 */
		/*=============================================================================================*/

		std::vector<double> & temp_cell_source_sink_solver_temp_export1 = this->get_cell_source_sink_solver_temp_export1();
		pack_buff(temp_cell_source_sink_solver_temp_export1, snd_buffer, len_buffer, position);

		std::vector<double> & temp_cell_source_sink_solver_temp_export2 = this->get_cell_source_sink_solver_temp_export2();
		pack_buff(temp_cell_source_sink_solver_temp_export2, snd_buffer, len_buffer, position);

		std::vector<double> & temp_previous_velocity = this->get_previous_velocity();
		pack_buff(temp_previous_velocity, snd_buffer, len_buffer, position);

		std::vector<double> & temp_total_extracellular_substrate_change = this->get_total_extracellular_substrate_change();
		pack_buff(temp_total_extracellular_substrate_change, snd_buffer, len_buffer, position);

		/* Now need to pack multiple std::vector<double> *ptr ; type pointers to vectors */

		pack_buff(*(this->secretion_rates), snd_buffer, len_buffer, position);

		pack_buff(*(this->saturation_densities), snd_buffer, len_buffer, position);

		pack_buff(*(this->uptake_rates), snd_buffer, len_buffer, position);
		/*==============================================================================================*/
		/* A new public data member of type pointer to a double vector named 'net_export_rates' in v1.7 */
		/* Now this will be packed 																																			*/
		/*==============================================================================================*/

		pack_buff(*(this->net_export_rates), snd_buffer, len_buffer, position);

		pack_buff(*(this->internalized_substrates), snd_buffer, len_buffer, position);

		pack_buff(*(this->fraction_released_at_death), snd_buffer, len_buffer, position);

		pack_buff(*(this->fraction_transferred_when_ingested), snd_buffer, len_buffer, position);
		/* 1 single int type; to be packed next, remember ID and position are already packed */

		pack_buff(this->type, snd_buffer, len_buffer, position);

		pack_buff(this->velocity, snd_buffer, len_buffer, position);
		//pCell->print_cell(world);
}

void Cell::unpack(std::vector<char>& rcv_buffer, int& len_buffer, int& position) {
    int len_int;
    std::string temp_str;
    int temp_int;

    unpack_buff(this->ID, rcv_buffer, len_buffer, position);

    // vector<double>
    unpack_buff(this->position, rcv_buffer, len_buffer, position);

    // string
    unpack_buff(this->type_name, rcv_buffer, len_buffer, position);

    // Unpack unordered_map<string, int>
    unpack_buff(len_int, rcv_buffer, len_buffer, position);
	if(len_int < 0) {
		std::cout << "Error unpack cell 1" << std::endl;
		return;
	}
    for (int i = 0; i < len_int; ++i) {
        std::string key;
        int value;
        unpack_buff(key, rcv_buffer, len_buffer, position);
        unpack_buff(value, rcv_buffer, len_buffer, position);
        this->custom_data.get_name_to_index_map()[key] = value;
    }

    // Variables
    unpack_buff(len_int, rcv_buffer, len_buffer, position);
	if(len_int < 0) {
		std::cout << "Error unpack cell 2" << std::endl;
		return;
	}
    this->custom_data.variables.resize(len_int);
    for (int i = 0; i < len_int; ++i) {
        unpack_buff(this->custom_data.variables[i].name, rcv_buffer, len_buffer, position);
        unpack_buff(this->custom_data.variables[i].value, rcv_buffer, len_buffer, position);
        unpack_buff(this->custom_data.variables[i].units, rcv_buffer, len_buffer, position);
        unpack_buff(this->custom_data.variables[i].conserved_quantity, rcv_buffer, len_buffer, position);
    }

    // Vector Variables
    unpack_buff(len_int, rcv_buffer, len_buffer, position);
	if(len_int < 0) {
		std::cout << "Error unpack cell 3" << std::endl;
		return;
	}
    this->custom_data.vector_variables.resize(len_int);
    for (int i = 0; i < len_int; ++i) {
        std::string name;
        unpack_buff(name, rcv_buffer, len_buffer, position); // redundant in original pack, likely error
        this->custom_data.vector_variables[i].name = name;

        unpack_buff(this->custom_data.vector_variables[i].value, rcv_buffer, len_buffer, position);
        unpack_buff(this->custom_data.vector_variables[i].units, rcv_buffer, len_buffer, position);
        unpack_buff(this->custom_data.vector_variables[i].conserved_quantity, rcv_buffer, len_buffer, position);
    }

    // Parameters
    unpack_buff(this->parameters.o2_hypoxic_threshold, rcv_buffer, len_buffer, position);
    unpack_buff(this->parameters.o2_hypoxic_response, rcv_buffer, len_buffer, position);
    unpack_buff(this->parameters.o2_hypoxic_saturation, rcv_buffer, len_buffer, position);
    unpack_buff(this->parameters.o2_proliferation_saturation, rcv_buffer, len_buffer, position);
    unpack_buff(this->parameters.o2_proliferation_threshold, rcv_buffer, len_buffer, position);
    unpack_buff(this->parameters.o2_reference, rcv_buffer, len_buffer, position);
    unpack_buff(this->parameters.o2_necrosis_threshold, rcv_buffer, len_buffer, position);
    unpack_buff(this->parameters.o2_necrosis_max, rcv_buffer, len_buffer, position);
    unpack_buff(this->parameters.max_necrosis_rate, rcv_buffer, len_buffer, position);
    unpack_buff(this->parameters.necrosis_type, rcv_buffer, len_buffer, position);

    // Functions.cycle_model.inverse_index_maps
    int inv_map_size;
    unpack_buff(inv_map_size, rcv_buffer, len_buffer, position);
    auto& inverse_index_maps = this->functions.cycle_model.get_inverse_index_maps();
    inverse_index_maps.resize(inv_map_size);
    for (int i = 0; i < inv_map_size; ++i) {
        int map_len;
        unpack_buff(map_len, rcv_buffer, len_buffer, position);
        for (int j = 0; j < map_len; ++j) {
            int key, value;
            unpack_buff(key, rcv_buffer, len_buffer, position);
            unpack_buff(value, rcv_buffer,len_buffer, position);
            inverse_index_maps[i][key] = value;
        }
    }

	unpack_buff(this->functions.cycle_model.name, rcv_buffer, len_buffer, position);
    unpack_buff(this->functions.cycle_model.code, rcv_buffer, len_buffer, position);

	//Phases
	int phase_count;
    unpack_buff(phase_count, rcv_buffer, len_buffer, position);
    this->functions.cycle_model.phases.resize(phase_count);
    for (int i = 0; i < phase_count; ++i) {
        unpack_buff(this->functions.cycle_model.phases[i].index, rcv_buffer, len_buffer, position);
        unpack_buff(this->functions.cycle_model.phases[i].code, rcv_buffer, len_buffer, position);
        unpack_buff(this->functions.cycle_model.phases[i].name, rcv_buffer, len_buffer, position);
        unpack_buff(this->functions.cycle_model.phases[i].division_at_phase_exit, rcv_buffer,len_buffer,  position);
        unpack_buff(this->functions.cycle_model.phases[i].removal_at_phase_exit, rcv_buffer, len_buffer, position);
    }

	// Phase Links
	int link_size;
	unpack_buff(link_size, rcv_buffer, len_buffer, position);
	this->functions.cycle_model.phase_links.resize(link_size);
	for (int i = 0; i < link_size; ++i) {
		int sub_size;
		unpack_buff(sub_size, rcv_buffer, len_buffer, position);
		this->functions.cycle_model.phase_links[i].resize(sub_size);
		for (int j = 0; j < sub_size; ++j) {
			unpack_buff(this->functions.cycle_model.phase_links[i][j].start_phase_index, rcv_buffer, len_buffer, position);
			unpack_buff(this->functions.cycle_model.phase_links[i][j].end_phase_index, rcv_buffer, len_buffer, position);
			unpack_buff(this->functions.cycle_model.phase_links[i][j].fixed_duration, rcv_buffer, len_buffer, position);
		}
	}

	unpack_buff(this->functions.cycle_model.default_phase_index, rcv_buffer, len_buffer, position);

	// Cycle_Data
	auto& cycle_data_maps = this->functions.cycle_model.data.get_inverse_index_maps();
	int data_map_size;
	unpack_buff(data_map_size, rcv_buffer, len_buffer, position);
	cycle_data_maps.resize(data_map_size);
	for (int i = 0; i < data_map_size; ++i) {
		int msize;
		unpack_buff(msize, rcv_buffer, len_buffer, position);
		for (int j = 0; j < msize; ++j) {
			int key, value;
			unpack_buff(key, rcv_buffer, len_buffer, position);
			unpack_buff(value, rcv_buffer, len_buffer, position);
			cycle_data_maps[i][key] = value;
		}
	}

	unpack_buff(this->functions.cycle_model.data.time_units, rcv_buffer, len_buffer, position);

	int tr_size;
    unpack_buff(tr_size, rcv_buffer, len_buffer, position);
    this->functions.cycle_model.data.transition_rates.resize(tr_size);
    for (int i = 0; i < tr_size; ++i) {
        unpack_buff(this->functions.cycle_model.data.transition_rates[i], rcv_buffer, len_buffer, position);
    }

	unpack_buff(this->functions.cycle_model.data.current_phase_index, rcv_buffer, len_buffer,position);
    unpack_buff(this->functions.cycle_model.data.elapsed_time_in_phase, rcv_buffer,  len_buffer,position);

    unpack_buff(this->state.orientation, rcv_buffer, len_buffer, position);
    unpack_buff(this->state.simple_pressure, rcv_buffer, len_buffer, position);
    unpack_buff(this->state.number_of_nuclei, rcv_buffer, len_buffer, position);
    unpack_buff(this->state.total_attack_time, rcv_buffer, len_buffer, position);
    unpack_buff(this->state.contact_with_basement_membrane, rcv_buffer, len_buffer, position);

    unpack_buff(this->phenotype.flagged_for_division, rcv_buffer, len_buffer, position);
    unpack_buff(this->phenotype.flagged_for_removal, rcv_buffer, len_buffer, position);


	std::vector<std::unordered_map<int,int>>& inverse_index_maps_cycle_data_1 = this->phenotype.cycle.data.get_inverse_index_maps();
	inverse_index_maps_cycle_data_1.clear(); // Ensure it's empty before filling

	int outer_size;
	unpack_buff(outer_size, rcv_buffer, len_buffer, position);
	inverse_index_maps_cycle_data_1.resize(outer_size);

	for (size_t i = 0; i < outer_size; i++)
	{
		int inner_size;
		unpack_buff(inner_size, rcv_buffer, len_buffer, position);

		for (size_t j = 0; j < inner_size; j++)
		{
			int key, value;
			unpack_buff(key, rcv_buffer, len_buffer, position);
			unpack_buff(value, rcv_buffer, len_buffer, position);

			inverse_index_maps_cycle_data_1[i][key] = value;
		}
	}

	// Unpack string time_units
	unpack_buff(this->phenotype.cycle.data.time_units, rcv_buffer, len_buffer, position);

	// Unpack transition_rates
	int transition_rates_size;
	unpack_buff(transition_rates_size, rcv_buffer, len_buffer, position);
	this->phenotype.cycle.data.transition_rates.resize(transition_rates_size);

	for (size_t i = 0; i < transition_rates_size; i++)
	{
		unpack_buff(this->phenotype.cycle.data.transition_rates[i], rcv_buffer, len_buffer, position);
	}

	// Unpack current_phase_index (int) and elapsed_time_in_phase (double)
	unpack_buff(this->phenotype.cycle.data.current_phase_index, rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.cycle.data.elapsed_time_in_phase, rcv_buffer, len_buffer, position);

	unpack_buff( this->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities, rcv_buffer, len_buffer, position);

	// Unpack death.rates (std::vector<double>)
	unpack_buff(this->phenotype.death.rates, rcv_buffer, len_buffer, position);

	// Unpack size of phenotype.death.parameters
	int death_parameters_size;
	unpack_buff(death_parameters_size, rcv_buffer, len_buffer, position);
	this->phenotype.death.parameters.resize(death_parameters_size);

	// Unpack each Death_Parameters object
	for (size_t i = 0; i < death_parameters_size; i++)
	{
		// Unpack string time_units
		unpack_buff(this->phenotype.death.parameters[i].time_units, rcv_buffer, len_buffer, position);

		// Unpack doubles
		unpack_buff(this->phenotype.death.parameters[i].unlysed_fluid_change_rate,       rcv_buffer, len_buffer, position);
		unpack_buff(this->phenotype.death.parameters[i].lysed_fluid_change_rate,         rcv_buffer, len_buffer, position);
		unpack_buff(this->phenotype.death.parameters[i].cytoplasmic_biomass_change_rate, rcv_buffer, len_buffer, position);
		unpack_buff(this->phenotype.death.parameters[i].nuclear_biomass_change_rate,     rcv_buffer, len_buffer, position);
		unpack_buff(this->phenotype.death.parameters[i].calcification_rate,              rcv_buffer, len_buffer, position);
		unpack_buff(this->phenotype.death.parameters[i].relative_rupture_volume,         rcv_buffer, len_buffer, position);
	}

	// Unpack bool phenotype.death.dead
	unpack_buff(this->phenotype.death.dead, rcv_buffer, len_buffer, position);

	// Unpack int phenotype.death.current_death_model_index
	unpack_buff(this->phenotype.death.current_death_model_index, rcv_buffer, len_buffer, position);

	// Unpack phenotype.volume fields (all doubles)
	unpack_buff(this->phenotype.volume.total,                              rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.volume.solid,                              rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.volume.fluid,                              rcv_buffer, len_buffer, position);

	unpack_buff(this->phenotype.volume.fluid_fraction,                     rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.volume.nuclear,                            rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.volume.nuclear_fluid,                      rcv_buffer, len_buffer, position);

	unpack_buff(this->phenotype.volume.nuclear_solid,                      rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.volume.cytoplasmic,                        rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.volume.cytoplasmic_fluid,                  rcv_buffer, len_buffer, position);

	unpack_buff(this->phenotype.volume.cytoplasmic_solid,                  rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.volume.calcified_fraction,                 rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.volume.cytoplasmic_to_nuclear_ratio,       rcv_buffer, len_buffer, position);

	unpack_buff(this->phenotype.volume.rupture_volume,                     rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.volume.cytoplasmic_biomass_change_rate,    rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.volume.nuclear_biomass_change_rate,        rcv_buffer, len_buffer, position);

	unpack_buff(this->phenotype.volume.fluid_change_rate,                  rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.volume.calcification_rate,                 rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.volume.target_solid_cytoplasmic,           rcv_buffer, len_buffer, position);

	unpack_buff(this->phenotype.volume.target_solid_nuclear,               rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.volume.target_fluid_fraction,              rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.volume.target_cytoplasmic_to_nuclear_ratio,rcv_buffer, len_buffer, position);

	unpack_buff(this->phenotype.volume.relative_rupture_volume,            rcv_buffer, len_buffer, position);

	// Geometry: 4 doubles
	unpack_buff(this->phenotype.geometry.radius,         rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.geometry.nuclear_radius, rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.geometry.surface_area,   rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.geometry.polarity,       rcv_buffer, len_buffer, position);

	// Mechanics: 5 original + 6 v1.8 new = 11 total (10 doubles + 1 int)
	unpack_buff(this->phenotype.mechanics.cell_cell_adhesion_strength,           rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.mechanics.cell_BM_adhesion_strength,            rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.mechanics.cell_cell_repulsion_strength,         rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.mechanics.cell_BM_repulsion_strength,           rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.mechanics.relative_maximum_adhesion_distance,   rcv_buffer, len_buffer, position);

	// New in v1.8
	unpack_buff(this->phenotype.mechanics.attachment_rate,                      rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.mechanics.detachment_rate,                     rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.mechanics.relative_maximum_attachment_distance, rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.mechanics.relative_detachment_distance,         rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.mechanics.attachment_elastic_constant,          rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.mechanics.maximum_attachment_rate,              rcv_buffer, len_buffer, position);

	// 1 int
	unpack_buff(this->phenotype.mechanics.maximum_number_of_attachments,        rcv_buffer, len_buffer, position);

	// Unpack Phenotype::mechanics::cell_adhesion_affinities
	unpack_buff(this->phenotype.mechanics.cell_adhesion_affinities, rcv_buffer, len_buffer, position);

	// Unpack Motility::is_motile (bool handled as int)
	unpack_buff(this->phenotype.motility.is_motile, rcv_buffer, len_buffer, position);

	// Unpack phenotype.motility.persistence_time and migration_speed
	unpack_buff(this->phenotype.motility.persistence_time, rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.motility.migration_speed,  rcv_buffer, len_buffer, position);

	// Unpack migration_bias_direction (likely a std::vector<double> of size 3)
	unpack_buff(this->phenotype.motility.migration_bias_direction, rcv_buffer, len_buffer, position);

	// Unpack restrict_to_2D (stored as int, unpack into bool)
	unpack_buff(this->phenotype.motility.restrict_to_2D, rcv_buffer, len_buffer, position);

	// Unpack motility_vector (likely a std::vector<double> of size 3)
	unpack_buff(this->phenotype.motility.motility_vector, rcv_buffer, len_buffer, position);

	// Motility
	unpack_buff(this->phenotype.motility.chemotaxis_index, rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.motility.chemotaxis_direction, rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.motility.chemotactic_sensitivities, rcv_buffer, len_buffer, position);

	// Secretion
	unpack_buff(this->phenotype.secretion.secretion_rates, rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.secretion.uptake_rates, rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.secretion.saturation_densities, rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.secretion.net_export_rates, rcv_buffer, len_buffer, position);

	// Molecular
	unpack_buff(this->phenotype.molecular.internalized_total_substrates, rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.molecular.fraction_released_at_death, rcv_buffer, len_buffer, position);
	unpack_buff(this->phenotype.molecular.fraction_transferred_when_ingested, rcv_buffer, len_buffer, position);

	// Unpack intracellular_type string first
	unpack_buff(temp_str, rcv_buffer, len_buffer, position);

	if (temp_str != "not-maboss") {
		// Allocate and create the intracellular object if not already created
		if (this->phenotype.intracellular == NULL) {
			//MaBoSSIntracellular intra = new MaBoSSIntracellular();
			//this->phenotype.intracellular = intra;
		}

	#ifdef ADDON_PHYSIBOSS
		if (this->phenotype.intracellular != NULL) {
			if (temp_str.compare("maboss") == 0) {
				MaBoSSIntracellular* t_intracellular = static_cast<MaBoSSIntracellular*>(this->phenotype.intracellular);
				t_intracellular->unpack(rcv_buffer, len_buffer, position);
			}
			else {
				// Handle other intracellular types or error as needed
			}
		}
	#endif
	}
	else {
		// No intracellular object, set pointer to NULL explicitly
		if (this->phenotype.intracellular != NULL) {
			delete this->phenotype.intracellular;
			this->phenotype.intracellular = NULL;
		}
	}

	unpack_buff(this->is_out_of_domain, rcv_buffer, len_buffer, position);

	unpack_buff(this->is_movable, rcv_buffer, len_buffer, position);

	unpack_buff(this->displacement, rcv_buffer, len_buffer, position);

	double total_volume;
	unpack_buff(total_volume, rcv_buffer, len_buffer, position);
	this->set_total_volume(total_volume);

	bool is_volume_changed;
	unpack_buff(is_volume_changed, rcv_buffer, len_buffer, position);
	this->set_is_volume_changed(is_volume_changed);

	std::vector<double> vec_double;
	unpack_buff(vec_double, rcv_buffer, len_buffer, position);
	this->set_cell_source_sink_solver_temp1(vec_double);

	unpack_buff(vec_double, rcv_buffer, len_buffer, position);
	this->set_cell_source_sink_solver_temp2(vec_double);

	unpack_buff(vec_double, rcv_buffer, len_buffer, position);
	this->set_cell_source_sink_solver_temp_export1(vec_double);

	unpack_buff(vec_double, rcv_buffer, len_buffer, position);
	this->set_cell_source_sink_solver_temp_export2(vec_double);

	unpack_buff(vec_double, rcv_buffer, len_buffer, position);
	//this->set_previous_velocity(vec_double);
	previous_velocity = vec_double;

	unpack_buff(vec_double, rcv_buffer, len_buffer, position);
	this->set_total_extracellular_substrate_change(vec_double);

	unpack_buff(*(this->secretion_rates), rcv_buffer, len_buffer, position);
	unpack_buff(*(this->saturation_densities), rcv_buffer, len_buffer, position);
	unpack_buff(*(this->uptake_rates), rcv_buffer, len_buffer, position);

	unpack_buff(*(this->net_export_rates), rcv_buffer, len_buffer, position);
	unpack_buff(*(this->internalized_substrates), rcv_buffer, len_buffer, position);
	unpack_buff(*(this->fraction_released_at_death), rcv_buffer, len_buffer, position);
	unpack_buff(*(this->fraction_transferred_when_ingested), rcv_buffer, len_buffer, position);
	unpack_buff(this->type, rcv_buffer, len_buffer, position);
	unpack_buff(this->velocity, rcv_buffer, len_buffer, position);
} 

#include <iostream>
#include <iomanip>  // for std::setprecision

//Jose: it does not print all parametes, but enought to know if pack/unpack are coherent
void Cell::print_parameters(std::ofstream& outFile) {
    outFile << std::fixed << std::setprecision(2); // for nicer floating-point printing

    // Basic variables
    outFile << "ID: " << this->ID << "\n";

    outFile << "Position: ";
    for (auto v : this->position) outFile << v << " ";
    outFile << "\n";

    outFile << "Type name: " << this->type_name << "\n";

    // Custom Data - name_to_index map
    /*
    
    outFile << "Custom Data - name_to_index:\n";
    for (const auto& kv : this->custom_data.get_name_to_index_map()) {
        outFile << "  [" << kv.first << "] = " << kv.second << "\n";
    }*/

    // Custom Data Variables (scalar)
    outFile << "Custom Data Variables:\n";
    for (const auto& var : this->custom_data.variables) {
        outFile << "  Name: " << var.name 
                  << ", Value: " << var.value 
                  << ", Units: " << var.units 
                  << ", Conserved: " << (var.conserved_quantity ? "Yes" : "No") << "\n";
    }

    // Custom Data Variables (vector)
    outFile << "Custom Data Vector Variables:\n";
    for (const auto& var : this->custom_data.vector_variables) {
        outFile << "  Name: " << var.name 
                  << ", Units: " << var.units 
                  << ", Conserved: " << (var.conserved_quantity ? "Yes" : "No") 
                  << ", Values: ";
        for (auto v : var.value) outFile << v << " ";
        outFile << "\n";
    }

    // Parameters
    outFile << "Parameters:\n"
              << "  o2_hypoxic_threshold: " << this->parameters.o2_hypoxic_threshold << "\n"
              << "  o2_hypoxic_response: " << this->parameters.o2_hypoxic_response << "\n"
              << "  o2_hypoxic_saturation: " << this->parameters.o2_hypoxic_saturation << "\n"
              << "  o2_proliferation_saturation: " << this->parameters.o2_proliferation_saturation << "\n"
              << "  o2_proliferation_threshold: " << this->parameters.o2_proliferation_threshold << "\n"
              << "  o2_reference: " << this->parameters.o2_reference << "\n"
              << "  o2_necrosis_threshold: " << this->parameters.o2_necrosis_threshold << "\n"
              << "  o2_necrosis_max: " << this->parameters.o2_necrosis_max << "\n"
              << "  max_necrosis_rate: " << this->parameters.max_necrosis_rate << "\n"
              << "  necrosis_type: " << this->parameters.necrosis_type << "\n";

    // Functions - cycle_model
    outFile << "Cycle Model:\n"
              << "  Name: " << this->functions.cycle_model.name << "\n"
              << "  Code: " << this->functions.cycle_model.code << "\n"
              << "  Default Phase Index: " << this->functions.cycle_model.default_phase_index << "\n";

    // Phases
    outFile << "  Phases (" << this->functions.cycle_model.phases.size() << "):\n";
    for (const auto& phase : this->functions.cycle_model.phases) {
        outFile << "    Index: " << phase.index 
                  << ", Name: " << phase.name 
                  << ", Division at phase exit: " << phase.division_at_phase_exit 
                  << ", Removal at phase exit: " << phase.removal_at_phase_exit << "\n";
    }

    // Phase links
    outFile << "  Phase Links (" << this->functions.cycle_model.phase_links.size() << "):\n";
    for (const auto& link : this->functions.cycle_model.phase_links) {
		for (const auto& phase : link)
        outFile << "    " << phase.start_phase_index << " -> " << phase.end_phase_index << "\n";
    }

    // Cycle data
    outFile << "Cycle Data:\n"
              << "  Time units: " << this->phenotype.cycle.data.time_units << "\n"
              << "  Transition rates (" << this->phenotype.cycle.data.transition_rates.size() << "): ";
    for (auto r : this->phenotype.cycle.data.transition_rates) outFile << r << " ";
    outFile << "\n  Current phase index: " << this->phenotype.cycle.data.current_phase_index << "\n"
              << "  Elapsed time in phase: " << this->phenotype.cycle.data.elapsed_time_in_phase << "\n"
			  << "  Asymetric division: ";
	for (const auto& p : phenotype.cycle.asymmetric_division.asymmetric_division_probabilities) {
		outFile << p << " ";
	}
	outFile << "\n";

    /*
    // Inverse index maps of cycle data
    outFile << "  Inverse index maps (" << this->phenotype.cycle.data.get_inverse_index_maps().size() << "):\n";
    for (size_t i = 0; i < this->phenotype.cycle.data.get_inverse_index_maps().size(); i++) {
        const auto& map_i = this->phenotype.cycle.data.get_inverse_index_maps()[i];
        outFile << "    Map " << i << " (" << map_i.size() << " entries):\n";
        for (const auto& kv : map_i) {
            outFile << "      " << kv.first << " => " << kv.second << "\n";
        }
    }*/

    // Death data
    outFile << "Death Data:\n"
              << "  Rates (" << this->phenotype.death.rates.size() << "): ";
    for (auto r : this->phenotype.death.rates) outFile << r << " ";
    outFile << "\n  Parameters (" << this->phenotype.death.parameters.size() << "):\n";

    for (const auto& p : this->phenotype.death.parameters) {
        outFile << "    Time units: " << p.time_units << "\n"
                  << "    Unlysed fluid change rate: " << p.unlysed_fluid_change_rate << "\n"
                  << "    Lysed fluid change rate: " << p.lysed_fluid_change_rate << "\n"
                  << "    Cytoplasmic biomass change rate: " << p.cytoplasmic_biomass_change_rate << "\n"
                  << "    Nuclear biomass change rate: " << p.nuclear_biomass_change_rate << "\n"
                  << "    Calcification rate: " << p.calcification_rate << "\n"
                  << "    Relative rupture volume: " << p.relative_rupture_volume << "\n";
    }

    outFile << "  Dead: " << this->phenotype.death.dead << "\n"
              << "  Current death model index: " << this->phenotype.death.current_death_model_index << "\n";

    // Volume
    outFile << "Volume:\n"
              << "  Total: " << this->phenotype.volume.total << "\n"
              << "  Solid: " << this->phenotype.volume.solid << "\n"
              << "  Fluid: " << this->phenotype.volume.fluid << "\n";

    outFile << "Geometry:\n"
              << "  Radius: " << this->phenotype.geometry.radius << "\n";

    // Mechanics
    outFile << "Mechanics:\n"
              << "  Cell-cell adhesion strength: " << this->phenotype.mechanics.cell_cell_adhesion_strength << "\n"
              << "  Cell adhesion affinities: " << this->phenotype.mechanics.cell_adhesion_affinities << "\n";

    // Motility
    outFile << "Motility:\n"
              << "  Is motile: " << (this->phenotype.motility.is_motile ? "Yes" : "No") << "\n"
              << "  Persistence time: " << this->phenotype.motility.persistence_time << "\n"
              << "  Migration speed: " << this->phenotype.motility.migration_speed << "\n"
              << "  Migration bias direction: " << this->phenotype.motility.migration_bias_direction << "\n"
              << "  Restrict to 2D: " << (this->phenotype.motility.restrict_to_2D ? "Yes" : "No") << "\n"
              << "  Motility vector: " << this->phenotype.motility.motility_vector << "\n";

    // Intracellular
    outFile << "Intracellular type: ";
    if (this->phenotype.intracellular != nullptr) {
        outFile << this->phenotype.intracellular->intracellular_type << "\n";
#ifdef ADDON_PHYSIBOSS
        if (this->phenotype.intracellular->intracellular_type == "maboss") {
            outFile << "  MaBoSS Intracellular object data (custom):\n";
            this->phenotype.intracellular->print();  // Assuming it has a print method
        }
#endif
    } else {
        outFile << "not-maboss\n";
    }

    // Other booleans and vectors
    outFile << "is_out_of_domain: " << (this->is_out_of_domain ? "true" : "false") << "\n";
    outFile << "is_movable: " << (this->is_movable ? "true" : "false") << "\n";

    outFile << "Displacement: ";
    for (auto v : this->displacement) outFile << v << " ";
    outFile << "\n";

    outFile << "Total volume: " << this->get_total_volume() << "\n";
    outFile << "Is volume changed: " << this->get_is_volume_changed() << "\n";

    // Vectors from getters
    outFile << "cell_source_sink_solver_temp_export1: ";
    for (auto v : this->get_cell_source_sink_solver_temp_export1()) outFile << v << " ";
    outFile << "\n";

    outFile << "cell_source_sink_solver_temp_export2: ";
    for (auto v : this->get_cell_source_sink_solver_temp_export2()) outFile << v << " ";
    outFile << "\n";

    outFile << "previous_velocity: ";
    for (auto v : this->get_previous_velocity()) outFile << v << " ";
    outFile << "\n";

    outFile << "total_extracellular_substrate_change: ";
    for (auto v : this->get_total_extracellular_substrate_change()) outFile << v << " ";
    outFile << "\n";

    // Secretion, saturation, uptake
    outFile << "Secretion rates: ";
    for (auto v : *(this->secretion_rates)) outFile << v << " ";
    outFile << "\n";

    outFile << "Saturation densities: ";
    for (auto v : *(this->saturation_densities)) outFile << v << " ";
    outFile << "\n";

    outFile << "Uptake rates: ";
    for (auto v : *(this->uptake_rates)) outFile << v << " ";
    outFile << "\n";

    // Type and velocity
    outFile << "Type: " << this->type << "\n";

    outFile << "Velocity: ";
    for (auto v : this->velocity) outFile << v << " ";
    outFile << "\n";

    // Custom data variables again, just to be exhaustive (optional)
    // ...

    outFile << "------ End of Cell parameters ------\n";
}

#include <cstdlib> // for rand()
#include <ctime>   // for seeding rand()

void Cell::initialize_random() {
    // Seed random number generator
    std::srand(static_cast<unsigned>(std::time(nullptr)));

    // Initialize ID
    this->ID = std::rand();

    // Initialize position (3D vector)
    this->position = {static_cast<double>(std::rand()) / RAND_MAX, 
                      static_cast<double>(std::rand()) / RAND_MAX, 
                      static_cast<double>(std::rand()) / RAND_MAX};

    // Initialize type_name
    this->type_name = "RandomType_" + std::to_string(std::rand() % 100);

    // Initialize custom_data name_to_index_map
    this->custom_data.get_name_to_index_map().clear();
    int map_size = std::rand() % 5 + 1; // Random size between 1 and 5
    for (int i = 0; i < map_size; ++i) {
        this->custom_data.get_name_to_index_map()["Key_" + std::to_string(i)] = std::rand();
    }

    // Initialize custom_data variables
    this->custom_data.variables.clear();
    int var_size = std::rand() % 5 + 1;
    for (int i = 0; i < var_size; ++i) {
        Variable var;
        var.name = "Var_" + std::to_string(i);
        var.value = static_cast<double>(std::rand()) / RAND_MAX;
        var.units = "unit_" + std::to_string(i);
        var.conserved_quantity = std::rand() % 2;
        this->custom_data.variables.push_back(var);
    }

    // Initialize custom_data vector_variables
    this->custom_data.vector_variables.clear();
    int vec_var_size = std::rand() % 5 + 1;
    for (int i = 0; i < vec_var_size; ++i) {
        Vector_Variable vec_var;
        vec_var.name = "VecVar_" + std::to_string(i);
        vec_var.units = "unit_" + std::to_string(i);
        vec_var.conserved_quantity = std::rand() % 2;
        int vec_size = std::rand() % 3 + 1;
        for (int j = 0; j < vec_size; ++j) {
            vec_var.value.push_back(static_cast<double>(std::rand()) / RAND_MAX);
        }
        this->custom_data.vector_variables.push_back(vec_var);
    }

    // Initialize parameters
    this->parameters.o2_hypoxic_threshold = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.o2_hypoxic_response = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.o2_hypoxic_saturation = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.o2_proliferation_saturation = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.o2_proliferation_threshold = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.o2_reference = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.o2_necrosis_threshold = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.o2_necrosis_max = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.max_necrosis_rate = static_cast<double>(std::rand()) / RAND_MAX;
    this->parameters.necrosis_type = std::rand();

    // Initialize functions.cycle_model
    this->functions.cycle_model.name = "CycleModel_" + std::to_string(std::rand() % 100);
    this->functions.cycle_model.code = std::rand();
    this->functions.cycle_model.default_phase_index = std::rand() % 10;
    this->functions.cycle_model.get_inverse_index_maps().clear();
    map_size = std::rand() % 5 + 1; // Random size between 1 and 5
    this->functions.cycle_model.get_inverse_index_maps().resize(map_size);
    for (int i = 0; i < map_size; ++i) {
        this->functions.cycle_model.get_inverse_index_maps()[i][i] = std::rand();
    }

    // Initialize phases
    this->functions.cycle_model.phases.clear();
    int phase_count = std::rand() % 5 + 1;
    for (int i = 0; i < phase_count; ++i) {
        Phase phase;
        phase.index = i;
        phase.code = std::rand();
        phase.name = "Phase_" + std::to_string(i);
        phase.division_at_phase_exit = std::rand() % 2;
        phase.removal_at_phase_exit = std::rand() % 2;
        this->functions.cycle_model.phases.push_back(phase);
    }

    // Initialize phase_links
    this->functions.cycle_model.phase_links.clear();
    int link_count = std::rand() % 3 + 1;
    for (int i = 0; i < link_count; ++i) {
        std::vector<Phase_Link> links;
        int sub_link_count = std::rand() % 3 + 1;
        for (int j = 0; j < sub_link_count; ++j) {
            Phase_Link link;
            link.start_phase_index = std::rand() % 10;
            link.end_phase_index = std::rand() % 10;
            link.fixed_duration = static_cast<double>(std::rand()) / RAND_MAX;
            links.push_back(link);
        }
        this->functions.cycle_model.phase_links.push_back(links);
    }

	std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_model = this->functions.cycle_model.get_inverse_index_maps();
    inverse_index_maps_cycle_model.clear();
    int map_count = std::rand() % 5 + 1; // Random number of maps between 1 and 5
    for (int i = 0; i < map_count; ++i) {
        std::unordered_map<int, int> map;
        int entry_count = std::rand() % 5 + 1; // Random number of entries per map
        for (int j = 0; j < entry_count; ++j) {
            int key = std::rand() % 100;
            int value = std::rand() % 100;
            map[key] = value;
        }
        inverse_index_maps_cycle_model.push_back(map);
    }

    this->functions.cycle_model.data.time_units = "seconds";
    this->functions.cycle_model.data.transition_rates.clear();
    int transition_rate_count = std::rand() % 5 + 1; // Random number of transition rates
    for (int i = 0; i < transition_rate_count; ++i) {
        std::vector<double> transition_rate;
        for (int j = 0; j < 3; ++j) { // Assuming 3 transition rates
            transition_rate.push_back(static_cast<double>(std::rand()) / RAND_MAX);
        }
        this->functions.cycle_model.data.transition_rates.push_back(transition_rate);
    }



    this->functions.cycle_model.data.current_phase_index = std::rand();
    this->functions.cycle_model.data.elapsed_time_in_phase = static_cast<double>(std::rand()) / RAND_MAX;


    this->state.orientation = {static_cast<double>(std::rand()) / RAND_MAX, 
        static_cast<double>(std::rand()) / RAND_MAX, 
        static_cast<double>(std::rand()) / RAND_MAX};
    
    this->state.simple_pressure = static_cast<double>(std::rand()) / RAND_MAX;
    this->state.number_of_nuclei = std::rand() % 5; // Random number of nuclei between 0 and 4
    this->state.total_attack_time = static_cast<double>(std::rand()) / RAND_MAX;
    this->state.contact_with_basement_membrane = std::rand() % 2; // Random boolean

   // Initialize phenotype
    this->phenotype.flagged_for_division = std::rand() % 2;
    this->phenotype.flagged_for_removal = std::rand() % 2;

    this->phenotype.cycle.data.get_inverse_index_maps().clear(); // Ensure it's empty before filling
    map_size = std::rand() % 5 + 1; // Random size between 1 and 5
    for (int i = 0; i < map_size; ++i) {
        std::unordered_map<int, int> map;
        int entry_count = std::rand() % 5 + 1; // Random number of entries per map
        for (int j = 0; j < entry_count; ++j) {
            int key = std::rand() % 100;
            int value = std::rand() % 100;
            map[key] = value;
        }
        this->phenotype.cycle.data.get_inverse_index_maps().push_back(map);
    }

    this->phenotype.cycle.data.time_units = "seconds";
    this->phenotype.cycle.data.transition_rates.clear();
    transition_rate_count = std::rand() % 5 + 1; // Random number of transition rates
    for (int i = 0; i < transition_rate_count; ++i) {
        std::vector<double> transition_rate;
        for (int j = 0; j < 3; ++j) { // Assuming 10 transition rates
            transition_rate.push_back(static_cast<double>(std::rand()) / RAND_MAX);
        }
        this->phenotype.cycle.data.transition_rates.push_back(transition_rate);
    }
    this->phenotype.cycle.data.current_phase_index = std::rand() % 10; // Random phase index
    this->phenotype.cycle.data.elapsed_time_in_phase = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities.clear();
    int asym_div_prob_count = std::rand() % 5 + 1; // Random number of probabilities
    for (int i = 0; i < asym_div_prob_count; ++i) {
        this->phenotype.cycle.asymmetric_division.asymmetric_division_probabilities.push_back(static_cast<double>(std::rand()) / RAND_MAX);
    }
    
    this->phenotype.death.rates.clear();
    int death_rate_count = std::rand() % 5 + 1; // Random number of death rates
    for (int i = 0; i < death_rate_count; ++i) {
        this->phenotype.death.rates.push_back(static_cast<double>(std::rand()) / RAND_MAX);
    }
    this->phenotype.death.parameters.clear();
    int death_param_count = std::rand() % 5 + 1; // Random number of death parameters
    for (int i = 0; i < death_param_count; ++i) {
        Death_Parameters param;
        param.time_units = "seconds";
        param.unlysed_fluid_change_rate = static_cast<double>(std::rand()) / RAND_MAX;
        param.lysed_fluid_change_rate = static_cast<double>(std::rand()) / RAND_MAX;
        param.cytoplasmic_biomass_change_rate = static_cast<double>(std::rand()) / RAND_MAX;
        param.nuclear_biomass_change_rate = static_cast<double>(std::rand()) / RAND_MAX;
        param.calcification_rate = static_cast<double>(std::rand()) / RAND_MAX;
        param.relative_rupture_volume = static_cast<double>(std::rand()) / RAND_MAX;
        this->phenotype.death.parameters.push_back(param);
    }

    this->phenotype.death.dead = std::rand() % 2; // Random boolean
    this->phenotype.death.current_death_model_index = std::rand() % 10; // Random index
    this->phenotype.volume.total = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.solid = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.fluid = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.fluid_fraction = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.nuclear = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.nuclear_fluid = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.nuclear_solid = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.cytoplasmic = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.cytoplasmic_fluid = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.cytoplasmic_solid = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.volume.calcified_fraction = static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.volume.cytoplasmic_to_nuclear_ratio= static_cast<double>(std::rand()) / RAND_MAX;

	this->phenotype.volume.rupture_volume= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.volume.cytoplasmic_biomass_change_rate= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.volume.nuclear_biomass_change_rate= static_cast<double>(std::rand()) / RAND_MAX;

	this->phenotype.volume.fluid_change_rate= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.volume.calcification_rate= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.volume.target_solid_cytoplasmic= static_cast<double>(std::rand()) / RAND_MAX;

	this->phenotype.volume.target_solid_nuclear= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.volume.target_fluid_fraction= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.volume.target_cytoplasmic_to_nuclear_ratio= static_cast<double>(std::rand()) / RAND_MAX;

	this->phenotype.volume.relative_rupture_volume= static_cast<double>(std::rand()) / RAND_MAX;


    this->phenotype.geometry.radius= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.geometry.nuclear_radius= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.geometry.surface_area= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.geometry.polarity= static_cast<double>(std::rand()) / RAND_MAX;

    this->phenotype.mechanics.cell_cell_adhesion_strength= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.mechanics.cell_BM_adhesion_strength= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.mechanics.cell_cell_repulsion_strength= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.mechanics.cell_BM_repulsion_strength= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.mechanics.relative_maximum_adhesion_distance= static_cast<double>(std::rand()) / RAND_MAX;

    this->phenotype.mechanics.attachment_rate= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.mechanics.detachment_rate= static_cast<double>(std::rand()) / RAND_MAX;
	
	this->phenotype.mechanics.relative_maximum_attachment_distance= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.mechanics.relative_detachment_distance= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.mechanics.attachment_elastic_constant= static_cast<double>(std::rand()) / RAND_MAX;
	this->phenotype.mechanics.maximum_attachment_rate= static_cast<double>(std::rand()) / RAND_MAX;
	
	this->phenotype.mechanics.maximum_number_of_attachments= std::rand();

    this->phenotype.mechanics.cell_adhesion_affinities.clear();
    int adhesion_affinity_count = std::rand() % 5 + 1; //
    for (int i = 0; i < adhesion_affinity_count; ++i) {
        this->phenotype.mechanics.cell_adhesion_affinities.push_back(static_cast<double>(std::rand()) / RAND_MAX);
    }

	this->phenotype.motility.is_motile = std::rand() % 2; // Random boolean
    this->phenotype.motility.persistence_time = static_cast<double>(std::rand()) / RAND_MAX;
    this->phenotype.motility.migration_speed = static_cast<double>(std::rand()) / RAND_MAX;


    this->phenotype.motility.migration_bias_direction = {static_cast<double>(std::rand()) / RAND_MAX, 
        static_cast<double>(std::rand()) / RAND_MAX, 
        static_cast<double>(std::rand()) / RAND_MAX};
    this->phenotype.motility.restrict_to_2D = std::rand() % 2; // Random boolean
    this->phenotype.motility.motility_vector = {static_cast<double>(std::rand()) / RAND_MAX, 
        static_cast<double>(std::rand()) / RAND_MAX, 
        static_cast<double>(std::rand()) / RAND_MAX};
    this->phenotype.motility.chemotaxis_index = std::rand() % 10; // Random index
    this->phenotype.motility.chemotaxis_direction = std::rand() % 10;
    this->phenotype.motility.chemotactic_sensitivities.clear();
    int chemotactic_sensitivity_count = std::rand() % 5 + 1; // Random number of sensitivities
    for (int i = 0; i < chemotactic_sensitivity_count; ++i) {
        this->phenotype.motility.chemotactic_sensitivities.push_back(static_cast<double>(std::rand()) / RAND_MAX);
    }

    /* Now packing data members of class Secretion as its object is a data member of class Phenotype */


	// Assume some random size for the vectors, e.g., 5 to 10 elements
    int substrate_count = std::rand() % 6 + 5; // random number between 5 and 10

    // Secretion rates
    this->phenotype.secretion.secretion_rates.clear();
    this->phenotype.secretion.secretion_rates.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        this->phenotype.secretion.secretion_rates[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // Uptake rates
    this->phenotype.secretion.uptake_rates.clear();
    this->phenotype.secretion.uptake_rates.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        this->phenotype.secretion.uptake_rates[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // Saturation densities
    this->phenotype.secretion.saturation_densities.clear();
    this->phenotype.secretion.saturation_densities.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        this->phenotype.secretion.saturation_densities[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // Net export rates
    this->phenotype.secretion.net_export_rates.clear();
    this->phenotype.secretion.net_export_rates.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        this->phenotype.secretion.net_export_rates[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // Internalized substrates
    this->phenotype.molecular.internalized_total_substrates.clear();
    this->phenotype.molecular.internalized_total_substrates.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        this->phenotype.molecular.internalized_total_substrates[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // Fraction released at death
    this->phenotype.molecular.fraction_released_at_death.clear();
    this->phenotype.molecular.fraction_released_at_death.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        this->phenotype.molecular.fraction_released_at_death[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // Fraction transferred when ingested
    this->phenotype.molecular.fraction_transferred_when_ingested.clear();
    this->phenotype.molecular.fraction_transferred_when_ingested.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        this->phenotype.molecular.fraction_transferred_when_ingested[i] = static_cast<double>(std::rand()) / RAND_MAX;

    //this->phenotype.intracellular = new Intracellular();
    //this->phenotype.intracellular->intracellular_type = "Random";

    // Random boolean for is_out_of_domain (true or false)
    this->is_out_of_domain = std::rand() % 2;

    // Random boolean for is_movable
    this->is_movable = std::rand() % 2;

    // Random 3D vector for displacement
    this->displacement = {
        static_cast<double>(std::rand()) / RAND_MAX,
        static_cast<double>(std::rand()) / RAND_MAX,
        static_cast<double>(std::rand()) / RAND_MAX
    };

    this->set_total_volume(static_cast<double>(std::rand()) / RAND_MAX);
    this->set_is_volume_changed(std::rand() % 2);

    int temp_vector_size = std::rand() % 6 + 5; // Random size between 5 and 10

    // Initialize temp1
    std::vector<double>& temp_cell_source_sink_solver_temp1 = this->get_cell_source_sink_solver_temp1();
    temp_cell_source_sink_solver_temp1.clear();
    temp_cell_source_sink_solver_temp1.resize(temp_vector_size);
    for (int i = 0; i < temp_vector_size; ++i) {
        temp_cell_source_sink_solver_temp1[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    // Initialize temp2
    std::vector<double>& temp_cell_source_sink_solver_temp2 = this->get_cell_source_sink_solver_temp2();
    temp_cell_source_sink_solver_temp2.clear();
    temp_cell_source_sink_solver_temp2.resize(temp_vector_size);
    for (int i = 0; i < temp_vector_size; ++i) {
        temp_cell_source_sink_solver_temp2[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    // Assume substrate-related vectors have N elements
    substrate_count = std::rand() % 6 + 5; // Or replace with known substrate count

    // temp_export1
    std::vector<double>& temp_cell_source_sink_solver_temp_export1 = this->get_cell_source_sink_solver_temp_export1();
    temp_cell_source_sink_solver_temp_export1.clear();
    temp_cell_source_sink_solver_temp_export1.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i) {
        temp_cell_source_sink_solver_temp_export1[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    // temp_export2
    std::vector<double>& temp_cell_source_sink_solver_temp_export2 = this->get_cell_source_sink_solver_temp_export2();
    temp_cell_source_sink_solver_temp_export2.clear();
    temp_cell_source_sink_solver_temp_export2.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i) {
        temp_cell_source_sink_solver_temp_export2[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    // previous_velocity (3D vector)
    std::vector<double>& temp_previous_velocity = this->get_previous_velocity();
    temp_previous_velocity.clear();
    temp_previous_velocity.resize(3);
    for (int i = 0; i < 3; ++i) {
        temp_previous_velocity[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    // total_extracellular_substrate_change
    std::vector<double>& temp_total_extracellular_substrate_change = this->get_total_extracellular_substrate_change();
    temp_total_extracellular_substrate_change.clear();
    temp_total_extracellular_substrate_change.resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i) {
        temp_total_extracellular_substrate_change[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    substrate_count = std::rand() % 6 + 5; // Or use a known substrate count

    // Initialize secretion_rates
    if (!this->secretion_rates) {
        this->secretion_rates = new std::vector<double>();
    }
    this->secretion_rates->clear();
    this->secretion_rates->resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i) {
        (*(this->secretion_rates))[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    // Initialize saturation_densities
    if (!this->saturation_densities) {
        this->saturation_densities = new std::vector<double>();
    }
    this->saturation_densities->clear();
    this->saturation_densities->resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i) {
        (*(this->saturation_densities))[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    // Initialize uptake_rates
    if (!this->uptake_rates) {
        this->uptake_rates = new std::vector<double>();
    }
    this->uptake_rates->clear();
    this->uptake_rates->resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i) {
        (*(this->uptake_rates))[i] = static_cast<double>(std::rand()) / RAND_MAX;
    }

    substrate_count = std::rand() % 6 + 5; // Or a known constant

    // net_export_rates
    if (!this->net_export_rates)
        this->net_export_rates = new std::vector<double>();
    this->net_export_rates->clear();
    this->net_export_rates->resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        (*this->net_export_rates)[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // internalized_substrates
    if (!this->internalized_substrates)
        this->internalized_substrates = new std::vector<double>();
    this->internalized_substrates->clear();
    this->internalized_substrates->resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        (*this->internalized_substrates)[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // fraction_released_at_death
    if (!this->fraction_released_at_death)
        this->fraction_released_at_death = new std::vector<double>();
    this->fraction_released_at_death->clear();
    this->fraction_released_at_death->resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        (*this->fraction_released_at_death)[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // fraction_transferred_when_ingested
    if (!this->fraction_transferred_when_ingested)
        this->fraction_transferred_when_ingested = new std::vector<double>();
    this->fraction_transferred_when_ingested->clear();
    this->fraction_transferred_when_ingested->resize(substrate_count);
    for (int i = 0; i < substrate_count; ++i)
        (*this->fraction_transferred_when_ingested)[i] = static_cast<double>(std::rand()) / RAND_MAX;

    // type (random int type identifier, 0–9)
    this->type = std::rand() % 10;

    // velocity (3D vector)
    this->velocity = {
        static_cast<double>(std::rand()) / RAND_MAX,
        static_cast<double>(std::rand()) / RAND_MAX,
        static_cast<double>(std::rand()) / RAND_MAX
    };

}

