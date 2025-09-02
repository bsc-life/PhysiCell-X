//To consider: USE MPI_Pack_size() to get the size of the packed buffer
void Cell::pack_cell(std::vector<char>& snd_buffer, int& len_buffer, int &position){
	int len_int;
	std::string temp_str;
	int temp_int;		
	
	pack_buff(this->ID, snd_buffer, len_buffer, position);

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

	pack_buff(this->custom_data.variables.size(), snd_buffer, len_buffer, position);
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
	pack_buff(this->custom_data.vector_variables.size(), snd_buffer, len_buffer, position);

    for(int i=0; i<this->custom_data.vector_variables.size(); i++)
    {
		pack_buff(this->custom_data.vector_variables.name, snd_buffer, len_buffer, position);
	

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
	

	pack_buff( inverse_index_maps_cycle_model.size(), snd_buffer, len_buffer, position);

    for(int i=0; i<inverse_index_maps_cycle_model.size(); i++)
    {
		pack_buff( inverse_index_maps_cycle_model[i].size(), snd_buffer, len_buffer, position);

        for(auto it = inverse_index_maps_cycle_model[i].cbegin(); it != inverse_index_maps_cycle_model[i].cend(); it++)
        {
            MPI_Pack(&(it->first),  1, MPI_INT, &snd_buffer[0], len_buffer, &position, MPI_COMM_WORLD);
            MPI_Pack(&(it->second), 1, MPI_INT, &snd_buffer[0], len_buffer, &position, MPI_COMM_WORLD );
        }
    }

	pack_buff( this->functions.cycle_model.name, snd_buffer, len_buffer, position);

	pack_buff( this->functions.cycle_model.code, snd_buffer, len_buffer, position);

	/* Packing of data members of Phases as its object is a member of Cycle_Model */
	/* std::vector<Phase> phases; is a member of it - its a vector */

	
	pack_buff( this->functions.cycle_model.phases.size(), snd_buffer, len_buffer, position);

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

	
	pack_buff( this->functions.cycle_model.phase_links.size(), snd_buffer, len_buffer, position);

    for(int i=0; i<this->functions.cycle_model.phase_links.size(); i++)
    {
		pack_buff(this->functions.cycle_model.phase_links[i].size(), snd_buffer, len_buffer, position);
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

	pack_buff( inverse_index_maps_cycle_data.size(), snd_buffer, len_buffer, position);
	
    for(int i=0; i<inverse_index_maps_cycle_data.size(); i++)
    {
		pack_buff( inverse_index_maps_cycle_data[i].size(), snd_buffer, len_buffer, position);

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

	pack_buff( this->functions.cycle_model.data.transition_rates.size(), snd_buffer, len_buffer, position);
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
	
    /* Packing of members of class Cycle as the class Phenotype contains a data member of this class */
    /* Not packing this->phenotype.cycle.pCycle_Model pointer data member *

    /*=====================================================================*/
    /* Cycle_Model* pCycle_Model;  IS NOT PACKED 													 */
    /*=====================================================================*/

    /* Packing data fields of Cycle_Data as it is an object of Cycle class which is    */
    /* a data member of the Phenotype class which is a member of class Cell						*/
    /* First field is std::vector<std::unordered_map<int,int>> inverse_index_maps;			*/

	std::vector<std::unordered_map<int,int>> &inverse_index_maps_cycle_data_1 = this->phenotype.cycle.data.get_inverse_index_maps();
	
	pack_buff( inverse_index_maps_cycle_data_1.size(), snd_buffer, len_buffer, position);

    for(int i=0; i<inverse_index_maps_cycle_data_1.size(); i++)
    {

		pack_buff( inverse_index_maps_cycle_data_1[i].size(), snd_buffer, len_buffer, position);

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
	pack_buff( this->phenotype.cycle.data.transition_rates.size(), snd_buffer, len_buffer, position);

    for(int i=0; i<this->phenotype.cycle.data.transition_rates.size(); i++)
    {
		pack_buff( this->phenotype.cycle.data.transition_rates[i], snd_buffer, len_buffer, position);
    }

    /* 1 int and 1 double are next to be packed in class Cycle_Data */
	pack_buff( this->phenotype.cycle.data.current_phase_index, snd_buffer, len_buffer, position);
	pack_buff( this->phenotype.cycle.data.elapsed_time_in_phase, snd_buffer, len_buffer, position);
    /* Now packing data members of class Death as there is an object of class Death in Phenotype class */

	pack_buff( this->phenotype.death.rates, snd_buffer, len_buffer, position);

    /* Next field is std::vector<Cycle_Model*> models ;  this is a vector of Cycle_Model * pointers */

    /*=====================================================================*/
    /* std::vector<Cycle_Model*> models ;  IS NOT PACKED  								 */
    /*=====================================================================*/

    /* Now packing members of class Death_Parameters as class Death contains */
    /* std::vector<Death_Parameters> parameters 														 */

	pack_buff( this->phenotype.death.parameters.size(), snd_buffer, len_buffer, position);

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
	pack_buff(this->phenotype.mechanics.dettachment_rate,                        snd_buffer, len_buffer, position);
	
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

		pack_buff(this->secretion_rates, snd_buffer, len_buffer, position);

		pack_buff(this->saturation_densities, snd_buffer, len_buffer, position);

		pack_buff(this->uptake_rates, snd_buffer, len_buffer, position);
		/*==============================================================================================*/
		/* A new public data member of type pointer to a double vector named 'net_export_rates' in v1.7 */
		/* Now this will be packed 																																			*/
		/*==============================================================================================*/

		pack_buff(this->net_export_rates, snd_buffer, len_buffer, position);

		pack_buff(this->internalized_substrates, snd_buffer, len_buffer, position);

		pack_buff(this->fraction_released_at_death, snd_buffer, len_buffer, position);

		pack_buff(this->fraction_transferred_when_ingested, snd_buffer, len_buffer, position);
		/* 1 single int type; to be packed next, remember ID and position are already packed */

		pack_buff(this->type, snd_buffer, len_buffer, position);

		pack_buff(this->velocity, snd_buffer, len_buffer, position);
		//pCell->print_cell(world);
}