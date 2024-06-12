	
			// std::cout<< "cell id "<< pCell->ID<<"Crossing from right to right rank with rank id "<<world.rank<<<<std::endl;
			//pCell->crossed_to_right_subdomain = false; //RESET IT BUT LATER REMOVE THIS LINE
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

			//std::cout<<"Type Name in Packing = "<<pCell->type_name<<std::endl;
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

				MPI_Pack(&len_str, 1, 					MPI_INT, 	&snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
				MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
				MPI_Pack(&temp_int, 1, 					MPI_INT, 	&snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD );
			}



		 /* Packing members of class 'Variable' (as we are going depth first i.e. Inorder traversal of data members) */

		 	len_vector = pCell->custom_data.variables.size();
		 	len_snd_buf_right = position_right + sizeof(len_vector);
		 	snd_buf_right.resize(len_snd_buf_right);
		 	MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		 	for(int i=0; i<pCell->custom_data.variables.size(); i++)
		 	{
		    temp_str = pCell->custom_data.variables[i].name;
		    len_str = temp_str.length();
		    len_snd_buf_right = position_right + len_str + sizeof(len_str);
		    snd_buf_right.resize(len_snd_buf_right);
		    MPI_Pack(&len_str, 1, 					MPI_INT,  &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		    MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		    /* There was a mistake in MPI_Pack here, instead of &snd_buf_right[0], I had &snd_buf_right */
		    /* A second mistake is the 'missing buffer resize statement ! 														*/

		    temp_double = pCell->custom_data.variables[i].value;
		    len_snd_buf_right = position_right + sizeof(temp_double);
		    snd_buf_right.resize(len_snd_buf_right);
		    MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);


		    temp_str = pCell->custom_data.variables[i].units;
		    len_str = temp_str.length();
		    len_snd_buf_right = position_right + len_str + sizeof(len_str);
		    snd_buf_right.resize(len_snd_buf_right);
		    MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		    MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
			}

			/* Packing members of class 'Vector_Variable' i.e. string name, a vector<double> value and string units */

			len_vector = pCell->custom_data.vector_variables.size();
		 	len_snd_buf_right = position_right + sizeof(len_vector);
		 	snd_buf_right.resize(len_snd_buf_right);
		 	MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);


			for(int i=0; i<pCell->custom_data.vector_variables.size(); i++)
			{
				temp_str = pCell->custom_data.vector_variables[i].name;
		    len_str = temp_str.length();
		    len_snd_buf_right = position_right + len_str + sizeof(len_str);
		    snd_buf_right.resize(len_snd_buf_right);
		    MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		    MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		    len_vector = pCell->custom_data.vector_variables[i].value.size();
		    len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		    snd_buf_right.resize(len_snd_buf_right);
		    MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
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

		 	/* Instead of packing elements one-by-one, I can pack 'len_int' doubles - check at other places also */
		 	for(int j=0; j<len_int; j++)
		 	{
		 	temp_double = pCell->functions.cycle_model.data.transition_rates[i][j];
		 	MPI_Pack(&temp_double, 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		 	}

			/* pCell->functions.cycle_data.data.transition_rates[i] - is INCORRECT, as MPI needs base address of a double array and NOT a vector address */
		 	//MPI_Pack(&(pCell->functions.cycle_model.data.transition_rates[i][0]), len_int, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		 }

		 /* 1 int and 1 double are next to be packed in class Cycle_Data */

		 len_snd_buf_right = position_right + sizeof(int) + sizeof(double);
		 snd_buf_right.resize(len_snd_buf_right);
		 MPI_Pack(&(pCell->functions.cycle_model.data.current_phase_index), 	1, MPI_INT, 	 &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		 MPI_Pack(&(pCell->functions.cycle_model.data.elapsed_time_in_phase), 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		/* Now packing members of class Cell_State as this is the next member in class Cell */
		/* has std::vector<double> orientation and double simple_pressure 									*/
		// std::cout<<"size snd_buf_right "<<len_snd_buf_right<<std::endl;; 
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
		// std::cout<<"size snd_buf_right "<<len_snd_buf_right<<std::endl;; 
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
		// std::cout<<"size snd_buf_right "<<len_snd_buf_right<<std::endl;; 
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
		/* 5 doubles are to be packed, now v1.8 has added 4 more doubles and 1 int */

		len_snd_buf_right = position_right + 9 * sizeof(double) + 1 * sizeof(int);
		snd_buf_right.resize(len_snd_buf_right);

		MPI_Pack(&(pCell->phenotype.mechanics.cell_cell_adhesion_strength), 			 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.mechanics.cell_BM_adhesion_strength), 				 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.mechanics.cell_cell_repulsion_strength), 			 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.mechanics.cell_BM_repulsion_strength), 				 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.mechanics.relative_maximum_adhesion_distance), 1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		MPI_Pack(&(pCell->phenotype.mechanics.relative_maximum_attachment_distance), 	1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.mechanics.relative_detachment_distance), 					1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.mechanics.attachment_elastic_constant), 					1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.mechanics.maximum_attachment_rate), 							1, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		MPI_Pack(&(pCell->phenotype.mechanics.maximum_number_of_attachments), 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);


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

		/*================================================================*/
		/* Version 1.7 introduces 2 new 'int' variables in class Motility */
		/* Now these will be packed 																			*/
		/*================================================================*/

		len_snd_buf_right = position_right + 2 * sizeof(int);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&(pCell->phenotype.motility.chemotaxis_index), 		1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.motility.chemotaxis_direction), 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

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

		/*==================================================================================*/
		/* Version 1.7 introduces 1 new double vector 'net_export_rates' in class Secretion */
		/* Now these will be packed 																												*/
		/*==================================================================================*/

		len_vector = pCell->phenotype.secretion.net_export_rates.size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->phenotype.secretion.net_export_rates[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);


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

		/* Now packing data members of class Intracellular as its object is a data member of class Phenotype 		  */
		/* In PhysiCell v1.9.0, 'type' has changed to 'intracellular_type' => different from physiboss-dev branch */
		
		/* Gaurav Saxena is removing the following code because there might be some problems in 			*/
		/* comparing the string type to NULL (maybe nullptr should be used). The intracellular_type 	*/
		/* need NOT be packed because in the constructor for MaBossIntracellular we can set the 			*/
		/* intracellular_type="maboss" (it is set in all constructors except the one where unpacking) */
		/* is taking place. Thus, I am shifting setting this value while unpacking when constructor		*/
		/* of this class is called 																																		*/

		
		if (pCell->phenotype.intracellular != NULL) 											//this is NULL and NOT nullptr
		{
			temp_str = pCell->phenotype.intracellular->intracellular_type; 				//This was the original code
		}
		else
			temp_str = "not-maboss";
						
		len_str = temp_str.length();
		len_snd_buf_right = position_right + sizeof(len_str) + len_str; 
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_str, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD); 
		MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
			
			
		
#ifdef ADDON_PHYSIBOSS

		if(pCell->phenotype.intracellular != NULL)			//Added by Gaurav Saxena
			if (pCell->phenotype.intracellular->intracellular_type.compare("maboss") == 0) 
			{
				MaBoSSIntracellular* t_intracellular = static_cast<MaBoSSIntracellular*>(pCell->phenotype.intracellular); 
				t_intracellular->pack(snd_buf_right, len_snd_buf_right, position_right);	
			}
				
#endif
		
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

		/*=============================================================================================*/
		/* v1.7 introduces 2 more new protected vector of doubles - helper functions have been written */
		/* Now they will be packed - similar to the way the other 3 protected vector doubles 					 */
		/*=============================================================================================*/

		std::vector<double> & temp_cell_source_sink_solver_temp_export1 = pCell->get_cell_source_sink_solver_temp_export1();
		len_vector = temp_cell_source_sink_solver_temp_export1.size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&temp_cell_source_sink_solver_temp_export1[0], len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		std::vector<double> & temp_cell_source_sink_solver_temp_export2 = pCell->get_cell_source_sink_solver_temp_export2();
		len_vector = temp_cell_source_sink_solver_temp_export2.size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&temp_cell_source_sink_solver_temp_export2[0], len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

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
		MPI_Pack(pCell->secretion_rates->data(), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		len_vector = pCell->saturation_densities->size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(pCell->saturation_densities->data(), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		len_vector = pCell->uptake_rates->size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(pCell->uptake_rates->data(), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		/*==============================================================================================*/
		/* A new public data member of type pointer to a double vector named 'net_export_rates' in v1.7 */
		/* Now this will be packed 																																			*/
		/*==============================================================================================*/

		len_vector = pCell->net_export_rates->size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(pCell->net_export_rates->data(), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);


		len_vector = pCell->internalized_substrates->size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(pCell->internalized_substrates->data(), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		len_vector = pCell->fraction_released_at_death->size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(pCell->fraction_released_at_death->data(), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		len_vector = pCell->fraction_transferred_when_ingested->size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(pCell->fraction_transferred_when_ingested->data(), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);

		/* 1 single int type; to be packed next, remember ID and position are already packed */

		len_snd_buf_right = position_right + sizeof(int);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&(pCell->type), 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);




		len_vector = pCell->velocity.size();
		len_snd_buf_right = position_right + sizeof(len_vector) + len_vector * sizeof(double);
		snd_buf_right.resize(len_snd_buf_right);
		MPI_Pack(&len_vector, 1, MPI_INT, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		MPI_Pack(&(pCell->velocity[0]), len_vector, MPI_DOUBLE, &snd_buf_right[0], len_snd_buf_right, &position_right, MPI_COMM_WORLD);
		if(pCell->phenotype.death.dead == true && pCell->phenotype.cycle.current_phase().code == 14)
		{
			std::cout<<"Packing right to right "<<"Rank= "<<world.rank<<" Cell ID= "<<pCell->ID<<" received IS ZOMBIE"<<std::endl;
			// pCell->print_cell(world);
		}
		int pcellsize=1332;
		if(len_snd_buf_right%pcellsize){
			pCell->print_cell(world);
			 (*all_cells)[i-1]->print_cell(world);
			std::cout<<"Packing phenotype.cycle.current_phase().code "<<pCell->phenotype.cycle.current_phase().code<<std::endl;
		}
		// pCell->print_cell(world);
	