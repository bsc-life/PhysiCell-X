#include "maboss_intracellular.h"
#include <mpi.h>
#include "../core/MPI_helper.h"

MaBoSSIntracellular::MaBoSSIntracellular() : Intracellular()
{
	intracellular_type = "maboss";
	initial_values.clear();
	mutations.clear();
	parameters.clear();
}

MaBoSSIntracellular::MaBoSSIntracellular(pugi::xml_node& node)
{
	intracellular_type = "maboss";
	initialize_intracellular_from_pugixml(node);
}

MaBoSSIntracellular::MaBoSSIntracellular(MaBoSSIntracellular* copy) 
{
	intracellular_type = copy->intracellular_type;
	bnd_filename = copy->bnd_filename;
	cfg_filename = copy->cfg_filename;
	time_step = copy->time_step;
	discrete_time = copy->discrete_time;
	time_tick = copy->time_tick;
	scaling = copy->scaling;
	time_stochasticity = copy->time_stochasticity;
	initial_values = copy->initial_values;
	mutations = copy->mutations;
	parameters = copy->parameters;
	
	if (copy->maboss.has_init()) {
		maboss.init_maboss(copy->bnd_filename, copy->cfg_filename);
		maboss.mutate(mutations);
		maboss.set_initial_values(initial_values);
		maboss.set_parameters(parameters);
		maboss.set_update_time_step(copy->time_step);
		maboss.set_discrete_time(copy->discrete_time, copy->time_tick);
		maboss.set_scaling(copy->scaling);
		maboss.set_time_stochasticity(copy->time_stochasticity);
		maboss.restart_node_values();
		//maboss.set_state(copy->maboss.get_maboss_state());
		//std::cout << get_state();
	}	
}

MaBoSSIntracellular::MaBoSSIntracellular(std::vector<char>& buffer, int& len_buffer, int& position)  
{
	
	intracellular_type = "maboss";

	// double len_str = 0;
	int temp_int;
	double temp_double;
	std::string temp_str;
	int len_str = 0;
		
// intracellular_type = "maboss"; //<---- Vincent ALSO suggested this, but am doing it in unpacking (Physicell_cell.cpp) 
	
	MPI_Unpack(&buffer[0], len_buffer, &position, &len_str, 1, MPI_INT, MPI_COMM_WORLD);
	temp_str.resize(len_str);
	MPI_Unpack(&buffer[0], len_buffer, &position, &temp_str[0], len_str, MPI_CHAR, MPI_COMM_WORLD);
	this->bnd_filename = temp_str;
	
	MPI_Unpack(&buffer[0], len_buffer, &position, &len_str, 1, MPI_INT, MPI_COMM_WORLD);
	temp_str.resize(len_str);
	MPI_Unpack(&buffer[0], len_buffer, &position, &temp_str[0], len_str, MPI_CHAR, MPI_COMM_WORLD);
	this->cfg_filename = temp_str;

	MPI_Unpack(&buffer[0], len_buffer, &position, &(this->time_step), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	
	MPI_Unpack(&buffer[0], len_buffer, &position, &temp_int, 1, MPI_INT, MPI_COMM_WORLD);
	this->discrete_time = (temp_int == 1 ? true : false);  //Gaurav Saxena added parenthesis for clarity

	MPI_Unpack(&buffer[0], len_buffer, &position, &(this->time_tick), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Unpack(&buffer[0], len_buffer, &position, &(this->scaling),  1, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Unpack(&buffer[0], len_buffer, &position, &(this->time_stochasticity),  1, MPI_DOUBLE, MPI_COMM_WORLD);


	this->maboss.init_maboss(this->bnd_filename, this->cfg_filename);
	
	MPI_Unpack(&buffer[0], len_buffer, &position, &temp_int, 1, MPI_INT, MPI_COMM_WORLD);
	
	for (int i=0; i < temp_int; i++) {
		MPI_Unpack(&buffer[0], len_buffer, &position, &len_str, 1, MPI_INT, MPI_COMM_WORLD);
		temp_str.resize(len_str);
		MPI_Unpack(&buffer[0], len_buffer, &position, &temp_str[0], len_str, MPI_CHAR, MPI_COMM_WORLD);
		MPI_Unpack(&buffer[0], len_buffer, &position, &(temp_double), 1, MPI_DOUBLE, MPI_COMM_WORLD);
		this->mutations[temp_str] = temp_double;
	}
	
	this->maboss.mutate(this->mutations);
	
	MPI_Unpack(&buffer[0], len_buffer, &position, &temp_int, 1, MPI_INT, MPI_COMM_WORLD);
	
	for (int i=0; i < temp_int; i++) {
		MPI_Unpack(&buffer[0], len_buffer, &position, &len_str, 1, MPI_INT, MPI_COMM_WORLD);
		temp_str.resize(len_str);
		MPI_Unpack(&buffer[0], len_buffer, &position, &temp_str[0], len_str, MPI_CHAR, MPI_COMM_WORLD);
		MPI_Unpack(&buffer[0], len_buffer, &position, &(temp_double), 1, MPI_DOUBLE, MPI_COMM_WORLD);
		this->initial_values[temp_str] = temp_double;
	}
	this->maboss.set_initial_values(this->initial_values);
	
	MPI_Unpack(&buffer[0], len_buffer, &position, &temp_int, 1, MPI_INT, MPI_COMM_WORLD);
	
	for (int i=0; i < temp_int; i++) {
		MPI_Unpack(&buffer[0], len_buffer, &position, &len_str, 1, MPI_INT, MPI_COMM_WORLD);
		temp_str.resize(len_str);
		MPI_Unpack(&buffer[0], len_buffer, &position, &temp_str[0], len_str, MPI_CHAR, MPI_COMM_WORLD);
		MPI_Unpack(&buffer[0], len_buffer, &position, &(temp_double), 1, MPI_DOUBLE, MPI_COMM_WORLD);
		this->parameters[temp_str] = temp_double;
	}
	this->maboss.set_parameters(this->parameters);
		
	MPI_Unpack(&buffer[0], len_buffer, &position, &(this->next_physiboss_run), 1, MPI_DOUBLE, MPI_COMM_WORLD);

	this->maboss.unpack(buffer, len_buffer, position);
	
}

void MaBoSSIntracellular::pack(std::vector<char>& buffer, int& len_buffer, int& position) 
{
	int len_snd_buf = 0; 
	int	temp_int;
	double temp_double;
	std::string temp_str;
	int len_str = 0;   		//This was unsigned int, Gaurav Saxena changed it to simple int. 
	
	//pack  bnd_filename
	pack_buff(this->bnd_filename, buffer, len_buffer, position);

	//pack cfg_filename
	pack_buff(this->cfg_filename, buffer, len_buffer, position);

	//pack time_step
	pack_buff(this->time_step, buffer, len_buffer, position);
	
	//pack discrete_time
	pack_buff(this->discrete_time, buffer, len_buffer, position);

	//pack time_tick
	pack_buff(this->time_tick, buffer, len_buffer, position);

	//pack scaling
	pack_buff(this->scaling, buffer, len_buffer, position);

	//pack time_stochasticity
	pack_buff(this->time_stochasticity, buffer, len_buffer, position);

	// pack inherit_state
	pack_buff(this->inherit_state, buffer, len_buffer, position);

	// pack inherit nodes
	
	temp_int = static_cast<int>(this->inherit_nodes.size());
	pack_buff(temp_int, buffer, len_buffer, position);

	// inherit nodes
	for (auto t_inherited: this->inherit_nodes) {
		
		pack_buff(t_inherited.first, buffer, len_buffer, position);

		pack_buff(t_inherited.second, buffer, len_buffer, position);

	}

	// pack start_time
	pack_buff(this->start_time, buffer, len_buffer, position);

	// Mutations
	temp_int = static_cast<int>(this->mutations.size());
	pack_buff(temp_int, buffer, len_buffer, position);

	//mutations
	for (auto t_mutation: this->mutations) {
	
		pack_buff(t_mutation.first, buffer, len_buffer, position);

		pack_buff(t_mutation.second, buffer, len_buffer, position);

	}
	
	// Initial values

	temp_int = static_cast<int>(this->initial_values.size());
	pack_buff(temp_int, buffer, len_buffer, position);

	for (auto t_initial_value: this->initial_values) {
		
		pack_buff(t_initial_value.first, buffer, len_buffer, position);

		pack_buff(t_initial_value.second, buffer, len_buffer, position);
	}
	
	// Parameters

	temp_int = static_cast<int>(this->parameters.size());
	pack_buff(temp_int, buffer, len_buffer, position);

	for (auto t_parameter: this->parameters) {
		pack_buff(t_parameter.first, buffer, len_buffer, position);

		pack_buff(t_parameter.second, buffer, len_buffer, position);
	}

	//pack listofinputs

	temp_int = static_cast<int>(this->listOfInputs.size());
	pack_buff(temp_int, buffer, len_buffer, position);

	for (auto t_input: this->listOfInputs) {

		pack_buff(t_input.first, buffer, len_buffer, position);
		
		t_input.second.pack(buffer, len_buffer, position); //class MaBoSS input
	}

	//pack indices of inputs
	
	pack_buff(indicesOfInputs, buffer, len_buffer, position);

	//pack lists of outputs
	temp_int = static_cast<int>(this->listOfOutputs.size());
	pack_buff(temp_int, buffer, len_buffer, position);

	for (auto t_output: this->listOfOutputs) {
		pack_buff(t_output.first, buffer, len_buffer, position);
		
		t_output.second.pack(buffer, len_buffer, position); //class MaBoSS output
	}

	//pack indices of outputs
	pack_buff(indicesOfOutputs, buffer, len_buffer, position);
	
	// Next run
	pack_buff(this->next_physiboss_run, buffer, len_buffer, position);

	//pack maboss that is from the class MaboSSNetwork
	maboss.pack(buffer, len_buffer, position);
	
}

void MaBoSSIntracellular::unpack(const std::vector<char>& buffer, int len_buffer, int& position)
{
    int temp_int;
    double temp_double;
    std::string temp_str;

    // Unpack bnd_filename
    unpack_buff(this->bnd_filename, buffer, len_buffer, position);

    // Unpack cfg_filename
    unpack_buff(this->cfg_filename, buffer, len_buffer, position);

    // Unpack time_step
    unpack_buff(this->time_step, buffer, len_buffer, position);

    // Unpack discrete_time
    unpack_buff(this->discrete_time, buffer, len_buffer, position);

    // Unpack time_tick
    unpack_buff(this->time_tick, buffer, len_buffer, position);

    // Unpack scaling
    unpack_buff(this->scaling, buffer, len_buffer, position);

    // Unpack time_stochasticity
    unpack_buff(this->time_stochasticity, buffer, len_buffer, position);

    // Unpack inherit_state
    unpack_buff(this->inherit_state, buffer, len_buffer, position);

    // Unpack size of inherit_nodes
    unpack_buff(temp_int, buffer, len_buffer, position);
    this->inherit_nodes.clear();
    for (int i = 0; i < temp_int; ++i) {
        std::string key;
        bool value;
        unpack_buff(key, buffer, len_buffer, position);
        unpack_buff(value, buffer, len_buffer, position);
        this->inherit_nodes[key] = value;
    }

    // Unpack start_time
    unpack_buff(this->start_time, buffer, len_buffer, position);

    // Unpack mutations
    unpack_buff(temp_int, buffer, len_buffer, position);
    this->mutations.clear();
    for (int i = 0; i < temp_int; ++i) {
        std::string key;
        int value;
        unpack_buff(key, buffer, len_buffer, position);
        unpack_buff(value, buffer, len_buffer, position);
        this->mutations[key] = value;
    }

    // Unpack initial_values
    unpack_buff(temp_int, buffer, len_buffer, position);
    this->initial_values.clear();
    for (int i = 0; i < temp_int; ++i) {
        std::string key;
        double value;
        unpack_buff(key, buffer, len_buffer, position);
        unpack_buff(value, buffer, len_buffer, position);
        this->initial_values[key] = value;
    }

    // Unpack parameters
    unpack_buff(temp_int, buffer, len_buffer, position);
    this->parameters.clear();
    for (int i = 0; i < temp_int; ++i) {
        std::string key;
        double value;
        unpack_buff(key, buffer, len_buffer, position);
        unpack_buff(value, buffer, len_buffer, position);
        this->parameters[key] = value;
    }

    // Unpack listOfInputs
    unpack_buff(temp_int, buffer, len_buffer, position);
    this->listOfInputs.clear();
    for (int i = 0; i < temp_int; ++i) {
        std::string key;
        MaBoSSInput input;
        unpack_buff(key, buffer, len_buffer, position);
        input.unpack(buffer, len_buffer, position);
        this->listOfInputs[key] = input;
    }

    // Unpack indicesOfInputs
    unpack_buff(this->indicesOfInputs, buffer, len_buffer, position);

    // Unpack listOfOutputs
    unpack_buff(temp_int, buffer, len_buffer, position);
    this->listOfOutputs.clear();
    for (int i = 0; i < temp_int; ++i) {
        std::string key;
        MaBoSSOutput output;
        unpack_buff(key, buffer, len_buffer, position);
        output.unpack(buffer, len_buffer, position);
        this->listOfOutputs[key] = output;
    }

    // Unpack indicesOfOutputs
    unpack_buff(this->indicesOfOutputs, buffer, len_buffer, position);

    // Unpack next_physiboss_run
    unpack_buff(this->next_physiboss_run, buffer, len_buffer, position);

    // Unpack maboss
    this->maboss.unpack(buffer, len_buffer, position); //todo
}

void MaBoSSIntracellular::initialize_intracellular_from_pugixml(pugi::xml_node& node)
{
	pugi::xml_node node_bnd = node.child( "bnd_filename" );
	if ( node_bnd )
	{ bnd_filename = PhysiCell::xml_get_my_string_value (node_bnd); }
	
	pugi::xml_node node_cfg = node.child( "cfg_filename" );
	if ( node_cfg )
	{ cfg_filename = PhysiCell::xml_get_my_string_value (node_cfg); }
	
	pugi::xml_node node_init_values = node.child( "initial_values" );
	if( node_init_values )
	{
		pugi::xml_node node_init_value = node_init_values.child( "initial_value" );
		while( node_init_value )
		{
			std::string node_name = node_init_value.attribute( "node" ).value(); 
			double node_value = PhysiCell::xml_get_my_double_value( node_init_value );
			
			initial_values[node_name] = node_value;
			
			node_init_value = node_init_value.next_sibling( "initial_value" ); 
		}
	}
	
	pugi::xml_node node_mutations = node.child( "mutations" );
	if( node_mutations )
	{
		pugi::xml_node node_mutation = node_mutations.child( "mutation" );
		while( node_mutation )
		{
			std::string node_name = node_mutation.attribute( "node" ).value(); 
			double node_value = PhysiCell::xml_get_my_double_value( node_mutation );
			
			mutations[node_name] = node_value;
			
			node_mutation = node_mutation.next_sibling( "mutation" ); 
		}
	}
	
	pugi::xml_node node_parameters = node.child( "parameters" );
	if( node_parameters )
	{
		pugi::xml_node node_parameter = node_parameters.child( "parameter" );
		while( node_parameter )
		{
			std::string param_name = node_parameter.attribute( "name" ).value(); 
			double param_value = PhysiCell::xml_get_my_double_value( node_parameter );
			
			parameters[param_name] = param_value;
			
			node_parameter = node_parameter.next_sibling( "parameter" ); 
		}
	}
	
	maboss.init_maboss(bnd_filename, cfg_filename);
	maboss.mutate(mutations);
	maboss.set_initial_values(initial_values);
	maboss.set_parameters(parameters);	
	
	pugi::xml_node node_timestep = node.child( "time_step" ); 
	if( node_timestep )
	{ 
		time_step = PhysiCell::xml_get_my_double_value( node_timestep );
		maboss.set_update_time_step(time_step);
	}
	
	pugi::xml_node node_discretetime = node.child( "discrete_time" ); 
	pugi::xml_node node_timetick = node.child( "time_tick" ); 

	if( node_discretetime && node_timetick )
	{ 
		discrete_time = PhysiCell::xml_get_my_bool_value( node_discretetime );		
		time_tick = PhysiCell::xml_get_my_double_value( node_timetick );
		maboss.set_discrete_time(discrete_time, time_tick);
	}

	pugi::xml_node node_scaling = node.child( "scaling" ); 
	if( node_scaling )
	{ 
		scaling = PhysiCell::xml_get_my_double_value( node_scaling );
		maboss.set_scaling(scaling);
	}
	
	pugi::xml_node node_time_stochasticity = node.child( "time_stochasticity" ); 
	if( node_time_stochasticity )
	{ 
		time_stochasticity = PhysiCell::xml_get_my_double_value( node_time_stochasticity );
		maboss.set_time_stochasticity(time_stochasticity);
	}
}

MaBoSSIntracellular* getMaBoSSModel(PhysiCell::Phenotype& phenotype) {
	return static_cast<MaBoSSIntracellular*>(phenotype.intracellular);
}

void MaBoSSIntracellular::save(std::string filename, std::vector<PhysiCell::Cell*>& cells)
{
					
	std::ofstream state_file( filename );
	
	state_file << "ID,state" << std::endl;

	for( auto cell : cells )
		state_file << cell->ID << "," << static_cast<MaBoSSIntracellular*>(cell->phenotype.intracellular)->get_state() << std::endl;
		
	state_file.close();

}