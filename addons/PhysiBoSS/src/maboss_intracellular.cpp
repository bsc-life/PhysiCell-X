#include "maboss_intracellular.h"
#include <mpi.h>

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
		maboss.restart_node_values();
		//maboss.set_state(copy->maboss.get_maboss_state());
		//std::cout << get_state();
	}	
}

MaBoSSIntracellular::MaBoSSIntracellular(std::vector<char>& buffer, int& len_buffer, int& position) {
	
	// double len_str = 0;
	int temp_int;
	double temp_double;
	std::string temp_str;
	int len_str = 0;
	
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
	this->discrete_time = temp_int == 1 ? true : false; 

	MPI_Unpack(&buffer[0], len_buffer, &position, &(this->time_tick), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Unpack(&buffer[0], len_buffer, &position, &(this->scaling),  1, MPI_DOUBLE, MPI_COMM_WORLD);


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
	
	MPI_Unpack(&buffer[0], len_buffer, &position, &(temp_double), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	this->maboss.set_time_to_update(temp_double);

	MPI_Unpack(&buffer[0], len_buffer, &position, &temp_int, 1, MPI_INT, MPI_COMM_WORLD);

	std::vector<Node*> t_nodes = this->maboss.getNetwork()->getNodes();

	for (unsigned int i=0; i < temp_int; i++) {
		int t_node = 0;
		MPI_Unpack(&buffer[0], len_buffer, &position, &t_node, 1, MPI_INT, MPI_COMM_WORLD);
		
		this->maboss.state.setNodeState(t_nodes[i], t_node == 1?true:false);
	}
	
}

void MaBoSSIntracellular::pack(std::vector<char>& buffer, int& len_buffer, int& position) 
{
	int len_snd_buf = 0; 
	int	temp_int;
	double temp_double;
	std::string temp_str;
	unsigned int len_str = 0;
	
	temp_str = this->bnd_filename;
	len_str = temp_str.length();
	len_buffer = position + sizeof(len_str) + len_str; 
	buffer.resize(len_buffer);
	MPI_Pack(&len_str, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD); 
	MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

	temp_str = this->cfg_filename;
	len_str = temp_str.length();
	len_buffer = position + sizeof(len_str) + len_str; 
	buffer.resize(len_buffer);
	MPI_Pack(&len_str, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD); 
	MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

	len_buffer = position + sizeof(double);		
	buffer.resize(len_buffer);
	MPI_Pack(&(this->time_step), 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
	
	len_buffer = position + sizeof(int);		
	buffer.resize(len_buffer);
	temp_int = (this->discrete_time == true) ? 1 : 0;
	// 	temp_int = 1; 
	// else
	// 	temp_int = 0; 
	MPI_Pack(&temp_int, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

	len_buffer = position + sizeof(double);		
	buffer.resize(len_buffer);
	MPI_Pack(&(this->time_tick), 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

	len_buffer = position + sizeof(double);		
	buffer.resize(len_buffer);
	MPI_Pack(&(this->scaling), 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

	// Mutations
	len_buffer = position + sizeof(int);		
	buffer.resize(len_buffer);
	temp_int = this->mutations.size();
	MPI_Pack(&(temp_int), 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

	for (auto t_mutation: this->mutations) {
		
		temp_str = t_mutation.first;
		len_str = temp_str.length();
		len_buffer = position + sizeof(len_str) + len_str; 
		buffer.resize(len_buffer);
		MPI_Pack(&len_str, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD); 
		MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

		len_buffer = position + sizeof(double);		
		buffer.resize(len_buffer);
		MPI_Pack(&(t_mutation.second), 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

	}
	
	// Initial values
	len_buffer = position + sizeof(int);		
	buffer.resize(len_buffer);
	temp_int = this->initial_values.size();
	MPI_Pack(&(temp_int), 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

	for (auto t_initial_value: this->initial_values) {
		
		temp_str = t_initial_value.first;
		len_str = temp_str.length();
		len_buffer = position + sizeof(len_str) + len_str; 
		buffer.resize(len_buffer);
		MPI_Pack(&len_str, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD); 
		MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

		len_buffer = position + sizeof(double);		
		buffer.resize(len_buffer);
		MPI_Pack(&(t_initial_value.second), 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

	}
	
	// Parameters
	len_buffer = position + sizeof(int);		
	buffer.resize(len_buffer);
	temp_int = this->parameters.size();
	MPI_Pack(&(temp_int), 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

	for (auto t_parameter: this->parameters) {
		
		temp_str = t_parameter.first;
		len_str = temp_str.length();
		len_buffer = position + sizeof(len_str) + len_str; 
		buffer.resize(len_buffer);
		MPI_Pack(&len_str, 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD); 
		MPI_Pack(&temp_str[0], len_str, MPI_CHAR, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

		len_buffer = position + sizeof(double);		
		buffer.resize(len_buffer);
		MPI_Pack(&(t_parameter.second), 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

	}
	
	// Next run
	len_buffer = position + sizeof(double);		
	buffer.resize(len_buffer);
	MPI_Pack(&(this->next_physiboss_run), 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

	// MaBoSS's time to update
	len_buffer = position + sizeof(double);		
	buffer.resize(len_buffer);
	temp_double = this->maboss.get_time_to_update();
	MPI_Pack(&(temp_double), 1, MPI_DOUBLE, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
	// this->maboss.set_time_to_update(time_to_update);
	
	std::vector<Node*> t_nodes = this->maboss.getNetwork()->getNodes();
	len_buffer = position + sizeof(int);		
	buffer.resize(len_buffer);
	temp_int = t_nodes.size();
	MPI_Pack(&(temp_int), 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);

	for (unsigned int i=0; i < t_nodes.size(); i++) {
		len_buffer = position + sizeof(int);		
		buffer.resize(len_buffer);
		temp_int = this->maboss.state.getNodeState(t_nodes[i]) == true ? 1 : 0;
		MPI_Pack(&(temp_int), 1, MPI_INT, &buffer[0], len_buffer, &position, MPI_COMM_WORLD);
	}

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