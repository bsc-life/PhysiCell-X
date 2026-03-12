#ifndef _MaBoSS_Intracellular_h_
#define _MaBoSS_Intracellular_h_

#include <string>
#include <map>
#include "../../../core/PhysiCell.h"
#include "../../../core/PhysiCell_phenotype.h"
#include "../../../core/PhysiCell_cell.h"
#include "../../../modules/PhysiCell_pugixml.h"
#include "maboss_network.h"
#include "utils.h"

static std::string PhysiBoSS_Version = "2.2.3"; 
static std::string PhysiBoSS_DOI = "10.1038/s41540-023-00314-4"; 
static std::string PhysiBoSS_URL = "https://github.com/PhysiBoSS/PhysiBoSS"; 

class MaBoSSIntracellular : public PhysiCell::Intracellular { //UPDATED
 private:
 public:

	static long counter;
	
	std::string bnd_filename; //packed
	std::string cfg_filename; //packed
	
	double time_step = 12; //packed
	bool discrete_time = false; //packed
	double time_tick = 0.5; //packed
	double scaling = 1.0; //packed
	double time_stochasticity = 0.0; //packed
	bool inherit_state = false; //new
	std::map<std::string, bool> inherit_nodes; //new
	double start_time = 0.0; //new
	
	std::map<std::string, double> initial_values; //packed
	std::map<std::string, double> mutations; //packed
	std::map<std::string, double> parameters; //packed
	
	std::map<std::string, MaBoSSInput> listOfInputs; //new
	std::vector<int> indicesOfInputs; //new
	std::map<std::string, MaBoSSOutput> listOfOutputs; //new
	std::vector<int> indicesOfOutputs; //new
	MaBoSSNetwork maboss; //packed

	double next_physiboss_run = 0; //packed

	MaBoSSIntracellular();
	
	MaBoSSIntracellular(pugi::xml_node& node);
	
	MaBoSSIntracellular(MaBoSSIntracellular* copy);
	
	MaBoSSIntracellular(std::vector<char>& buffer, int& len_buffer, int& position);
	void pack(std::vector<char>& buffer, int& len_buffer, int& position);
	void unpack(std::vector<char>& buffer, int len_buffer, int& position);
	
	Intracellular* clone() {
		return static_cast<Intracellular*>(new MaBoSSIntracellular(this));
	}
	Intracellular* getIntracellularModel() {
		return static_cast<Intracellular*>(this);
	}
	
	void initialize_intracellular_from_pugixml(pugi::xml_node& node);
	
	void start() {
		this->maboss.restart_node_values();
		this->next_physiboss_run = std::max(this->start_time, PhysiCell::PhysiCell_globals.current_time);
	}
	
	void update() {
		this->maboss.run_simulation();
		this->next_physiboss_run += this->maboss.get_time_to_update();
	}

	void update(PhysiCell::Cell * cell, PhysiCell::Phenotype& phenotype, double dt) {
		this->update_inputs(cell, phenotype, dt);
		this->maboss.run_simulation();
		this->update_outputs(cell, phenotype, dt);
		this->next_physiboss_run += this->maboss.get_time_to_update();
	}
	
	bool need_update() {
		return PhysiCell::PhysiCell_globals.current_time >= this->next_physiboss_run;
	}
	
	void inherit(PhysiCell::Cell * cell) {
		maboss.inherit_state(
			static_cast<MaBoSSIntracellular*>(cell->phenotype.intracellular)->maboss.get_maboss_state(), 
			inherit_state, inherit_nodes
		);
	}
	void update_inputs(PhysiCell::Cell* cell, PhysiCell::Phenotype& phenotype, double dt);
	void update_outputs(PhysiCell::Cell * cell, PhysiCell::Phenotype& phenotype, double dt);

	bool has_variable(std::string name) {
		return this->maboss.has_node(name);
	}
	
	bool get_boolean_variable_value(std::string name) {
		return this->maboss.get_node_value(name);
	}
	
	//The following 2 functions MAY not be needed but since Vincent defined them in 
	//physiboss-dev branch, I am incorporating them
	
	void set_double_variable_value(std::string name, double value) {}
	
	double get_double_variable_value(std::string name) {
		return 0.0;
	}

	void set_boolean_variable_value(std::string name, bool value) {
		this->maboss.set_node_value(name, value);
	}
	
	double get_parameter_value(std::string name) {
		return this->maboss.get_parameter_value(name);
	}
	
	void set_parameter_value(std::string name, double value) {
		this->maboss.set_parameter_value(name, value);
	}
	
	std::string get_state() {
		return this->maboss.get_state();
	}
	
	void print_current_nodes(){
		this->maboss.print_nodes();
	}

	void display(std::ostream& os);
	
	std::string get_bnd_filename() { return this->bnd_filename; }
	std::string get_cfg_filename() { return this->cfg_filename; }	
	static void save(std::string filename);

    // unneeded for this type
    int update_phenotype_parameters(PhysiCell::Phenotype& phenotype) {return 0;}
    int validate_PhysiCell_tokens(PhysiCell::Phenotype& phenotype) {return 0; }
    int validate_SBML_species() {return 0; }
    int create_custom_data_for_SBML(PhysiCell::Phenotype& phenotype) {return 0; }

};

#endif
