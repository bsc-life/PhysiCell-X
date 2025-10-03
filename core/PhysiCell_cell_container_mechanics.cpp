#include "../BioFVM/BioFVM_agent_container.h"
#include "PhysiCell_constants.h"
#include "../BioFVM/BioFVM_vector.h"
#include "PhysiCell_cell.h"
#include "PhysiCell_standard_models.h"
#include "PhysiCell_cell_container.h"


#include <algorithm>
#include <iterator> 

using namespace BioFVM;


namespace PhysiCell{

std::vector<Interacting_Cell_Info> Cell_Container::get_neighbour_interacting_cells(Cell* pCell, mpi_Environment &world, mpi_Cartesian &cart_topo){
    std::vector<Interacting_Cell_Info> neighbours;
    int x_dim = pCell->get_container()->underlying_mesh.x_coordinates.size(); 
	int y_dim = pCell->get_container()->underlying_mesh.y_coordinates.size();
	int z_dim = pCell->get_container()->underlying_mesh.z_coordinates.size();
	
	int local_vxl_inex = pCell->get_current_mechanics_voxel_index();
	int position_voxel = (local_vxl_inex / (z_dim * y_dim));

    if (position_voxel == 0  && world.rank > 0) { //left edge
        int yx_index = local_vxl_inex %(y_dim*z_dim);
		std::vector<int> moore_list = pCell->get_container()->underlying_mesh.moore_connected_voxel_global_indices_left[yx_index];
        for(int i=0; i<moore_list.size(); i++)
        {
            Interacting_Voxel &ivi = pCell->get_container()->um_ivfl.at(moore_list[i]);
                
            for(int cell_ctr=0; cell_ctr<ivi.cells.size(); cell_ctr++) {
                //pCell->add_potentials(mvi.cells[cell_ctr], world, cart_topo); 
                neighbours.push_back(ivi.cells[cell_ctr]);
            }
        } 
    } 
    else if (position_voxel == (x_dim-1) && (world.rank < world.size-1)) { //right edge
        int yx_index = local_vxl_inex %(y_dim*z_dim);
		std::vector<int> moore_list = pCell->get_container()->underlying_mesh.moore_connected_voxel_global_indices_right[yx_index];
        for(int i=0; i<moore_list.size(); i++)
        {
            Interacting_Voxel &ivi = pCell->get_container()->um_ivfr.at(moore_list[i]);
                
            for(int cell_ctr=0; cell_ctr<ivi.cells.size(); cell_ctr++) {
                //pCell->add_potentials(mvi.cells[cell_ctr], world, cart_topo); 
                neighbours.push_back(ivi.cells[cell_ctr]);
            }
        } 
    }
    return neighbours;
}

void Cell_Container::evaluate_cell_cell_interactions(double time_since_last_mechanics,  mpi_Environment &world, mpi_Cartesian &cart_topo) {

   
    pack_cell_interact_info(world, cart_topo);

    #pragma omp parallel for 
    for( int i=0; i < (*all_cells).size(); i++ )
    {
        Cell* pC = (*all_cells)[i];

        standard_cell_cell_interactions( pC, pC->phenotype, time_since_last_mechanics , world, cart_topo );
    }

    unpack_cell_interact_info(world, cart_topo);

}

void Cell_Container::update_cell_potentials(double time_since_last_mechanics,  mpi_Environment &world, mpi_Cartesian &cart_topo){

    pack_moore_info(world, cart_topo);
    //pack your subdomains and send
    //pack_cell_interact_info(world, cart_topo);
    //Main loop to traverse all the cells
    #pragma omp parallel for 
    for( int i=0; i < (*all_cells).size(); i++ )
    {
        Cell* pC = (*all_cells)[i];
        
        if( pC->functions.contact_function && pC->is_out_of_domain == false )
			{ evaluate_interactions( pC,pC->phenotype,time_since_last_mechanics ); } //No parallelism necessary yet


        if( pC->functions.custom_cell_rule && pC->is_out_of_domain == false )
			{ pC->functions.custom_cell_rule( pC,pC->phenotype,time_since_last_mechanics ); } //No parallelism necessary yet

        if( PhysiCell_settings.disable_automated_spring_adhesions == false ) 
		{
            if( pC->is_movable )
            {
                if(pC->functions.update_velocity_parallel && pC->is_out_of_domain == false && pC->is_movable)
                {
                    pC->functions.update_velocity_parallel( pC, pC->phenotype, time_since_last_mechanics, world, cart_topo);
                }

                dynamic_spring_attachments(pC,pC->phenotype,time_since_last_mechanics); 
                for( int j=0; j < pC->state.spring_attachments.size(); j++ )
                {
                    Cell* pC1 = pC->state.spring_attachments[j]; 
                    // standard_elastic_contact_function_confluent_rest_length(pC,pC->phenotype,pC1,pC1->phenotype,time_since_last_mechanics);  
                    standard_elastic_contact_function(pC,pC->phenotype,pC1,pC1->phenotype,time_since_last_mechanics);  
                }

                if (std::isnan(pC->position[0])) {
                std::cout << "I fucked here 3 | " << pC->ID << std::endl;
                continue;
                }
                int max_elastic = pC->phenotype.mechanics.maximum_number_of_attachments;
                int done_elastic = pC->state.spring_attachments.size();
                /*
                vector<Interacting_Cell_Info> neighbours = get_neighbour_interacting_cells(pC, world, cart_topo);
                int i = 0;
                while (done_elastic < max_elastic && i < neighbours.size()) {
                    standard_elastic_contact_function(pC,pC->phenotype, neighbours[i].position ,time_since_last_mechanics, world, cart_topo);
                    ++i;
                    ++done_elastic;
                    }*/
                }
        }
    
    }	

}
};