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


// Constructor that unpacks a cell
Moore_Cell_Info::Moore_Cell_Info(std::vector<char>& buffer, int size, int &pos)
{
    MPI_Unpack(buffer.data(), size, &pos, &ID, 1, MPI_INT, MPI_COMM_WORLD);
    position.resize(3);
    MPI_Unpack(buffer.data(), size, &pos, position.data(), 3, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer.data(), size, &pos, &radius, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer.data(), size, &pos, &nuclear_radius, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer.data(), size, &pos, &cell_cell_repulsion_strength, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer.data(), size, &pos, &relative_maximum_adhesion_distance, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer.data(), size, &pos, &cell_cell_adhesion_strength, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer.data(), size, &pos, &type, 1, MPI_INT, MPI_COMM_WORLD);
}

Moore_Voxel_Info::Moore_Voxel_Info(std::vector<char>& buffer, int size, int &pos)
{
    MPI_Unpack(buffer.data(), size, &pos, &global_mesh_index, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buffer.data(), size, &pos, &no_of_cells_in_vxl, 1, MPI_INT, MPI_COMM_WORLD);
    center.resize(3);
    MPI_Unpack(buffer.data(), size, &pos, center.data(), 3, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Unpack(buffer.data(), size, &pos, &max_cell_interactive_distance_in_voxel, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    moore_cells.reserve(no_of_cells_in_vxl);
    for (int i = 0; i < no_of_cells_in_vxl; ++i)
    {
        moore_cells.emplace_back(buffer, size, pos);
    }
}


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

void Cell_Container::evaluate_cell_elastic_interactions( PhysiCell::Cell *pCell , double dt_mec,  mpi_Environment &world, mpi_Cartesian &cart_topo){

    dynamic_spring_attachments(pCell,pCell->phenotype,dt_mec); 
    for( int j=0; j < pCell->state.spring_attachments.size(); j++ )
    {
        Cell* pC1 = pCell->state.spring_attachments[j]; 
        // standard_elastic_contact_function_confluent_rest_length(pC,pC->phenotype,pC1,pC1->phenotype,time_since_last_mechanics);  
        standard_elastic_contact_function(pCell,pCell->phenotype,pC1,pC1->phenotype,dt_mec);  
    }

    //elastic interactions with other subdomains
    //check if voxel is left or right before!

    int max_elastic = pCell->phenotype.mechanics.maximum_number_of_attachments;
    int done_elastic = pCell->state.spring_attachments.size();
    vector<Interacting_Cell_Info> neighbours = get_neighbour_interacting_cells(pCell, world, cart_topo);

    int i = 0;
    while (done_elastic < max_elastic && i < neighbours.size()) {
        standard_elastic_contact_function(pCell,pCell->phenotype, neighbours[i].position ,dt_mec, world, cart_topo);
        ++i;
        ++done_elastic;
    }
}

void Cell_Container::exchange_mechanics_halos(mpi_Environment &world, mpi_Cartesian &cart_topo) {
    // Exchange ghost cells for mechanical potentials and cell-cell interactions once per mechanics step.
    pack_moore_info(world, cart_topo);
    pack_cell_interact_info(world, cart_topo);
}

void Cell_Container::evaluate_cell_cell_interactions(double time_since_last_mechanics,  mpi_Environment &world, mpi_Cartesian &cart_topo) {
    // Boundary ghost data must be available prior to calling (exchange_mechanics_halos).

    #pragma omp parallel for 
    for( int i=0; i < (*all_cells).size(); i++ )
    {
        Cell* pC = (*all_cells)[i];

        standard_cell_cell_interactions( pC, pC->phenotype, time_since_last_mechanics , world, cart_topo );
    }

    unpack_cell_interact_info(world, cart_topo);

}

void Cell_Container::update_cell_potentials(double time_since_last_mechanics,  mpi_Environment &world, mpi_Cartesian &cart_topo){

    // Boundary data must be exchanged before calling this routine.
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

                evaluate_cell_elastic_interactions( pC , time_since_last_mechanics, world, cart_topo);
                
                
            }
        }
    
    }	
}
};
