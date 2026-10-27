#include "PhysiCell_standard_models.h" 
#include "PhysiCell_cell.h" 

namespace PhysiCell{

namespace
{
inline int current_death_model_code( const Death& death )
{
	const int index = death.current_death_model_index;
	if( index < 0 || index >= static_cast<int>(death.models.size()) || death.models[index] == NULL )
	{ return -1; }
	return death.models[index]->code;
}

inline double dead_phagocytosis_rate( const Cell_Interactions& interactions, int death_model_code )
{
	if( death_model_code == PhysiCell_constants::apoptosis_death_model )
	{ return interactions.apoptotic_phagocytosis_rate; }
	if( death_model_code == PhysiCell_constants::necrosis_death_model )
	{ return interactions.necrotic_phagocytosis_rate; }
	return interactions.other_dead_phagocytosis_rate;
}

inline void interact_with_neighboring_subdomain(
	Cell* pCell, Phenotype& phenotype, double dt, bool& attacked,
	const std::vector<int>& moore_list,
	const std::unordered_map<int, std::size_t>& ghost_voxel_indices,
	std::vector<Interacting_Voxel>& interacting_voxels,
	mpi_Environment& world, mpi_Cartesian& cart_topo )
{
	for( const int voxel_index : moore_list )
	{
		const auto voxel = ghost_voxel_indices.find( voxel_index );
		if( voxel == ghost_voxel_indices.end() )
		{ continue; }

		Interacting_Voxel& interacting_voxel = interacting_voxels[voxel->second];
		for( Interacting_Cell_Info& interacting_cell : interacting_voxel.cells )
		{
			if( std::isnan(interacting_cell.position[0]) )
			{
				std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;
				continue;
			}

			const int type = interacting_cell.type;
			if( interacting_cell.dead )
			{
				const double probability = dead_phagocytosis_rate(
					phenotype.cell_interactions, interacting_cell.death_model_code ) * dt;
				if( probability > 0.0 && UniformRandom() < probability )
				{ pCell->ingest_cell(&interacting_cell); }
				continue;
			}

			if( UniformRandom() <= phenotype.cell_interactions.live_phagocytosis_rates[type] * dt )
			{ pCell->ingest_cell(&interacting_cell); }

			if( UniformRandom() <= phenotype.cell_interactions.attack_rates[type] * dt && !attacked )
			{
				pCell->attack_cell(&interacting_cell,dt);
				attacked = true;
			}

			if( UniformRandom() <= phenotype.cell_interactions.fusion_rates[type] * dt )
			{ pCell->fuse_cell(&interacting_cell, world, cart_topo); }
		}
	}
}
}

void standard_elastic_contact_function( Cell* pC1, Phenotype& p1, std::vector<double> neighbour, double, mpi_Environment&, mpi_Cartesian& )
{
	if( pC1->position.size() != 3 || neighbour.size() != 3 )
	{ return; }

	std::vector<double> displacement = neighbour;
	displacement -= pC1->position;

	axpy( &(pC1->velocity) , p1.mechanics.attachment_elastic_constant , displacement );

	return;
}

//MPI parallel version
void standard_cell_cell_interactions( Cell* pCell, Phenotype& phenotype, double dt , mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	if( phenotype.death.dead == true )
	{ return; }

	Cell* pC = NULL; 
	int type = -1; 
	double probability = 0.0; 

	// std::cout << "testing against " << pCell->state.neighbors.size() << " cells " << std::endl; 

	bool attacked = false; 
	for( int n=0; n < pCell->state.neighbors.size(); n++ )
	{
		pC = pCell->state.neighbors[n]; 
		type = pC->type; 
		if (std::isnan(pC->position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;
		if( pC->phenotype.death.dead == true )
		{
			// dead phagocytosis 
			probability = dead_phagocytosis_rate(
				phenotype.cell_interactions, current_death_model_code(pC->phenotype.death) ) * dt;
			if( probability > 0.0 && UniformRandom() < probability ) 
			{ pCell->ingest_cell(pC); 
			  //std::cout << "		Rank " << world.rank << " Cell with ID " << pC->ID << " is ingested " << std::endl;
			  if (std::isnan(pCell->position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;
			} 
		}
		else
		{
			// live phagocytosis
			probability = phenotype.cell_interactions.live_phagocytosis_rates[type] * dt;
			if (std::isnan(pC->position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;
			if (std::isnan(pCell->position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;
			if( UniformRandom() <= probability ) 
			{ pCell->ingest_cell(pC); 
			  //std::cout << "		Rank " << world.rank << " Cell with ID " << pC->ID << " is ingested " << std::endl;
			 if (std::isnan(pCell->position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;
			} 

			// attack 

			// assume you can only attack one cell at a time 
			probability = phenotype.cell_interactions.attack_rates[type] * dt;  
			if( UniformRandom() <= probability && attacked == false ) 
			{
				pCell->attack_cell(pC,dt);
				if (std::isnan(pCell->position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl; 
				attacked = true;
				//std::cout << "		Rank " << world.rank << " Cell with ID " << pC->ID << " is attacked" << std::endl;
			} 

			// fusion 
			probability = phenotype.cell_interactions.fusion_rates[type] * dt;  
			if( UniformRandom() <= probability ) 
			{ pCell->fuse_cell(pC); 
				if (std::isnan(pCell->position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;
			  //std::cout << "		Rank " << world.rank << " Cell with ID " << pC->ID << " is fused " << std::endl;
			} 
		}

	

	}
//	std::cout << std::endl; 
	
	
	Cell_Container* container = pCell->get_container();
	int x_dim = container->underlying_mesh.x_coordinates.size();
	int y_dim = container->underlying_mesh.y_coordinates.size();
	int z_dim = container->underlying_mesh.z_coordinates.size();
	
	int local_vxl_inex = pCell->get_current_mechanics_voxel_index();
	//int position_voxel = (local_vxl_inex % (x_dim * y_dim)) % x_dim;  Layout warning!
	int position_voxel = (local_vxl_inex / (z_dim * y_dim));

	if (std::isnan(pCell->position[0])) std::cout << "Error in line " << __func__ << " " << __LINE__ << std::endl;

	if( position_voxel == 0 && world.rank > 0 )
	{
		const int yx_index = local_vxl_inex % (y_dim*z_dim);
		interact_with_neighboring_subdomain(
			pCell, phenotype, dt, attacked,
			container->underlying_mesh.moore_connected_voxel_global_indices_left[yx_index],
			container->ghost_voxel_from_left, container->ivfl, world, cart_topo );
	}
	if( position_voxel == x_dim-1 && world.rank < world.size-1 )
	{
		const int yx_index = local_vxl_inex % (y_dim*z_dim);
		interact_with_neighboring_subdomain(
			pCell, phenotype, dt, attacked,
			container->underlying_mesh.moore_connected_voxel_global_indices_right[yx_index],
			container->ghost_voxel_from_right, container->ivfr, world, cart_topo );
	}
}

};
