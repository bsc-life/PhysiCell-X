#include "../BioFVM/BioFVM_agent_container.h"
#include "PhysiCell_constants.h"
#include "../BioFVM/BioFVM_vector.h"
#include "PhysiCell_cell.h"
#include "PhysiCell_standard_models.h"
#include "PhysiCell_cell_container.h"


#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstring>
#include <iterator> 
#include <stdexcept>

using namespace BioFVM;


namespace PhysiCell{

namespace
{

template <typename T>
inline void ghost_write(char*& cursor, const T& value)
{
    std::memcpy(cursor, &value, sizeof(T));
    cursor += sizeof(T);
}

template <typename T>
inline void ghost_read(const char*& cursor, T& value)
{
    std::memcpy(&value, cursor, sizeof(T));
    cursor += sizeof(T);
}

inline void ghost_write_doubles(char*& cursor, const double* values, std::size_t count)
{
    if (count > 0)
    {
        std::memcpy(cursor, values, count * sizeof(double));
        cursor += count * sizeof(double);
    }
}

inline void ghost_read_doubles(const char*& cursor, double* values, std::size_t count)
{
    if (count > 0)
    {
        std::memcpy(values, cursor, count * sizeof(double));
        cursor += count * sizeof(double);
    }
}

inline std::size_t ghost_voxel_header_bytes()
{
    return 2 * sizeof(std::int32_t) + 4 * sizeof(double);
}

inline std::size_t ghost_cell_bytes(int substrates, bool interactions_active)
{
    // Shared mechanics fields: ID, type, position and five mechanics scalars.
    std::size_t bytes = 2 * sizeof(std::int32_t) + 8 * sizeof(double);
    if (interactions_active)
    {
        // dead, number_of_nuclei, seven volume fields and two substrate vectors.
        bytes += sizeof(std::uint8_t) + sizeof(std::int32_t);
        bytes += (7 + 2 * static_cast<std::size_t>(substrates)) * sizeof(double);
    }
    return bytes;
}

inline int checked_mpi_count(std::size_t bytes)
{
    if (bytes > static_cast<std::size_t>(INT_MAX))
    {
        throw std::runtime_error("ghost exchange payload exceeds the MPI int count limit");
    }
    return static_cast<int>(bytes);
}

bool interaction_rates_active(const Cell_Interactions& interactions)
{
    if (interactions.apoptotic_phagocytosis_rate > 0.0 ||
        interactions.necrotic_phagocytosis_rate > 0.0 ||
        interactions.other_dead_phagocytosis_rate > 0.0)
    {
        return true;
    }

    for (std::size_t i = 0; i < interactions.live_phagocytosis_rates.size(); ++i)
    {
        if (interactions.live_phagocytosis_rates[i] > 0.0)
        { return true; }
    }
    for (std::size_t i = 0; i < interactions.attack_rates.size(); ++i)
    {
        if (interactions.attack_rates[i] > 0.0)
        { return true; }
    }
    for (std::size_t i = 0; i < interactions.fusion_rates.size(); ++i)
    {
        if (interactions.fusion_rates[i] > 0.0)
        { return true; }
    }
    return false;
}

void serialize_ghost_voxel(Cell_Container& container, int voxel_index,
    int substrates, bool interactions_active, char* destination)
{
    char* cursor = destination;
    const BioFVM::Voxel& voxel = container.underlying_mesh.voxels[voxel_index];
    const std::vector<Cell*>& cells = container.agent_grid[voxel_index];
    const std::int32_t global_index = static_cast<std::int32_t>(voxel.global_mesh_index);
    const std::int32_t cell_count = static_cast<std::int32_t>(cells.size());

    ghost_write(cursor, global_index);
    ghost_write(cursor, cell_count);
    ghost_write_doubles(cursor, voxel.center.data(), 3);
    ghost_write(cursor, container.max_cell_interactive_distance_in_voxel[voxel_index]);

    for (std::size_t i = 0; i < cells.size(); ++i)
    {
        const Cell* cell = cells[i];
        const std::int32_t id = static_cast<std::int32_t>(cell->ID);
        const std::int32_t type = static_cast<std::int32_t>(cell->type);
        ghost_write(cursor, id);
        ghost_write(cursor, type);
        ghost_write_doubles(cursor, cell->position.data(), 3);
        ghost_write(cursor, cell->phenotype.geometry.radius);
        ghost_write(cursor, cell->phenotype.geometry.nuclear_radius);
        ghost_write(cursor, cell->phenotype.mechanics.cell_cell_repulsion_strength);
        ghost_write(cursor, cell->phenotype.mechanics.relative_maximum_adhesion_distance);
        ghost_write(cursor, cell->phenotype.mechanics.cell_cell_adhesion_strength);

        if (interactions_active)
        {
            const std::uint8_t dead = cell->phenotype.death.dead ? 1u : 0u;
            const std::int32_t nuclei = static_cast<std::int32_t>(cell->state.number_of_nuclei);
            ghost_write(cursor, dead);
            ghost_write(cursor, nuclei);
            ghost_write(cursor, cell->phenotype.volume.total);
            ghost_write(cursor, cell->phenotype.volume.cytoplasmic_fluid);
            ghost_write(cursor, cell->phenotype.volume.nuclear_fluid);
            ghost_write(cursor, cell->phenotype.volume.cytoplasmic_solid);
            ghost_write(cursor, cell->phenotype.volume.nuclear_solid);
            ghost_write(cursor, cell->phenotype.volume.target_solid_cytoplasmic);
            ghost_write(cursor, cell->phenotype.volume.target_solid_nuclear);
            ghost_write_doubles(cursor, cell->internalized_substrates->data(), substrates);
            ghost_write_doubles(cursor,
                cell->phenotype.molecular.fraction_transferred_when_ingested.data(), substrates);
        }
    }
}

void scan_ghost_offsets(const std::vector<char>& buffer, int voxel_count,
    std::size_t cell_bytes, std::vector<std::size_t>& offsets)
{
    offsets.resize(static_cast<std::size_t>(voxel_count) + 1);
    std::size_t cursor = 0;
    for (int voxel = 0; voxel < voxel_count; ++voxel)
    {
        if (cursor + ghost_voxel_header_bytes() > buffer.size())
        { throw std::runtime_error("truncated ghost voxel header"); }

        offsets[voxel] = cursor;
        const char* header = buffer.data() + cursor;
        std::int32_t global_index = 0;
        std::int32_t cell_count = 0;
        ghost_read(header, global_index);
        ghost_read(header, cell_count);
        (void) global_index;
        if (cell_count < 0)
        { throw std::runtime_error("negative ghost-cell count"); }

        const std::size_t record_bytes = ghost_voxel_header_bytes() +
            static_cast<std::size_t>(cell_count) * cell_bytes;
        if (record_bytes > buffer.size() - cursor)
        { throw std::runtime_error("truncated ghost-cell payload"); }
        cursor += record_bytes;
    }
    offsets[voxel_count] = cursor;
    if (cursor != buffer.size())
    { throw std::runtime_error("ghost payload has trailing bytes"); }
}

void deserialize_ghost_voxel(const std::vector<char>& buffer, std::size_t offset,
    int substrates, bool interactions_active, Interacting_Voxel& voxel)
{
    const char* cursor = buffer.data() + offset;
    std::int32_t global_index = 0;
    std::int32_t cell_count = 0;
    ghost_read(cursor, global_index);
    ghost_read(cursor, cell_count);
    voxel.global_mesh_index = global_index;
    voxel.no_of_cells_in_vxl = cell_count;
    voxel.center.resize(3);
    ghost_read_doubles(cursor, voxel.center.data(), 3);
    ghost_read(cursor, voxel.max_cell_interactive_distance_in_voxel);
    voxel.cells.resize(static_cast<std::size_t>(cell_count));

    for (int i = 0; i < cell_count; ++i)
    {
        Interacting_Cell_Info& cell = voxel.cells[i];
        std::int32_t id = 0;
        std::int32_t type = 0;
        ghost_read(cursor, id);
        ghost_read(cursor, type);
        cell.ID = id;
        cell.type = type;
        cell.position.resize(3);
        ghost_read_doubles(cursor, cell.position.data(), 3);
        ghost_read(cursor, cell.radius);
        ghost_read(cursor, cell.nuclear_radius);
        ghost_read(cursor, cell.cell_cell_repulsion_strength);
        ghost_read(cursor, cell.relative_maximum_adhesion_distance);
        ghost_read(cursor, cell.cell_cell_adhesion_strength);

        cell.attacked = false;
        cell.damage_suffered = 0.0;
        cell.time_attacked = 0.0;
        cell.fused = false;
        cell.ingested = false;

        if (interactions_active)
        {
            std::uint8_t dead = 0;
            std::int32_t nuclei = 0;
            ghost_read(cursor, dead);
            ghost_read(cursor, nuclei);
            cell.dead = (dead != 0);
            cell.number_of_nuclei = nuclei;
            ghost_read(cursor, cell.phenotype_volume);
            ghost_read(cursor, cell.cytoplasmic_fluid);
            ghost_read(cursor, cell.nuclear_fluid);
            ghost_read(cursor, cell.cytoplasmic_solid);
            ghost_read(cursor, cell.nuclear_solid);
            ghost_read(cursor, cell.target_solid_cytoplasmic);
            ghost_read(cursor, cell.target_solid_nuclear);
            cell.internalized_substrates.resize(substrates);
            ghost_read_doubles(cursor, cell.internalized_substrates.data(), substrates);
            cell.fraction_transferred_when_ingested.resize(substrates);
            ghost_read_doubles(cursor, cell.fraction_transferred_when_ingested.data(), substrates);
        }
    }
}

inline std::uint8_t ghost_event_mask(const Interacting_Cell_Info& cell)
{
    return static_cast<std::uint8_t>((cell.attacked ? 1u : 0u) |
        (cell.fused ? 2u : 0u) | (cell.ingested ? 4u : 0u));
}

inline std::size_t ghost_event_bytes()
{
    return sizeof(std::int32_t) + sizeof(std::uint8_t) + 2 * sizeof(double);
}

void serialize_ghost_events(const std::vector<Interacting_Voxel>& voxels,
    std::vector<char>& buffer)
{
    std::size_t event_count = 0;
    for (std::size_t v = 0; v < voxels.size(); ++v)
    {
        for (std::size_t c = 0; c < voxels[v].cells.size(); ++c)
        {
            if (ghost_event_mask(voxels[v].cells[c]) != 0)
            { ++event_count; }
        }
    }

    buffer.resize(event_count * ghost_event_bytes());
    char* cursor = buffer.empty() ? NULL : buffer.data();
    for (std::size_t v = 0; v < voxels.size(); ++v)
    {
        for (std::size_t c = 0; c < voxels[v].cells.size(); ++c)
        {
            const Interacting_Cell_Info& cell = voxels[v].cells[c];
            const std::uint8_t mask = ghost_event_mask(cell);
            if (mask == 0)
            { continue; }
            const std::int32_t id = static_cast<std::int32_t>(cell.ID);
            ghost_write(cursor, id);
            ghost_write(cursor, mask);
            ghost_write(cursor, cell.damage_suffered);
            ghost_write(cursor, cell.time_attacked);
        }
    }
}

void apply_ghost_events(const std::vector<char>& buffer,
    const std::unordered_map<int, Cell*>& owned_cells)
{
    if (buffer.size() % ghost_event_bytes() != 0)
    { throw std::runtime_error("invalid sparse ghost-event payload"); }

    const char* cursor = buffer.empty() ? NULL : buffer.data();
    const std::size_t count = buffer.size() / ghost_event_bytes();
    for (std::size_t i = 0; i < count; ++i)
    {
        std::int32_t id = 0;
        std::uint8_t mask = 0;
        double damage = 0.0;
        double attack_time = 0.0;
        ghost_read(cursor, id);
        ghost_read(cursor, mask);
        ghost_read(cursor, damage);
        ghost_read(cursor, attack_time);

        std::unordered_map<int, Cell*>::const_iterator found = owned_cells.find(id);
        if (found == owned_cells.end())
        { continue; }
        Cell* cell = found->second;
        if ((mask & 1u) != 0)
        { cell->was_attacked(damage, attack_time); }
        if ((mask & 2u) != 0)
        { cell->was_fused(); }
        if ((mask & 4u) != 0)
        { cell->was_ingested(); }
    }
}

} // anonymous namespace


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

//TODO: 1: encapsulate the packing of each voxel's data into a separate function to avoid code repetition
std::vector<Moore_Cell_Info> Cell_Container::get_moore_neighbour_cells(Cell* pCell, mpi_Environment &world, mpi_Cartesian &cart_topo){
    std::vector<Moore_Cell_Info> neighbours;
    int x_dim = pCell->get_container()->underlying_mesh.x_coordinates.size(); 
    int y_dim = pCell->get_container()->underlying_mesh.y_coordinates.size();
    int z_dim = pCell->get_container()->underlying_mesh.z_coordinates.size();
    
    int local_vxl_index = pCell->get_current_mechanics_voxel_index();
    int position_voxel = (local_vxl_index / (z_dim * y_dim));

    if (position_voxel == 0  && world.rank > 0) { //left edge
        int yx_index = local_vxl_index %(y_dim*z_dim);
        std::vector<int> moore_list = pCell->get_container()->underlying_mesh.moore_connected_voxel_global_indices_left[yx_index];
        for(int i=0; i<moore_list.size(); i++)
        {
            auto mvi_iter = ghost_voxel_from_left.find(moore_list[i]);
            if( mvi_iter == ghost_voxel_from_left.end() )
            { continue; }
            Interacting_Voxel &mvi = ivfl[mvi_iter->second];
            for(int cell_ctr=0; cell_ctr<mvi.cells.size(); cell_ctr++) {
                neighbours.push_back(static_cast<Moore_Cell_Info&>(mvi.cells[cell_ctr]));
            }
        } 
    } 
    if (position_voxel == (x_dim-1) && (world.rank < world.size-1)) { //right edge
        int yx_index = local_vxl_index %(y_dim*z_dim);
        std::vector<int> moore_list = pCell->get_container()->underlying_mesh.moore_connected_voxel_global_indices_right[yx_index];
        for(int i=0; i<moore_list.size(); i++)
        {
            auto mvi_iter = ghost_voxel_from_right.find(moore_list[i]);
            if( mvi_iter == ghost_voxel_from_right.end() )
            { continue; }
            Interacting_Voxel &mvi = ivfr[mvi_iter->second];
                
            for(int cell_ctr=0; cell_ctr<mvi.cells.size(); cell_ctr++) {
                neighbours.push_back(static_cast<Moore_Cell_Info&>(mvi.cells[cell_ctr]));
            }
        } 
    }
    return neighbours;
}

//TODO: 1: encapsulate the packing of each voxel's data into a separate function to avoid code repetition
std::vector<Interacting_Cell_Info> Cell_Container::get_neighbour_interacting_cells(Cell* pCell, mpi_Environment &world, mpi_Cartesian &cart_topo){
    std::vector<Interacting_Cell_Info> neighbours;
    int x_dim = pCell->get_container()->underlying_mesh.x_coordinates.size(); 
    int y_dim = pCell->get_container()->underlying_mesh.y_coordinates.size();
    int z_dim = pCell->get_container()->underlying_mesh.z_coordinates.size();
    
    int local_vxl_index = pCell->get_current_mechanics_voxel_index();
    int position_voxel = (local_vxl_index / (z_dim * y_dim));

    if (position_voxel == 0  && world.rank > 0) { //left edge
        int yx_index = local_vxl_index %(y_dim*z_dim);
        std::vector<int> moore_list = pCell->get_container()->underlying_mesh.moore_connected_voxel_global_indices_left[yx_index];
        for(int i=0; i<moore_list.size(); i++)
        {
            auto ivi_iter = pCell->get_container()->ghost_voxel_from_left.find(moore_list[i]);
            if( ivi_iter == pCell->get_container()->ghost_voxel_from_left.end() )
            { continue; }
            Interacting_Voxel &ivi = pCell->get_container()->ivfl[ivi_iter->second];
                
            for(int cell_ctr=0; cell_ctr<ivi.cells.size(); cell_ctr++) {
                neighbours.push_back(ivi.cells[cell_ctr]);
            }
        } 
    } 
    else if (position_voxel == (x_dim-1) && (world.rank < world.size-1)) { //right edge
        int yx_index = local_vxl_index %(y_dim*z_dim);
        std::vector<int> moore_list = pCell->get_container()->underlying_mesh.moore_connected_voxel_global_indices_right[yx_index];
        for(int i=0; i<moore_list.size(); i++)
        {
            auto ivi_iter = pCell->get_container()->ghost_voxel_from_right.find(moore_list[i]);
            if( ivi_iter == pCell->get_container()->ghost_voxel_from_right.end() )
            { continue; }
            Interacting_Voxel &ivi = pCell->get_container()->ivfr[ivi_iter->second];
                
            for(int cell_ctr=0; cell_ctr<ivi.cells.size(); cell_ctr++) {
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
		standard_elastic_contact_function(pCell,pCell->phenotype,pC1,pC1->phenotype,dt_mec);  
    }

    //elastic interactions with other subdomains
    //check if voxel is left or right before!

    int max_elastic = pCell->phenotype.mechanics.maximum_number_of_attachments;
    int done_elastic = pCell->state.spring_attachments.size();
    vector<Moore_Cell_Info> neighbours = get_moore_neighbour_cells(pCell, world, cart_topo);
    
     double attachment_probability = pCell->phenotype.mechanics.attachment_rate * dt_mec; 
    int i = 0;
    while (done_elastic < max_elastic && i < neighbours.size()) {
        double probability_of_attachment = attachment_probability * pCell->phenotype.mechanics.cell_adhesion_affinities[neighbours[i].type];
        if (uniform_random() < probability_of_attachment) {
            standard_elastic_contact_function(pCell,pCell->phenotype, neighbours[i].position ,dt_mec, world, cart_topo);
			++done_elastic;
		}
        ++i;
        
    }
}

void Cell_Container::invalidate_ghost_interaction_activity()
{
    ghost_interaction_activity_initialized = false;
}

bool Cell_Container::ghost_cell_interactions_enabled() const
{
    return ghost_interactions_active;
}

void Cell_Container::ghost_exchange(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
    // Rules can change interaction rates per cell, so dynamic rules require a
    // fresh local scan. Static configurations are scanned once and cached.
    if (!ghost_interaction_activity_initialized || PhysiCell_settings.rules_enabled)
    {
        bool local_active = false;
        for (std::size_t i = 0; i < cell_definitions_by_index.size() && !local_active; ++i)
        {
            local_active = interaction_rates_active(
                cell_definitions_by_index[i]->phenotype.cell_interactions);
        }

        // The first scan also catches setup-time, cell-specific overrides that
        // differ from the cell definitions. Rules-enabled models repeat it.
        if (!ghost_interaction_activity_initialized || PhysiCell_settings.rules_enabled)
        {
            for (std::size_t i = 0; i < all_cells->size() && !local_active; ++i)
            {
                local_active = interaction_rates_active(
                    (*all_cells)[i]->phenotype.cell_interactions);
            }
        }

        int local_flag = local_active ? 1 : 0;
        int global_flag = 0;
        MPI_Allreduce(&local_flag, &global_flag, 1, MPI_INT, MPI_MAX,
            cart_topo.mpi_cart_comm);
        ghost_interactions_active = (global_flag != 0);
        ghost_interaction_activity_initialized = true;
    }

    const int z_dim = static_cast<int>(underlying_mesh.z_coordinates.size());
    const int y_dim = static_cast<int>(underlying_mesh.y_coordinates.size());
    const int x_dim = static_cast<int>(underlying_mesh.x_coordinates.size());
    const int voxel_count = y_dim * z_dim;
    const int substrates = underlying_mesh.n_substrates;
    const std::size_t cell_bytes = ghost_cell_bytes(substrates, ghost_interactions_active);
    const std::size_t header_bytes = ghost_voxel_header_bytes();

    ghost_offsets_left.resize(static_cast<std::size_t>(voxel_count) + 1);
    ghost_offsets_right.resize(static_cast<std::size_t>(voxel_count) + 1);
    std::size_t total_left = 0;
    std::size_t total_right = 0;
    std::size_t boundary_cell_count = 0;

    for (int j = 0; j < y_dim; ++j)
    {
        for (int k = 0; k < z_dim; ++k)
        {
            const int plane_index = j * z_dim + k;
            ghost_offsets_left[plane_index] = total_left;
            ghost_offsets_right[plane_index] = total_right;

            if (world.rank > 0)
            {
                const int voxel = underlying_mesh.voxel_index(0, j, k);
                const std::size_t cells = agent_grid[voxel].size();
                total_left += header_bytes + cells * cell_bytes;
                boundary_cell_count += cells;
            }
            if (world.rank < world.size - 1)
            {
                const int voxel = underlying_mesh.voxel_index(x_dim - 1, j, k);
                const std::size_t cells = agent_grid[voxel].size();
                total_right += header_bytes + cells * cell_bytes;
                boundary_cell_count += cells;
            }
        }
    }
    ghost_offsets_left[voxel_count] = total_left;
    ghost_offsets_right[voxel_count] = total_right;

    // Sparse replies identify owning cells by ID. Build this only when the
    // interaction protocol is active; duplicate insertion is harmless when a
    // one-voxel-wide subdomain exposes the same cell on both faces.
    owned_boundary_cells_by_id.clear();
    if (ghost_interactions_active)
    {
        owned_boundary_cells_by_id.reserve(boundary_cell_count);
        for (int j = 0; j < y_dim; ++j)
        {
            for (int k = 0; k < z_dim; ++k)
            {
                if (world.rank > 0)
                {
                    const int voxel = underlying_mesh.voxel_index(0, j, k);
                    for (std::size_t i = 0; i < agent_grid[voxel].size(); ++i)
                    {
                        Cell* cell = agent_grid[voxel][i];
                        owned_boundary_cells_by_id[cell->ID] = cell;
                    }
                }
                if (world.rank < world.size - 1)
                {
                    const int voxel = underlying_mesh.voxel_index(x_dim - 1, j, k);
                    for (std::size_t i = 0; i < agent_grid[voxel].size(); ++i)
                    {
                        Cell* cell = agent_grid[voxel][i];
                        owned_boundary_cells_by_id[cell->ID] = cell;
                    }
                }
            }
        }
    }

    ghost_snd_left.resize(total_left);
    ghost_snd_right.resize(total_right);

    if (world.rank > 0)
    {
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < y_dim; ++j)
        {
            for (int k = 0; k < z_dim; ++k)
            {
                const int plane_index = j * z_dim + k;
                const int voxel = underlying_mesh.voxel_index(0, j, k);
                serialize_ghost_voxel(*this, voxel, substrates,
                    ghost_interactions_active,
                    ghost_snd_left.data() + ghost_offsets_left[plane_index]);
            }
        }
    }
    if (world.rank < world.size - 1)
    {
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < y_dim; ++j)
        {
            for (int k = 0; k < z_dim; ++k)
            {
                const int plane_index = j * z_dim + k;
                const int voxel = underlying_mesh.voxel_index(x_dim - 1, j, k);
                serialize_ghost_voxel(*this, voxel, substrates,
                    ghost_interactions_active,
                    ghost_snd_right.data() + ghost_offsets_right[plane_index]);
            }
        }
    }

    int send_left_bytes = checked_mpi_count(total_left);
    int send_right_bytes = checked_mpi_count(total_right);
    int receive_right_bytes = 0;
    int receive_left_bytes = 0;
    MPI_Request size_requests[4];
    MPI_Irecv(&receive_right_bytes, 1, MPI_INT, cart_topo.X_RIGHT, 7101,
        cart_topo.mpi_cart_comm, &size_requests[0]);
    MPI_Irecv(&receive_left_bytes, 1, MPI_INT, cart_topo.X_LEFT, 7102,
        cart_topo.mpi_cart_comm, &size_requests[1]);
    MPI_Isend(&send_left_bytes, 1, MPI_INT, cart_topo.X_LEFT, 7101,
        cart_topo.mpi_cart_comm, &size_requests[2]);
    MPI_Isend(&send_right_bytes, 1, MPI_INT, cart_topo.X_RIGHT, 7102,
        cart_topo.mpi_cart_comm, &size_requests[3]);
    MPI_Waitall(4, size_requests, MPI_STATUSES_IGNORE);

    if (receive_right_bytes < 0 || receive_left_bytes < 0)
    { throw std::runtime_error("negative ghost exchange receive size"); }
    ghost_rcv_right.resize(static_cast<std::size_t>(receive_right_bytes));
    ghost_rcv_left.resize(static_cast<std::size_t>(receive_left_bytes));

    MPI_Request payload_requests[4];
    MPI_Irecv(ghost_rcv_right.empty() ? NULL : ghost_rcv_right.data(),
        receive_right_bytes, MPI_BYTE, cart_topo.X_RIGHT, 7103,
        cart_topo.mpi_cart_comm, &payload_requests[0]);
    MPI_Irecv(ghost_rcv_left.empty() ? NULL : ghost_rcv_left.data(),
        receive_left_bytes, MPI_BYTE, cart_topo.X_LEFT, 7104,
        cart_topo.mpi_cart_comm, &payload_requests[1]);
    MPI_Isend(ghost_snd_left.empty() ? NULL : ghost_snd_left.data(),
        send_left_bytes, MPI_BYTE, cart_topo.X_LEFT, 7103,
        cart_topo.mpi_cart_comm, &payload_requests[2]);
    MPI_Isend(ghost_snd_right.empty() ? NULL : ghost_snd_right.data(),
        send_right_bytes, MPI_BYTE, cart_topo.X_RIGHT, 7104,
        cart_topo.mpi_cart_comm, &payload_requests[3]);
    MPI_Waitall(4, payload_requests, MPI_STATUSES_IGNORE);

    if (world.rank < world.size - 1)
    {
        scan_ghost_offsets(ghost_rcv_right, voxel_count, cell_bytes,
            ghost_offsets_from_right);
        ivfr.resize(voxel_count);
        #pragma omp parallel for
        for (int voxel = 0; voxel < voxel_count; ++voxel)
        {
            deserialize_ghost_voxel(ghost_rcv_right,
                ghost_offsets_from_right[voxel], substrates,
                ghost_interactions_active, ivfr[voxel]);
        }

        if (!ghost_right_index_initialized)
        {
            ghost_voxel_from_right.clear();
            ghost_voxel_from_right.reserve(voxel_count);
            for (int voxel = 0; voxel < voxel_count; ++voxel)
            {
                ghost_voxel_from_right[ivfr[voxel].global_mesh_index] = voxel;
            }
            ghost_right_index_initialized = true;
        }
    }

    if (world.rank > 0)
    {
        scan_ghost_offsets(ghost_rcv_left, voxel_count, cell_bytes,
            ghost_offsets_from_left);
        ivfl.resize(voxel_count);
        #pragma omp parallel for
        for (int voxel = 0; voxel < voxel_count; ++voxel)
        {
            deserialize_ghost_voxel(ghost_rcv_left,
                ghost_offsets_from_left[voxel], substrates,
                ghost_interactions_active, ivfl[voxel]);
        }

        if (!ghost_left_index_initialized)
        {
            ghost_voxel_from_left.clear();
            ghost_voxel_from_left.reserve(voxel_count);
            for (int voxel = 0; voxel < voxel_count; ++voxel)
            {
                ghost_voxel_from_left[ivfl[voxel].global_mesh_index] = voxel;
            }
            ghost_left_index_initialized = true;
        }
    }
}

void Cell_Container::return_ghost_interaction_events(
    mpi_Environment &world, mpi_Cartesian &cart_topo)
{
    if (!ghost_interactions_active)
    { return; }

    if (world.rank > 0)
    { serialize_ghost_events(ivfl, ghost_event_snd_left); }
    else
    { ghost_event_snd_left.clear(); }

    if (world.rank < world.size - 1)
    { serialize_ghost_events(ivfr, ghost_event_snd_right); }
    else
    { ghost_event_snd_right.clear(); }

    int send_left_bytes = checked_mpi_count(ghost_event_snd_left.size());
    int send_right_bytes = checked_mpi_count(ghost_event_snd_right.size());
    int receive_right_bytes = 0;
    int receive_left_bytes = 0;
    MPI_Request size_requests[4];
    MPI_Irecv(&receive_right_bytes, 1, MPI_INT, cart_topo.X_RIGHT, 7201,
        cart_topo.mpi_cart_comm, &size_requests[0]);
    MPI_Irecv(&receive_left_bytes, 1, MPI_INT, cart_topo.X_LEFT, 7202,
        cart_topo.mpi_cart_comm, &size_requests[1]);
    MPI_Isend(&send_left_bytes, 1, MPI_INT, cart_topo.X_LEFT, 7201,
        cart_topo.mpi_cart_comm, &size_requests[2]);
    MPI_Isend(&send_right_bytes, 1, MPI_INT, cart_topo.X_RIGHT, 7202,
        cart_topo.mpi_cart_comm, &size_requests[3]);
    MPI_Waitall(4, size_requests, MPI_STATUSES_IGNORE);

    if (receive_right_bytes < 0 || receive_left_bytes < 0)
    { throw std::runtime_error("negative ghost-event receive size"); }
    ghost_event_rcv_right.resize(static_cast<std::size_t>(receive_right_bytes));
    ghost_event_rcv_left.resize(static_cast<std::size_t>(receive_left_bytes));

    MPI_Request payload_requests[4];
    MPI_Irecv(ghost_event_rcv_right.empty() ? NULL : ghost_event_rcv_right.data(),
        receive_right_bytes, MPI_BYTE, cart_topo.X_RIGHT, 7203,
        cart_topo.mpi_cart_comm, &payload_requests[0]);
    MPI_Irecv(ghost_event_rcv_left.empty() ? NULL : ghost_event_rcv_left.data(),
        receive_left_bytes, MPI_BYTE, cart_topo.X_LEFT, 7204,
        cart_topo.mpi_cart_comm, &payload_requests[1]);
    MPI_Isend(ghost_event_snd_left.empty() ? NULL : ghost_event_snd_left.data(),
        send_left_bytes, MPI_BYTE, cart_topo.X_LEFT, 7203,
        cart_topo.mpi_cart_comm, &payload_requests[2]);
    MPI_Isend(ghost_event_snd_right.empty() ? NULL : ghost_event_snd_right.data(),
        send_right_bytes, MPI_BYTE, cart_topo.X_RIGHT, 7204,
        cart_topo.mpi_cart_comm, &payload_requests[3]);
    MPI_Waitall(4, payload_requests, MPI_STATUSES_IGNORE);

    apply_ghost_events(ghost_event_rcv_right, owned_boundary_cells_by_id);
    apply_ghost_events(ghost_event_rcv_left, owned_boundary_cells_by_id);
}

void Cell_Container::exchange_mechanics_halos(mpi_Environment &world, mpi_Cartesian &cart_topo) {
    // One merged exchange supplies both mechanical and cell-interaction ghosts.
    auto start = std::chrono::high_resolution_clock::now();
    ghost_exchange(world, cart_topo);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    //std::cout << "Time taken to ghost exchange: " << elapsed.count() << " seconds" << std::endl;
}

void Cell_Container::evaluate_cell_cell_interactions(double time_since_last_mechanics,  mpi_Environment &world, mpi_Cartesian &cart_topo) {
    // Boundary ghost data must be available prior to calling (exchange_mechanics_halos).
	if (!ghost_interactions_active)
	{ return; }

    #pragma omp parallel for 
    for( int i=0; i < (*all_cells).size(); i++ )
    {
        Cell* pC = (*all_cells)[i];
		if (!pC->is_out_of_domain)
        	standard_cell_cell_interactions( pC, pC->phenotype, time_since_last_mechanics , world, cart_topo );
    }

	return_ghost_interaction_events(world, cart_topo);

}

void Cell_Container::update_cell_potentials(double time_since_last_mechanics,  mpi_Environment &world, mpi_Cartesian &cart_topo){

    // Boundary data must be exchanged before calling this routine.
    #pragma omp parallel for schedule(guided)
    for( int i=0; i < (*all_cells).size(); i++ )
    {
        Cell* pC = (*all_cells)[i];
        
        if( pC->functions.contact_function && pC->is_out_of_domain == false )
			{ evaluate_interactions( pC,pC->phenotype,time_since_last_mechanics ); } //No parallelism necessary yet


        if( pC->functions.custom_cell_rule && pC->is_out_of_domain == false )
			{ pC->functions.custom_cell_rule( pC,pC->phenotype,time_since_last_mechanics ); } //No parallelism necessary yet

        if( pC->is_movable && pC->is_out_of_domain == false )
        {
            if(pC->functions.update_velocity_parallel && pC->is_out_of_domain == false)
            {
                pC->functions.update_velocity_parallel( pC, pC->phenotype, time_since_last_mechanics, world, cart_topo);
            }
            if( PhysiCell_settings.disable_automated_spring_adhesions == false ) 
            {
                evaluate_cell_elastic_interactions( pC , time_since_last_mechanics, world, cart_topo);
            }
        }
    }	
}

void Cell_Container::pack(std::vector<Cell*> *all_cells, mpi_Environment &world, mpi_Cartesian &cart_topo)
{

	int len_snd_buf_left  	 = 0;
	int len_snd_buf_right 	 = 0;

	int 	temp_int;
	double 	temp_double;
	std::string temp_str;
	std::unordered_map<std::string, int> :: iterator it;
	Cell *pCell;


	int len_int  = 0;
	int len_double  = 0;
	int len_str = 0;
	int len_vector = 0;

	position_left 			= 0;					//Must be initialized to 0
	position_right 			= 0;					//Must be initialized to 0
	no_cells_cross_left  	= 0;					//Must be initialized to 0
	no_cells_cross_right	= 0;					//Must be initialized to 0
	no_of_cells_from_right 	= 0;
	no_of_cells_from_left 	= 0;

	snd_buf_left.resize(0); 					//When we enter function again, this is reset
	snd_buf_right.resize(0);					//Reset the others also
	rcv_buf_left.resize(0);
	rcv_buf_right.resize(0);
    snd_pos_left.resize(0);
    snd_pos_right.resize(0);
    rcv_pos_left.resize(0);
    rcv_pos_right.resize(0);


	/*--------------------------------------------------------------------------------*/
	/* First count the number of cells crossing to left subdomain and right subdomain */
	/*--------------------------------------------------------------------------------*/

	for(int i=0; i<(*all_cells).size();i++)
	{
		pCell = (*all_cells)[i];
		if(pCell->crossed_to_left_subdomain == true)
			no_cells_cross_left++;
		if(pCell->crossed_to_right_subdomain == true)
			no_cells_cross_right++;
	}

    snd_pos_left.resize(no_cells_cross_left);
    snd_pos_right.resize(no_cells_cross_right);
	/*-----------------------------------------------------------------------------------*/
	/* There is no need to pack crossed_to_left_subdomain and crossed_to_right_subdomain */
	/* as we send these two fields as the First Communication between MPI processes. 	 */
	/* Also try to send the length of the buffer to be allocated when sending.			 */
	/*-----------------------------------------------------------------------------------*/

  //std::cout<<"Total cells crossing to left in Rank "<<world.rank<<":"<<no_cells_cross_left<<std::endl;
	//std::cout<<"Total cells crossing to right in Rank "<<world.rank<<":"<<no_cells_cross_right<<std::endl;

	/* IMPORTANT: CANNOT USE #pragma omp for HERE AS ALL THREADS WILL WRITE TO THE SAME SHARED BUFFER */
	int left = 0;
	int right = 0;
	int right_prev, left_prev;
	left_prev = 0;
	right_prev = 0;
	for(int i=0; i<(*all_cells).size();i++)
	{
		pCell = (*all_cells)[i];

		if(pCell->crossed_to_left_subdomain == true)
		{
            snd_pos_left[left] = position_left;
			//pCell->pack(snd_buf_left, len_snd_buf_left, position_left);
			pCell->pack_light(snd_buf_left, len_snd_buf_left, position_left);
			left++;
		}
		else if(pCell->crossed_to_right_subdomain == true)
		{
            snd_pos_right[right] = position_right;
			//pCell->pack(snd_buf_right, len_snd_buf_right, position_right);
			pCell->pack_light(snd_buf_right, len_snd_buf_right, position_right);
			right++;
		}
	}
 }

void Cell_Container::unpack(mpi_Environment &world, mpi_Cartesian &cart_topo)
{

		/*---------------------------------------------------------------------------*/
		/* position_left/right are members of class Cell_Container and contain the 	 */
		/* size of the buffers to be sent to the left and right neighbour, respect.	 */
		/* Since their purpose gets over after packing and sending/receiving, they	 */
		/* can be used for storing the position we are at in the recv buffers				 */
		/* (See beginning of Cell_Container::pack() function for 10 data members		 */
		/* relevant to packing, sending/receiving, unpacking data)  								 */
		/*---------------------------------------------------------------------------*/

		position_left  = 0;
		position_right = 0;

		/* Temporary variables to help in unpacking */

		
		int size_left, size_right;
		int cell_ID = 0;
		double cell_position[3];

		int 	 				temp_int;
		double 				temp_double;
		std::string 	temp_str;
		char 					*temp_char_array;
		double 				*temp_double_array;
		int 					temp_key_value[2];

		std::unordered_map<std::string, int> :: iterator it;
		std::vector<double> temp_double_vector;

		int len_int 		 		= 0;
		int len_double 	 		= 0;
		int len_str 		 		= 0;
		int len_vector	 		= 0;
		int len_vector_nest = 0;
		
		int right = 0;
		int left = 0;
		auto start_time = std::chrono::high_resolution_clock::now();
		if(no_of_cells_from_right > 0)
		{
			size_right = rcv_buf_right.size();

			/*-------------------------------------------------------------------------------------------------*/
			/* IMPORTANT: need to start a for(int i=0; i<no_of_cells_from_right; i++) here BUT cannot put this */
			/* loop BEFORE completing the unpacking, as the next cell would get filled with wrong values 			 */
			/* so complete unpacking of 1 cell first THEN put the loop 																				 */
			/*-------------------------------------------------------------------------------------------------*/
                //#pragma omp parallel for
				for(int loop_ctr = 0; loop_ctr < no_of_cells_from_right; loop_ctr++)
				{
					Cell *pCell = create_cell(cell_ID);
                    
                    //int  pos = rcv_pos_right[loop_ctr]; 

					//pCell->unpack(rcv_buf_right, size_right, position_right);
					pCell->unpack_light(rcv_buf_right, size_right, position_right);
					++right;

					pCell->assign_position(pCell->position[0], pCell->position[1], pCell->position[2], world, cart_topo);
					/*
					std::cout << "\t\t\t\t[Rank " << world.rank << " ] Unpacked cell ID: " << pCell->ID << " at position: (" 
					          << pCell->position[0] << "," << pCell->position[1] << "," << pCell->position[2] << ")  velocity: (" 
					          << pCell->velocity[0] << "," << pCell->velocity[1] << "," << pCell->velocity[2] << ") from RIGHT" 
					          << std::endl;*/
		 		}	 
		}

		if(no_of_cells_from_left > 0)
		{
			 size_left = rcv_buf_left.size();

			/*-------------------------------------------------------------------------------------------------*/
			/* IMPORTANT: need to start a for(int i=0; i<no_of_cells_from_left; i++) here BUT cannot put this */
			/* loop BEFORE completing the unpacking, as the next cell would get filled with wrong values 			 */
			/* so complete unpacking of 1 cell first THEN put the loop 																				 */
			/*-------------------------------------------------------------------------------------------------*/
                //#pragma omp parallel for
				for(int loop_ctr = 0; loop_ctr < no_of_cells_from_left; loop_ctr++)
				{
					Cell *pCell = create_cell(cell_ID);

                    int  pos = rcv_pos_left[loop_ctr];

					//pCell->unpack(rcv_buf_left, size_left, position_left);
					pCell->unpack_light(rcv_buf_left, size_left, position_left);
					
					pCell->assign_position(pCell->position[0], pCell->position[1], pCell->position[2], world, cart_topo);
					++left;
					/*
					std::cout << "\t\t\t\t[Rank " << world.rank << " ] Unpacked cell ID: " << pCell->ID << " at position: (" 
					          << pCell->position[0] << "," << pCell->position[1] << "," << pCell->position[2] << ")  velocity: (" 
					          << pCell->velocity[0] << "," << pCell->velocity[1] << "," << pCell->velocity[2] << ") from LEFT" 
					          << std::endl;*/
				}
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
		/*
		if (no_of_cells_from_left + no_of_cells_from_right > 0) {
			std::cout << "\t [Rank " << world.rank << " ] Unpacking completed in avg time per cell " << duration.count()/(no_of_cells_from_left+no_of_cells_from_right) << " microseconds. " << std::endl;
			std::cout << "\t [Rank " << world.rank << " ] Unpacking completed total time " << duration.count() << " microseconds for cells: " << no_of_cells_from_left + no_of_cells_from_right << std::endl;
		}*/
}

void Cell_Container::unpack_parallel(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	int len_rcv_buf_left = rcv_buf_left.size();
	int len_rcv_buf_right = rcv_buf_right.size();

	std::vector<Cell*> new_cells_from_left(no_of_cells_from_left, nullptr);
	std::vector<Cell*> new_cells_from_right(no_of_cells_from_right, nullptr);
	auto start_time = std::chrono::high_resolution_clock::now();
	#pragma omp parallel for schedule(static)
	for( int i = 0; i < no_of_cells_from_left; i++ )
	{
		// Use the received byte offsets to unpack independent cells concurrently.
		Cell* pCell = new Cell(0);
		int pos = rcv_pos_left[i];
		//pCell->unpack(rcv_buf_left, len_rcv_buf_left, pos);
		pCell->unpack_light(rcv_buf_left, len_rcv_buf_left, pos);
		new_cells_from_left[i] = pCell;
	}

	#pragma omp parallel for schedule(static)
	for( int i = 0; i < no_of_cells_from_right; i++ )
	{
		Cell* pCell = new Cell(0);
		int pos = rcv_pos_right[i];
		//pCell->unpack(rcv_buf_right, len_rcv_buf_right, pos);
		pCell->unpack_light(rcv_buf_right, len_rcv_buf_right, pos);
		new_cells_from_right[i] = pCell;
	}
	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration_part = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
	/*
	if (no_of_cells_from_left + no_of_cells_from_right > 0) {
		std::cout << "\t [Rank " << world.rank << " ] Unpacking completed in avg time per cell " << duration_part.count()/(no_of_cells_from_left+no_of_cells_from_right) << " microseconds. " << std::endl;
		std::cout << "\t [Rank " << world.rank << " ] Unpacking completed total time " << duration_part.count() << " microseconds for cells: " << no_of_cells_from_left + no_of_cells_from_right << std::endl;
	}*/
	// Integration into global cell storage and mechanics voxels remains serialized
	// because all_cells and agent_grid are shared mutable containers.
	start_time = std::chrono::high_resolution_clock::now();
	for( int i = 0; i < no_of_cells_from_right; i++ )
	{
		Cell* pCell = new_cells_from_right[i];
		(*all_cells).push_back( pCell );
		pCell->index = (*all_cells).size() - 1;
		pCell->assign_position(pCell->position[0], pCell->position[1], pCell->position[2], world, cart_topo);
	}

	for( int i = 0; i < no_of_cells_from_left; i++ )
	{
		Cell* pCell = new_cells_from_left[i];
		(*all_cells).push_back( pCell );
		pCell->index = (*all_cells).size() - 1;
		pCell->assign_position(pCell->position[0], pCell->position[1], pCell->position[2], world, cart_topo);
	}
	end_time = std::chrono::high_resolution_clock::now();
	auto duration_integration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
	//if (no_of_cells_from_left + no_of_cells_from_right > 0) {
	//	std::cout << "\t [Rank " << world.rank << " ] Integration of unpacked cells completed in avg time per cell " << duration_integration.count()/(no_of_cells_from_left+no_of_cells_from_right) << " microseconds. " << std::endl;
	//	std::cout << "\t [Rank " << world.rank << " ] Integration of unpacked cells completed total time " << duration_integration.count() << " microseconds for cells: " << no_of_cells_from_left + no_of_cells_from_right << std::endl;
	//}
}

} // namespace PhysiCell
