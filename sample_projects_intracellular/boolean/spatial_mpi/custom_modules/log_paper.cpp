#include "log_paper.h"

#include "../core/PhysiCell.h"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <unistd.h>
#include <omp.h>

namespace
{

double process_memory_gb()
{
	std::ifstream statm( "/proc/self/statm" );
	long total_pages = 0;
	long resident_pages = 0;
	if( !( statm >> total_pages >> resident_pages ) )
	{ return 0.0; }

	const long page_size = sysconf( _SC_PAGESIZE );
	if( page_size <= 0 )
	{ return 0.0; }

	return static_cast<double>( resident_pages ) * static_cast<double>( page_size ) /
		( 1024.0 * 1024.0 * 1024.0 );
}

}

namespace PhysiCell
{

void log_paper( const std::string& csv_path, double simulation_time,
	DistPhy::mpi::mpi_Environment& world,
	DistPhy::mpi::mpi_Cartesian& cart_topo )
{
	using clock = std::chrono::steady_clock;
	static bool initialized = false;
	static std::string active_path;
	static std::ofstream csv;
	static clock::time_point previous_checkpoint;

	const clock::time_point now = clock::now();
	double local_elapsed_seconds = initialized
		? std::chrono::duration<double>( now - previous_checkpoint ).count() : 0.0;
	previous_checkpoint = now;

	double global_elapsed_seconds = 0.0;
	double global_diffusion_us = 0.0;
	double global_mechanics_us = 0.0;
	double global_phenotype_us = 0.0;
	double global_intracellular_us = 0.0;
	double local_memory_gb = process_memory_gb();
	double global_memory_gb = 0.0;
	unsigned long long local_cells = static_cast<unsigned long long>( all_cells->size() );
	unsigned long long global_cells = 0;

	MPI_Reduce( &local_elapsed_seconds, &global_elapsed_seconds, 1, MPI_DOUBLE,
		MPI_MAX, 0, cart_topo.mpi_cart_comm );
	MPI_Reduce( &time_diff, &global_diffusion_us, 1, MPI_DOUBLE,
		MPI_MAX, 0, cart_topo.mpi_cart_comm );
	MPI_Reduce( &time_mechs, &global_mechanics_us, 1, MPI_DOUBLE,
		MPI_MAX, 0, cart_topo.mpi_cart_comm );
	MPI_Reduce( &time_pheno, &global_phenotype_us, 1, MPI_DOUBLE,
		MPI_MAX, 0, cart_topo.mpi_cart_comm );
	MPI_Reduce( &time_intracell, &global_intracellular_us, 1, MPI_DOUBLE,
		MPI_MAX, 0, cart_topo.mpi_cart_comm );
	MPI_Reduce( &local_memory_gb, &global_memory_gb, 1, MPI_DOUBLE,
		MPI_SUM, 0, cart_topo.mpi_cart_comm );
	MPI_Reduce( &local_cells, &global_cells, 1, MPI_UNSIGNED_LONG_LONG,
		MPI_SUM, 0, cart_topo.mpi_cart_comm );

	int output_ok = 1;
	if( world.rank == 0 )
	{
		if( !initialized )
		{
			active_path = csv_path;
			csv.open( csv_path.c_str(), std::ios::out | std::ios::trunc );
			if( csv )
			{
				csv << "simulation_time,elapsed_wall_seconds,diffusion_seconds,mechanics_seconds,"
					"phenotype_seconds,intracellular_seconds,diffusion_percent,mechanics_percent,"
					"phenotype_percent,intracellular_percent,memory_gb,mpi_ranks,"
					"omp_threads_per_rank,total_cells\n";
			}
			else
			{ output_ok = 0; }
		}
		else if( csv_path != active_path )
		{ output_ok = 0; }

		if( output_ok )
		{
			const double diffusion_seconds = global_diffusion_us * 1e-6;
			const double mechanics_seconds = global_mechanics_us * 1e-6;
			const double phenotype_seconds = global_phenotype_us * 1e-6;
			const double intracellular_seconds = global_intracellular_us * 1e-6;
			const double timed_seconds = diffusion_seconds + mechanics_seconds +
				phenotype_seconds + intracellular_seconds;
			const double percent_scale = timed_seconds > 0.0 ? 100.0 / timed_seconds : 0.0;

			csv << std::setprecision( 17 )
				<< simulation_time << ',' << global_elapsed_seconds << ','
				<< diffusion_seconds << ',' << mechanics_seconds << ','
				<< phenotype_seconds << ',' << intracellular_seconds << ','
				<< diffusion_seconds * percent_scale << ','
				<< mechanics_seconds * percent_scale << ','
				<< phenotype_seconds * percent_scale << ','
				<< intracellular_seconds * percent_scale << ','
				<< global_memory_gb << ',' << world.size << ','
				<< omp_get_max_threads() << ',' << global_cells << '\n';
			csv.flush();
			if( !csv )
			{ output_ok = 0; }
		}
	}

	MPI_Bcast( &output_ok, 1, MPI_INT, 0, cart_topo.mpi_cart_comm );
	if( !output_ok )
	{ throw std::runtime_error( "Cannot write paper performance CSV: " + csv_path ); }

	initialized = true;
	time_diff = 0.0;
	time_mechs = 0.0;
	time_pheno = 0.0;
	time_intracell = 0.0;
	time_maboss_update = 0.0;
	count_maboss_updates = 0;
	count_mechanics_cells = 0;
	count_mechanics_neighbor_candidates = 0;
	count_mechanics_neighbor_interactions = 0;
#ifdef MECHS_TIME
	time_mechs_gradient_total = 0.0;
	time_mechs_halo_total = 0.0;
	time_mechs_potential_total = 0.0;
	time_mechs_interactions_total = 0.0;
	time_mechs_clear_dummy_total = 0.0;
	time_mechs_update_position_total = 0.0;
	time_mechs_pack_total = 0.0;
	time_mechs_transfer_total = 0.0;
	time_mechs_unpack_total = 0.0;
	time_mechs_voxels_update_total = 0.0;
	count_mechanics_timing_steps = 0;
#endif
}

}
