#ifndef __log_paper_h__
#define __log_paper_h__

#include <string>

#include "../DistPhy/DistPhy.h"

namespace PhysiCell
{

// Collect one full-save interval across all ranks, append it on rank zero, and
// reset the local phase timers.
void log_paper( const std::string& csv_path, double simulation_time,
	DistPhy::mpi::mpi_Environment& world,
	DistPhy::mpi::mpi_Cartesian& cart_topo );

}

#endif
