#include<mpi.h>
#include "DistPhy_Collective.h"

namespace DistPhy
{
	namespace mpi
	{	
		double distribute_global_sum(double x, mpi_Cartesian &cart_topo)
		{
			double global_sum;
			int count = 1;
			MPI_Allreduce(&x, &global_sum, count, MPI_DOUBLE, MPI_SUM, cart_topo.mpi_cart_comm);
			return (global_sum);			
		}
		
		int distribute_global_sum(int x, mpi_Cartesian &cart_topo)
		{
			int global_sum;
			int count = 1;
			MPI_Allreduce(&x, &global_sum, count, MPI_INT, MPI_SUM, cart_topo.mpi_cart_comm);
			return (global_sum);			
		}
		
		double distribute_global_max(double x, mpi_Cartesian &cart_topo)
		{
		  double global_max;
		  int count = 1;
		  MPI_Allreduce(&x, &global_max, count, MPI_DOUBLE, MPI_MAX, cart_topo.mpi_cart_comm); 
		  return (global_max);  
		}
		
		double distribute_global_min(double x, mpi_Cartesian &cart_topo)
		{
		  double global_min;
		  int count = 1;
		  MPI_Allreduce(&x, &global_min, count, MPI_DOUBLE, MPI_MIN, cart_topo.mpi_cart_comm); 
		  return (global_min);  
		}		
	}
}