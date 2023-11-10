/*
#############################################################################
# If you use BioFVM in your project, please cite BioFVM and the version     #
# number, such as below:                                                    #
#                                                                           #
# We solved the diffusion equations using BioFVM (Version 1.1.7) [1]        #
#                                                                           #
# [1] A. Ghaffarizadeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730 #
#                                                                           #
#############################################################################
#                                                                           #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)   #
#                                                                           #
# Copyright (c) 2015-2017, Paul Macklin and the BioFVM Project              #
# All rights reserved.                                                      #
#                                                                           #
# Redistribution and use in source and binary forms, with or without        #
# modification, are permitted provided that the following conditions are    #
# met:                                                                      #
#                                                                           #
# 1. Redistributions of source code must retain the above copyright notice, #
# this list of conditions and the following disclaimer.                     #
#                                                                           #
# 2. Redistributions in binary form must reproduce the above copyright      #
# notice, this list of conditions and the following disclaimer in the       #
# documentation and/or other materials provided with the distribution.      #
#                                                                           #
# 3. Neither the name of the copyright holder nor the names of its          #
# contributors may be used to endorse or promote products derived from this #
# software without specific prior written permission.                       #
#                                                                           #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED #
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           #
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER #
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  #
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       #
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              #
#                                                                           #
#############################################################################
*/

#include "BioFVM_solvers.h" 
#include "BioFVM_vector.h" 

#include <iostream>
#include <omp.h>

namespace BioFVM{

// do I even need this? 
void diffusion_decay_solver__constant_coefficients_explicit( Microenvironment& M, double dt )
{
	static bool precomputations_and_constants_done = false; 
	if( !precomputations_and_constants_done )
	{
		std::cout	<< std::endl << "Using solver: " << __FUNCTION__ << std::endl 
					<< "     (constant diffusion coefficient with explicit stepping, implicit decay) ... " << std::endl << std::endl;  

		if( M.mesh.uniform_mesh == true )
		{
			std::cout << "Uniform mesh detected! Consider switching to a more efficient method, such as " << std::endl  
			<< "     diffusion_decay_solver__constant_coefficients_explicit_uniform_mesh" << std::endl  
			<< std::endl; 
		}

		precomputations_and_constants_done = true; 
	}

	return; 
}

void diffusion_decay_solver__constant_coefficients_explicit_uniform_mesh( Microenvironment& M, double dt )
{
	static bool precomputations_and_constants_done = false; 
	if( !precomputations_and_constants_done )
	{
		std::cout	<< std::endl << "Using solver: " << __FUNCTION__ << std::endl 
					<< "     (constant diffusion coefficient with explicit stepping, implicit decay, uniform mesh) ... " << std::endl << std::endl;  

		if( M.mesh.regular_mesh == false )
		{ std::cout << "Error. This code is only supported for regular meshes." << std::endl; }

		precomputations_and_constants_done = true; 
	}

	return; 
}

/*
void diffusion_decay_solver__constant_coefficients_LOD_3D( Microenvironment& M, double dt )
{
	if( M.mesh.regular_mesh == false || M.mesh.Cartesian_mesh == false )
	{
		std::cout << "Error: This algorithm is written for regular Cartesian meshes. Try: other solvers!" << std::endl << std::endl; 
	return; 
	}

	// define constants and pre-computed quantities 
	
	if( !M.diffusion_solver_setup_done )
	{
		std::cout << std::endl << "Using method " << __FUNCTION__ << " (implicit 3-D LOD with Thomas Algorithm) ... " 
		<< std::endl << std::endl;  
		
		M.thomas_denomx.resize( M.mesh.x_coordinates.size() , M.zero );
		M.thomas_cx.resize( M.mesh.x_coordinates.size() , M.zero );

		M.thomas_denomy.resize( M.mesh.y_coordinates.size() , M.zero );
		M.thomas_cy.resize( M.mesh.y_coordinates.size() , M.zero );
		
		M.thomas_denomz.resize( M.mesh.z_coordinates.size() , M.zero );
		M.thomas_cz.resize( M.mesh.z_coordinates.size() , M.zero );

		M.thomas_i_jump = 1; 
		M.thomas_j_jump = M.mesh.x_coordinates.size(); 
		M.thomas_k_jump = M.thomas_j_jump * M.mesh.y_coordinates.size(); 

		M.thomas_constant1 =  M.diffusion_coefficients; // dt*D/dx^2 
		M.thomas_constant1a = M.zero; // -dt*D/dx^2; 
		M.thomas_constant2 =  M.decay_rates; // (1/3)* dt*lambda 
		M.thomas_constant3 = M.one; // 1 + 2*constant1 + constant2; 
		M.thomas_constant3a = M.one; // 1 + constant1 + constant2; 		
			
		M.thomas_constant1 *= dt; 
		M.thomas_constant1 /= M.mesh.dx; 
		M.thomas_constant1 /= M.mesh.dx; 

		M.thomas_constant1a = M.thomas_constant1; 
		M.thomas_constant1a *= -1.0; 

		M.thomas_constant2 *= dt; 
		M.thomas_constant2 /= 3.0; // for the LOD splitting of the source 

		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant2; 

		M.thomas_constant3a += M.thomas_constant1; 
		M.thomas_constant3a += M.thomas_constant2; 

		// Thomas solver coefficients 

		M.thomas_cx.assign( M.mesh.x_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomx.assign( M.mesh.x_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomx[0] = M.thomas_constant3a; 
		M.thomas_denomx[ M.mesh.x_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.x_coordinates.size() == 1 )
		{ M.thomas_denomx[0] = M.one; M.thomas_denomx[0] += M.thomas_constant2; } 

		M.thomas_cx[0] /= M.thomas_denomx[0]; 
		for( unsigned int i=1 ; i <= M.mesh.x_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomx[i] , M.thomas_constant1 , M.thomas_cx[i-1] ); 
			M.thomas_cx[i] /= M.thomas_denomx[i]; // the value at  size-1 is not actually used  
		}

		M.thomas_cy.assign( M.mesh.y_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomy.assign( M.mesh.y_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomy[0] = M.thomas_constant3a; 
		M.thomas_denomy[ M.mesh.y_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.y_coordinates.size() == 1 )
		{ M.thomas_denomy[0] = M.one; M.thomas_denomy[0] += M.thomas_constant2; } 

		M.thomas_cy[0] /= M.thomas_denomy[0]; 
		for( unsigned int i=1 ; i <= M.mesh.y_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomy[i] , M.thomas_constant1 , M.thomas_cy[i-1] ); 
			M.thomas_cy[i] /= M.thomas_denomy[i]; // the value at  size-1 is not actually used  
		}

		M.thomas_cz.assign( M.mesh.z_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomz.assign( M.mesh.z_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomz[0] = M.thomas_constant3a; 
		M.thomas_denomz[ M.mesh.z_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.z_coordinates.size() == 1 )
		{ M.thomas_denomz[0] = M.one; M.thomas_denomz[0] += M.thomas_constant2; } 

		M.thomas_cz[0] /= M.thomas_denomz[0]; 
		for( unsigned int i=1 ; i <= M.mesh.z_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomz[i] , M.thomas_constant1 , M.thomas_cz[i-1] ); 
			M.thomas_cz[i] /= M.thomas_denomz[i]; // the value at  size-1 is not actually used  
		}	

		M.diffusion_solver_setup_done = true; 
	}

	// x-diffusion 
	
	M.apply_dirichlet_conditions();
	#pragma omp parallel for 
	for( unsigned int k=0; k < M.mesh.z_coordinates.size() ; k++ )
	{
		for( unsigned int j=0; j < M.mesh.y_coordinates.size() ; j++ )
		{
			// Thomas solver, x-direction

			// remaining part of forward elimination, using pre-computed quantities 
			int n = M.voxel_index(0,j,k);
			(*M.p_density_vectors)[n] /= M.thomas_denomx[0]; 

			for( unsigned int i=1; i < M.mesh.x_coordinates.size() ; i++ )
			{
				n = M.voxel_index(i,j,k); 
				axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_i_jump] ); 
				(*M.p_density_vectors)[n] /= M.thomas_denomx[i]; 
			}

			for( int i = M.mesh.x_coordinates.size()-2 ; i >= 0 ; i-- )
			{
				n = M.voxel_index(i,j,k); 
				naxpy( &(*M.p_density_vectors)[n] , M.thomas_cx[i] , (*M.p_density_vectors)[n+M.thomas_i_jump] ); 
			}

		}
	}

	// y-diffusion 

	M.apply_dirichlet_conditions();
	#pragma omp parallel for 
	for( unsigned int k=0; k < M.mesh.z_coordinates.size() ; k++ )
	{
		for( unsigned int i=0; i < M.mesh.x_coordinates.size() ; i++ )
		{
   // Thomas solver, y-direction

	// remaining part of forward elimination, using pre-computed quantities 

	int n = M.voxel_index(i,0,k);
	(*M.p_density_vectors)[n] /= M.thomas_denomy[0]; 

	for( unsigned int j=1; j < M.mesh.y_coordinates.size() ; j++ )
	{
		n = M.voxel_index(i,j,k); 
		axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_j_jump] ); 
		(*M.p_density_vectors)[n] /= M.thomas_denomy[j]; 
	}

	// back substitution 
	// n = voxel_index( mesh.x_coordinates.size()-2 ,j,k); 

	for( int j = M.mesh.y_coordinates.size()-2 ; j >= 0 ; j-- )
	{
		n = M.voxel_index(i,j,k); 
		naxpy( &(*M.p_density_vectors)[n] , M.thomas_cy[j] , (*M.p_density_vectors)[n+M.thomas_j_jump] ); 
	}

  }
 }

 // z-diffusion 

	M.apply_dirichlet_conditions();
 #pragma omp parallel for 
 for( unsigned int j=0; j < M.mesh.y_coordinates.size() ; j++ )
 {
	 
  for( unsigned int i=0; i < M.mesh.x_coordinates.size() ; i++ )
  {
   // Thomas solver, y-direction

	// remaining part of forward elimination, using pre-computed quantities 

	int n = M.voxel_index(i,j,0);
	(*M.p_density_vectors)[n] /= M.thomas_denomz[0]; 

	// should be an empty loop if mesh.z_coordinates.size() < 2  
	for( unsigned int k=1; k < M.mesh.z_coordinates.size() ; k++ )
	{
		n = M.voxel_index(i,j,k); 
		axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_k_jump] ); 
		(*M.p_density_vectors)[n] /= M.thomas_denomz[k]; 
	}

	// back substitution 

	// should be an empty loop if mesh.z_coordinates.size() < 2 
	for( int k = M.mesh.z_coordinates.size()-2 ; k >= 0 ; k-- )
	{
		n = M.voxel_index(i,j,k); 
		naxpy( &(*M.p_density_vectors)[n] , M.thomas_cz[k] , (*M.p_density_vectors)[n+M.thomas_k_jump] ); 
		// n -= i_jump; 
	}
  }
 }
 
	M.apply_dirichlet_conditions();
	
	// reset gradient vectors 
//	M.reset_all_gradient_vectors(); 

	return; 
}
*/
/*--------------------------------------------------------------------------------------*/
/* Parallel version of the 3-D Thomas Solver, uses X-decomposition i.e. the X-direction */
/* of BioFVM going from left to right is divided into slices, each slice given to 1 MPI */
/* process. 
/*--------------------------------------------------------------------------------------*/
//Jose
void diffusion_decay_solver__constant_coefficients_LOD_3D( Microenvironment& M, double dt, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
        int granurality = 64;
        int size = world.Size();
        int rank = world.Rank();
	    MPI_Request send_req[granurality+1];
        MPI_Request recv_req[granurality+1];
        double t_strt_set, t_end_set;
        double t_strt_x, t_end_x;
        double t_strt_y, t_end_y;
        double t_strt_z, t_end_z;

        if (M.mesh.uniform_mesh == false || M.mesh.Cartesian_mesh == false)
        {
            std::cout << "Error: This algorithm is written for uniform Cartesian meshes. Try: other solvers!" << std::endl
                      << std::endl;
            return;
        }

        // define constants and pre-computed quantities

        if (!M.diffusion_solver_setup_done)
        {
            // t_strt_set = MPI_Wtime();
            // std::cout << std::endl << "Using method " << __FUNCTION__ << " (implicit 3-D LOD with Thomas Algorithm) ... "
            //<< std::endl << std::endl;

            /*-------------------------------------------------------------*/
            /* x_coordinates are of size local_x_nodes                     */
            /* (see function resize() of class Cartesian Mesh in           */
            /* BioFVM_parallel.cpp.                                        */
            /* Each line of Voxels going from left to right forms          */
            /* a tridiagonal system of Equations                           */
            /* Now these lines are going to split in the X decomposition   */
            /*-------------------------------------------------------------*/

            M.thomas_denomx.resize(M.mesh.x_coordinates.size(), M.zero); // sizeof(x_coordinates) = local_x_nodes, denomx is the main diagonal elements
            M.thomas_cx.resize(M.mesh.x_coordinates.size(), M.zero);     // Both b and c of tridiagonal matrix are equal, hence just one array needed

            /*-------------------------------------------------------------*/
            /* y_coordinates are of size of local_y_nodes.                 */
            /* Each line of Voxels going                                   */
            /* from bottom to top forms a tridiagonal system of Equations  */
            /*-------------------------------------------------------------*/

            M.thomas_denomy.resize(M.mesh.y_coordinates.size(), M.zero);
            M.thomas_cy.resize(M.mesh.y_coordinates.size(), M.zero);

            /*-------------------------------------------------------------*/
            /* z_coordinates are of size of local_z_nodes.                 */
            /* Each line of Voxels going                                   */
            /* from front to back forms a tridiagonal system of Equations  */
            /*-------------------------------------------------------------*/

            M.thomas_denomz.resize(M.mesh.z_coordinates.size(), M.zero);
            M.thomas_cz.resize(M.mesh.z_coordinates.size(), M.zero);

            /*-------------------------------------------------------------*/
            /* For X-decomposition thomas_i_jump - 1 can be in the previous*/
            /* process and thomas_i_jump+1 can be in the next processs     */
            /* hence we can use thomas_j_jump and thomas_k_jump safely     */
            /* but we CANNOT use thomas_i_jump safely                      */
            /*-------------------------------------------------------------*/

            M.thomas_i_jump = M.number_of_densities() * M.mesh.z_coordinates.size() * M.mesh.y_coordinates.size();
            M.thomas_j_jump = M.number_of_densities() * M.mesh.z_coordinates.size();
            M.thomas_k_jump = M.number_of_densities(); // M.thomas_j_jump * M.mesh.y_coordinates.size();

            /*-------------------------------------------------------------*/
            /* This part below of defining constants SHOULD typically      */
            /* not change during parallelization.                          */
            /*-------------------------------------------------------------*/

            M.thomas_constant1 = M.diffusion_coefficients; // dt*D/dx^2
            M.thomas_constant1a = M.zero;                  // -dt*D/dx^2;
            M.thomas_constant2 = M.decay_rates;            // (1/3)* dt*lambda
            M.thomas_constant3 = M.one;                    // 1 + 2*constant1 + constant2;
            M.thomas_constant3a = M.one;                   // 1 + constant1 + constant2;

            M.thomas_constant1 *= dt;
            M.thomas_constant1 /= M.mesh.dx;
            M.thomas_constant1 /= M.mesh.dx;

            M.thomas_constant1a = M.thomas_constant1;
            M.thomas_constant1a *= -1.0;

            M.thomas_constant2 *= dt;
            M.thomas_constant2 /= 3.0; // for the LOD splitting of the source, division by 3 is for 3-D

            M.thomas_constant3 += M.thomas_constant1;
            M.thomas_constant3 += M.thomas_constant1;
            M.thomas_constant3 += M.thomas_constant2;

            M.thomas_constant3a += M.thomas_constant1;
            M.thomas_constant3a += M.thomas_constant2;

            // Thomas solver coefficients

            /*--------------------------------------------------------------------*/
            /* In 1-D X decomposition, y and z-lines are contiguous and typically */
            /* the assignments below for y,z should not be changed                */
            /*--------------------------------------------------------------------*/

            M.thomas_cx.assign(M.mesh.x_coordinates.size(), M.thomas_constant1a);    // Fill b and c elements with -D * dt/dx^2
            M.thomas_denomx.assign(M.mesh.x_coordinates.size(), M.thomas_constant3); // Fill diagonal elements with (1 + 1/3 * lambda * dt + 2*D*dt/dx^2)

            if (rank == 0)
                M.thomas_denomx[0] = M.thomas_constant3a; // First diagonal element is   (1 + 1/3 * lambda * dt + 1*D*dt/dx^2)

            if (rank == (size - 1))
                M.thomas_denomx[M.mesh.x_coordinates.size() - 1] = M.thomas_constant3a; // Last diagonal element  is   (1 + 1/3 * lambda * dt + 1*D*dt/dx^2)

            if (rank == 0)
                if (M.mesh.x_coordinates.size() == 1) // This is an extreme case, won't exist, still if it does
                {                                     // then this must be at rank 0
                    M.thomas_denomx[0] = M.one;
                    M.thomas_denomx[0] += M.thomas_constant2;
                }
            if (rank == 0)
                M.thomas_cx[0] /= M.thomas_denomx[0]; // The first c element of tridiagonal matrix is div by first diagonal el.

            // axpy(1st, 2nd, 3rd) => 1st = 1st + 2nd * 3rd
            // the value at  size-1 is not actually used
            // Since value of size-1 is not used, it means it is the value after the last Diagonal element
            // cout << "Rank " << rank << endl;
            for (int ser_ctr = 0; ser_ctr <= size - 1; ser_ctr++)
            {
                if (rank == ser_ctr)
                {
                    if (rank == 0 && rank <= size - 1) // If size=1, then this process does not send data
                    {

                        for (int i = 1; i <= M.mesh.x_coordinates.size() - 1; i++)
                        {
                            axpy(&M.thomas_denomx[i], M.thomas_constant1, M.thomas_cx[i - 1]);
                            M.thomas_cx[i] /= M.thomas_denomx[i]; // the value at  size-1 is not actually used
                        }
                    }
                    else
                    {
                        for (int i = 1; i <= M.mesh.x_coordinates.size() - 1; i++)
                        {
                            axpy(&M.thomas_denomx[i], M.thomas_constant1, M.thomas_cx[i - 1]);
                            M.thomas_cx[i] /= M.thomas_denomx[i]; // the value at  size-1 is not actually used
                        }
                    }

                    if (rank < (size - 1))
                    {
                        MPI_Isend(&(M.thomas_cx[M.mesh.x_coordinates.size() - 1][0]), M.thomas_cx[M.mesh.x_coordinates.size() - 1].size(), MPI_DOUBLE, ser_ctr + 1, 1111, mpi_Cart_comm, &send_req[0]);
                    }
                }

                if (rank == (ser_ctr + 1) && (ser_ctr + 1) <= (size - 1))
                {

                    std::vector<double> temp_cx(M.thomas_cx[0].size());

                    MPI_Irecv(&temp_cx[0], temp_cx.size(), MPI_DOUBLE, ser_ctr, 1111, mpi_Cart_comm, &recv_req[0]);
                    MPI_Wait(&recv_req[0], MPI_STATUS_IGNORE);

                    axpy(&M.thomas_denomx[0], M.thomas_constant1, temp_cx); // CHECK IF &temp_cz[0] is OK, axpy() in BioFVM_vector.cpp
                    M.thomas_cx[0] /= M.thomas_denomx[0];                   // the value at  size-1 is not actually used
                }

                MPI_Barrier(mpi_Cart_comm);
            }
            cout << "Diffusion set up is done" << endl;

            /*--------------------------------------------------------------------*/
            /* In 1-D X decomposition, z and y-lines are contiguous adn typically */
            /* the assignments below for z,y should not be changed                */
            /* Both the first voxel i.e. index 0 and last voxel i.e. index=       */
            /* y_coordinates.size()-1 are on the same process                     */
            /*--------------------------------------------------------------------*/

            M.thomas_cy.assign(M.mesh.y_coordinates.size(), M.thomas_constant1a);
            M.thomas_denomy.assign(M.mesh.y_coordinates.size(), M.thomas_constant3);
            M.thomas_denomy[0] = M.thomas_constant3a;
            M.thomas_denomy[M.mesh.y_coordinates.size() - 1] = M.thomas_constant3a;
            if (M.mesh.y_coordinates.size() == 1)
            {
                M.thomas_denomy[0] = M.one;
                M.thomas_denomy[0] += M.thomas_constant2;
            }
            M.thomas_cy[0] /= M.thomas_denomy[0];
            for (int i = 1; i <= M.mesh.y_coordinates.size() - 1; i++)
            {
                axpy(&M.thomas_denomy[i], M.thomas_constant1, M.thomas_cy[i - 1]);
                M.thomas_cy[i] /= M.thomas_denomy[i]; // the value at  size-1 is not actually used
            }

            M.thomas_cz.assign(M.mesh.z_coordinates.size(), M.thomas_constant1a);
            M.thomas_denomz.assign(M.mesh.z_coordinates.size(), M.thomas_constant3);
            M.thomas_denomz[0] = M.thomas_constant3a;
            M.thomas_denomz[M.mesh.z_coordinates.size() - 1] = M.thomas_constant3a;
            if (M.mesh.z_coordinates.size() == 1)
            {
                M.thomas_denomz[0] = M.one;
                M.thomas_denomz[0] += M.thomas_constant2;
            }
            M.thomas_cz[0] /= M.thomas_denomz[0];
            for (int i = 1; i <= M.mesh.z_coordinates.size() - 1; i++)
            {
                axpy(&M.thomas_denomz[i], M.thomas_constant1, M.thomas_cz[i - 1]);
                M.thomas_cz[i] /= M.thomas_denomz[i]; // the value at  size-1 is not actually used
            }

            M.diffusion_solver_setup_done = true;
            // t_end_set = MPI_Wtime();
            // std::cout<<"Set-up time = "<<(t_end_set-t_strt_set)<<std::endl;

            
            //if (rank == 0) file << "X-diffusion,Y-diffusion,Z-diffusion,Apply Dirichlet" << std::endl;
        }


        // x-diffusion
        // cout << "Rank " << rank << " starting x-diffusion" << endl;
        //auto start_time = std::chrono::high_resolution_clock::now();
        M.apply_dirichlet_conditions(rank, size);
        //auto end_time = std::chrono::high_resolution_clock::now();
        //auto apply_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();


        // cout << "Rank " << rank << " apply dirichlet condition done" << endl;
        /*-----------------------------------------------------------------------------------*/
        /*                        FORWARD ELIMINATION - x DIRECTION/DECOMPOSITION            */
        /*-----------------------------------------------------------------------------------*/

        /* For data packing...                                                                                 */
        /* My direction of traversing is go up up up i.e. y direction points then go in i.e. Z-direction       */
        /* Remember to visualize 3-D as 2-D plates kept after one another. Hence Z-direction data is farther   */
        /* apart than X/Y direction                                                                            */

        int y_size = M.mesh.y_coordinates.size();
        int z_size = M.mesh.z_coordinates.size();
        int p_size = M.number_of_densities(); 

        int step_size = (z_size * y_size) / granurality;

        int snd_data_size = step_size * p_size; // Number of data elements to be sent
        int rcv_data_size = step_size * p_size; // All p_density_vectors elements have same size, use anyone

        int snd_data_size_last = ((z_size * y_size) % granurality) * p_size; // Number of data elements to be sent
        int rcv_data_size_last = ((z_size * y_size) % granurality) * p_size;
        bool last_iteration = ((z_size * y_size) % granurality) > 0;
        // cout << "Rank " << rank << " snd_data_size: " << snd_data_size << " rcv_data_size: " << rcv_data_size << endl;

        /* So row is along Z axis, column of each row is along Y-axis and each element has p_density_vector*/

        std::vector<double> block3d(z_size * y_size * p_size);

        // t_strt_x = MPI_Wtime();
        // cout << "Rank " << rank << " starting forward" << endl;
        //start_time = std::chrono::high_resolution_clock::now();
        
        if (rank == 0)
        {
            // cout << "Rank " << rank << "is computing" << endl;
            for (int step = 0; step < granurality; ++step)
            {
                int initial_index = step * snd_data_size;
                #pragma omp parallel for
                for (int index = initial_index; index < initial_index + snd_data_size; index += p_size)
                {
                    int index_dec = index; 
                    for (int d = 0; d < M.thomas_denomx[0].size(); d++)
                    {
                        M.p_density_vectors[index + d] /= M.thomas_denomx[0][d];
                    }

                    for (int i = 1; i < M.mesh.x_coordinates.size(); i++)
                    {
                        
                        int index_inc = index_dec + M.thomas_i_jump;
                        // axpy(&(*M.p_density_vectors)[n], M.thomas_constant1, (*M.p_density_vectors)[n - M.thomas_i_jump]);
                        for (int d = 0; d < M.thomas_k_jump; d++)
                        {
                            M.p_density_vectors[index_inc + d] += M.thomas_constant1[d] * M.p_density_vectors[index_dec + d];
                        }

                        //(*M.p_density_vectors)[n] /= M.thomas_denomx[i];
                        for (int d = 0; d < M.thomas_denomx[i].size(); d++)
                        {
                            M.p_density_vectors[index_inc + d] /= M.thomas_denomx[i][d];
                        }
                        index_dec = index_inc;
                    }
                }

                if (size > 1) {
                    int x_end = M.mesh.x_coordinates.size() - 1;
                    int offset = step * snd_data_size;
                    MPI_Status status;
                    MPI_Isend(&(M.p_density_vectors[x_end * M.thomas_i_jump + offset]), snd_data_size, MPI_DOUBLE, rank + 1, step, mpi_Cart_comm, &send_req[step]);
                }
            }
            //Last iteration
            if (last_iteration) {
                int initial_index = granurality * snd_data_size;
                #pragma omp parallel for
                for (int index = initial_index; index < initial_index + snd_data_size_last; index += p_size)
                {
                    int index_dec = index; 
                    for (int d = 0; d < M.thomas_denomx[0].size(); d++)
                    {
                        M.p_density_vectors[index + d] /= M.thomas_denomx[0][d];
                    }

                    for (int i = 1; i < M.mesh.x_coordinates.size(); i++)
                    {
                        int index_inc = index_dec + M.thomas_i_jump;
                        // axpy(&(*M.p_density_vectors)[n], M.thomas_constant1, (*M.p_density_vectors)[n - M.thomas_i_jump]);
                        for (int d = 0; d < M.thomas_k_jump; d++)
                        {
                            M.p_density_vectors[index_inc + d] += M.thomas_constant1[d] * M.p_density_vectors[index_dec + d];
                        }

                        //(*M.p_density_vectors)[n] /= M.thomas_denomx[i];
                        for (int d = 0; d < M.thomas_denomx[i].size(); d++)
                        {
                            M.p_density_vectors[index_inc + d] /= M.thomas_denomx[i][d];
                        }
                        index_dec = index_inc;
                    }
                }

                if (size > 1) {
                    int x_end = M.mesh.x_coordinates.size() - 1;
                    int offset = granurality * snd_data_size;
                    MPI_Status status;
                    MPI_Isend(&(M.p_density_vectors[x_end * M.thomas_i_jump + offset]), snd_data_size_last, MPI_DOUBLE, rank + 1, granurality, mpi_Cart_comm, &send_req[granurality]);
                    
                }
            }
        }
        else
        {
            if (rank >= 1 && rank <= (size - 1))
            {
                for (int step = 0; step < granurality; ++step)
                {
                    int initial_index = step * snd_data_size;
                    MPI_Irecv(&(block3d[initial_index]), rcv_data_size, MPI_DOUBLE, rank-1, step, mpi_Cart_comm, &(recv_req[step]));
                }
                if (last_iteration)
                    MPI_Irecv(&(block3d[granurality*snd_data_size]), rcv_data_size_last, MPI_DOUBLE, rank-1, granurality, mpi_Cart_comm, &(recv_req[granurality]));
                for (int step = 0; step < granurality; ++step)
                {
                    int initial_index = step * snd_data_size;
                    MPI_Wait(&recv_req[step], MPI_STATUS_IGNORE);
                    #pragma omp parallel for
                    for (int index = initial_index; index < initial_index + snd_data_size; index += p_size)
                    {
                        // axpy(&(*M.p_density_vectors)[n], M.thomas_constant1, block3d[k][j]);
                        int index_dec = index;
                        for (int d = 0; d < M.thomas_k_jump; d++)
                        {
                            M.p_density_vectors[index + d] += M.thomas_constant1[d] * block3d[index + d];
                        }

                        //(*M.p_density_vectors)[n] /= M.thomas_denomx[0];
                        for (int d = 0; d < M.thomas_denomx[0].size(); d++)
                        {
                            M.p_density_vectors[index + d] /= M.thomas_denomx[0][d];
                        }

                        for (int i = 1; i < M.mesh.x_coordinates.size(); i++)
                        {

                            int index_inc = index_dec + M.thomas_i_jump;
                            // axpy(&(*M.p_density_vectors)[n], M.thomas_constant1, (*M.p_density_vectors)[n - M.thomas_i_jump]);
                            for (int d = 0; d < M.thomas_k_jump; d++)
                            {
                                M.p_density_vectors[index_inc + d] += M.thomas_constant1[d] * M.p_density_vectors[index_dec + d];
                            }
                            //(*M.p_density_vectors)[n] /= M.thomas_denomx[i];
                            for (int d = 0; d < M.thomas_denomx[i].size(); d++)
                            {
                                M.p_density_vectors[index_inc + d] /= M.thomas_denomx[i][d];
                            }

                            index_dec = index_inc;
                        }
                        
                    }
                    if (rank < (size - 1))
                    {
                        int x_end = M.mesh.x_coordinates.size() - 1;
                        MPI_Isend(&(M.p_density_vectors[x_end * M.thomas_i_jump + initial_index]), snd_data_size, MPI_DOUBLE, rank + 1, step, mpi_Cart_comm, &send_req[step]);
                    }
                }
                if (last_iteration)
                {
                    int initial_index = granurality * snd_data_size;
                    MPI_Wait(&recv_req[granurality], MPI_STATUS_IGNORE); //Need to change
                    #pragma omp parallel for
                    for (int index = initial_index; index < initial_index + snd_data_size_last; index += p_size)
                    {
                        // axpy(&(*M.p_density_vectors)[n], M.thomas_constant1, block3d[k][j]);
                        int index_dec = index;
                        for (int d = 0; d < M.thomas_k_jump; d++)
                        {
                            M.p_density_vectors[index + d] += M.thomas_constant1[d] * block3d[index + d];
                        }

                        //(*M.p_density_vectors)[n] /= M.thomas_denomx[0];
                        for (int d = 0; d < M.thomas_denomx[0].size(); d++)
                        {
                            M.p_density_vectors[index + d] /= M.thomas_denomx[0][d];
                        }

                        for (int i = 1; i < M.mesh.x_coordinates.size(); i++)
                        {

                            int index_inc = index_dec + M.thomas_i_jump;
                            // axpy(&(*M.p_density_vectors)[n], M.thomas_constant1, (*M.p_density_vectors)[n - M.thomas_i_jump]);
                            for (int d = 0; d < M.thomas_k_jump; d++)
                            {
                                M.p_density_vectors[index_inc + d] += M.thomas_constant1[d] * M.p_density_vectors[index_dec + d];
                            }
                            //(*M.p_density_vectors)[n] /= M.thomas_denomx[i];
                            for (int d = 0; d < M.thomas_denomx[i].size(); d++)
                            {
                                M.p_density_vectors[index_inc + d] /= M.thomas_denomx[i][d];
                            }

                            index_dec = index_inc;
                        }
                        
                    }
                    // End of computation region
                    if (rank < (size - 1))
                    {
                        int x_end = M.mesh.x_coordinates.size() - 1;
                        MPI_Request aux;
                        MPI_Isend(&(M.p_density_vectors[x_end * M.thomas_i_jump + initial_index]), snd_data_size_last, MPI_DOUBLE, rank + 1, granurality, mpi_Cart_comm, &send_req[granurality]);
                      
                    }
                }
            }
        }
        
        /*-----------------------------------------------------------------------------------*/
        /*                         CODE FOR BACK SUBSITUTION                                 */
        /*-----------------------------------------------------------------------------------*/
        
        //cout << "Rank " << rank << " starting backward substitution" << endl;

        if (rank == (size - 1))
        {
            for (int step = 0; step < granurality; ++step)
            {
                int initial_index = ((M.mesh.x_coordinates.size() - 1)*M.thomas_i_jump) + (step * snd_data_size);
                #pragma omp parallel for

                for (int index = initial_index; index < initial_index + snd_data_size; index += p_size)
                {
                    int index_aux = index;
                    //int index = j * M.thomas_j_jump + k * M.thomas_k_jump + (M.mesh.x_coordinates.size() - 1) * M.thomas_i_jump;
                    for (int i = M.mesh.x_coordinates.size() - 2; i >= 0; i--)
                    {

                        int index_dec = index_aux - M.thomas_i_jump;
                        // naxpy(&(*M.p_density_vectors)[n], M.thomas_cx[i], (*M.p_density_vectors)[n + M.thomas_i_jump]);
                        for (int d = 0; d < M.thomas_k_jump; d++)
                        {
                            M.p_density_vectors[index_dec + d] -= M.thomas_cx[i][d] * M.p_density_vectors[index_aux + d];
                        }
                        index_aux = index_dec;
                    }
                }
                if (size > 1) {
                    MPI_Request aux;
                    MPI_Isend(&(M.p_density_vectors[step * snd_data_size]), snd_data_size, MPI_DOUBLE, rank - 1, step, mpi_Cart_comm, &send_req[step]);
                }
            }

            //Last iteration
            if (last_iteration) {
                int initial_index = ((M.mesh.x_coordinates.size() - 1)*M.thomas_i_jump) + (granurality * snd_data_size);
                #pragma omp parallel for
                for (int index = initial_index; index < initial_index + snd_data_size_last; index += p_size)
                {
                    int index_aux = index;
                    for (int i = M.mesh.x_coordinates.size() - 2; i >= 0; i--)
                    {

                        int index_dec = index_aux - M.thomas_i_jump;
                        // naxpy(&(*M.p_density_vectors)[n], M.thomas_cx[i], (*M.p_density_vectors)[n + M.thomas_i_jump]);
                        for (int d = 0; d < M.thomas_k_jump; d++)
                        {
                            M.p_density_vectors[index_dec + d] -= M.thomas_cx[i][d] * M.p_density_vectors[index_aux + d];
                        }
                        index_aux = index_dec;
                    }
                }
                if (size > 1) {
                    MPI_Request aux;
                    MPI_Isend(&(M.p_density_vectors[granurality * snd_data_size]), snd_data_size_last, MPI_DOUBLE, rank - 1, granurality, mpi_Cart_comm, &send_req[granurality]);
                    //cout << "Rank " << rank << " has send" << endl;
                }
            
            }
        }
        else
        {
            MPI_Request status[granurality]; 
            MPI_Request status_last;
            for (int step = 0; step < granurality; ++step)
            {
                MPI_Irecv(&(block3d[step*snd_data_size]), rcv_data_size, MPI_DOUBLE, rank+1, step, mpi_Cart_comm, &recv_req[step]);
            }
            if (last_iteration)
                MPI_Irecv(&(block3d[granurality*snd_data_size]), rcv_data_size_last, MPI_DOUBLE, rank+1, granurality, mpi_Cart_comm, &recv_req[granurality]);

            
            for (int step = 0; step < granurality; ++step)
            {
                int initial_index = ((M.mesh.x_coordinates.size() - 1)*M.thomas_i_jump) + (step * snd_data_size);
                int index_3d_initial = (step * snd_data_size);
                MPI_Wait(&recv_req[step], MPI_STATUS_IGNORE);
                #pragma omp parallel for
                for (int offset = 0; offset < snd_data_size; offset += p_size)
                {
                    int index_aux = initial_index + offset;
                    int index_3d = index_3d_initial + offset;
                    for (int d = 0; d < M.thomas_k_jump; d++)
                    {
                        M.p_density_vectors[index_aux + d] -= M.thomas_cx[M.mesh.x_coordinates.size() - 1][d] * block3d[index_3d + d];
                    }

                    for (int i = M.mesh.x_coordinates.size() - 2; i >= 0; i--)
                    {

                        int index_dec = index_aux - M.thomas_i_jump;
                        // naxpy(&(*M.p_density_vectors)[n], M.thomas_cx[i], (*M.p_density_vectors)[n + M.thomas_i_jump]);
                        for (int d = 0; d < M.thomas_k_jump; d++)
                        {
                            M.p_density_vectors[index_dec + d] -= M.thomas_cx[i][d] * M.p_density_vectors[index_aux + d];
                        }
                        index_aux = index_dec;
                        
                    }
                }
                if (rank > 0)
                {
                    MPI_Request aux;
                    MPI_Isend(&(M.p_density_vectors[step * snd_data_size]), snd_data_size, MPI_DOUBLE, rank - 1, step, mpi_Cart_comm, &send_req[step]);
                    // cout << "Rank " << rank << " has send" << endl;
                }
            }
            if (last_iteration)
            {
                int initial_index = ((M.mesh.x_coordinates.size() - 1)*M.thomas_i_jump) + (granurality * snd_data_size);
                int index_3d_initial = (granurality * snd_data_size);
                MPI_Wait(&recv_req[granurality], MPI_STATUS_IGNORE);
                #pragma omp parallel for
                for (int offset = 0; offset < snd_data_size_last; offset += p_size)
                {
                    int index_aux = initial_index + offset;
                    //int index = j * M.thomas_j_jump + k * M.thomas_k_jump + (M.mesh.x_coordinates.size() - 1) * M.thomas_i_jump;
                    int index_3d = index_3d_initial + offset;
                    // naxpy(&(*M.p_density_vectors)[n], M.thomas_cx[M.mesh.x_coordinates.size() - 1], block3d[k][j]);
                    for (int d = 0; d < M.thomas_k_jump; d++)
                    {
                        M.p_density_vectors[index_aux + d] -= M.thomas_cx[M.mesh.x_coordinates.size() - 1][d] * block3d[index_3d + d];
                    }

                    for (int i = M.mesh.x_coordinates.size() - 2; i >= 0; i--)
                    {

                        int index_dec = index_aux - M.thomas_i_jump;
                        // naxpy(&(*M.p_density_vectors)[n], M.thomas_cx[i], (*M.p_density_vectors)[n + M.thomas_i_jump]);
                        for (int d = 0; d < M.thomas_k_jump; d++)
                        {
                            M.p_density_vectors[index_dec + d] -= M.thomas_cx[i][d] * M.p_density_vectors[index_aux + d];
                        }
                        index_aux = index_dec;
                        
                    }
                }
                if (rank > 0)
                {
                    MPI_Request aux;
                    MPI_Isend(&(M.p_density_vectors[granurality * snd_data_size]), snd_data_size_last, MPI_DOUBLE, rank - 1, granurality, mpi_Cart_comm, &send_req[granurality]);
                }
            }
        }
        MPI_Barrier(mpi_Cart_comm);

        //end_time = std::chrono::high_resolution_clock::now();
        //auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
        //if (rank == 0)
            //std::cout << "   X diffusion: " << duration_ms << "ms" << std::endl;
         //   file << duration_us << ",";

       
        //start_time = std::chrono::high_resolution_clock::now();
        M.apply_dirichlet_conditions(rank, size);
        //end_time = std::chrono::high_resolution_clock::now();
        //apply_us += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
        //if (rank == 0)
            //std::cout << "   Apply dirichlet: " << duration_ms << "ms" << std::endl;
            //std::ofstream << duration_us << ",";

        // cout << "Rank " << rank << " Y diffusion" << endl;
        // t_strt_y = MPI_Wtime();
        //start_time = std::chrono::high_resolution_clock::now();

#pragma omp parallel for collapse(2)
        for (int k = 0; k < M.mesh.z_coordinates.size(); k++)
        {
            for (int i = 0; i < M.mesh.x_coordinates.size(); i++)
            {

                int index = i * M.thomas_i_jump + k * M.thomas_k_jump;
                //(*M.p_density_vectors)[n] /= M.thomas_denomy[0];
                for (int d = 0; d < M.thomas_denomy[0].size(); d++)
                {
                    M.p_density_vectors[index + d] /= M.thomas_denomy[0][d];
                }

                for (int j = 1; j < M.mesh.y_coordinates.size(); j++)
                {

                    int index_inc = index + M.thomas_j_jump;
                    // axpy(&(*M.p_density_vectors)[n], M.thomas_constant1, (*M.p_density_vectors)[n - M.thomas_j_jump]);
                    for (int d = 0; d < M.thomas_k_jump; d++)
                    {
                        M.p_density_vectors[index_inc + d] += M.thomas_constant1[d] * M.p_density_vectors[index + d];
                    }
                    //(*M.p_density_vectors)[n] /= M.thomas_denomy[j];
                    for (int d = 0; d < M.thomas_denomy[j].size(); d++)
                    {
                        M.p_density_vectors[index_inc + d] /= M.thomas_denomy[j][d];
                    }
                    index = index_inc;
                }

                // back substitution

                index = i * M.thomas_i_jump + k * M.thomas_k_jump + (M.thomas_j_jump * (M.mesh.y_coordinates.size() - 1));
                for (int j = M.mesh.y_coordinates.size() - 2; j >= 0; j--)
                {

                    int index_dec = index - M.thomas_j_jump;
                    // naxpy(&(*M.p_density_vectors)[n], M.thomas_cy[j], (*M.p_density_vectors)[n + M.thomas_j_jump]);
                    for (int d = 0; d < M.thomas_k_jump; d++)
                    {
                        M.p_density_vectors[index_dec + d] -= M.thomas_cy[j][d] * M.p_density_vectors[index + d];
                    }
                    index = index_dec;
                }
            }
        }
        //end_time = std::chrono::high_resolution_clock::now();
        //duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
        //if (rank == 0)
            //std::cout << "   Y diffussion: " << duration_ms << "ms" << std::endl;
        //    file << duration_us << ",";
        // t_end_y = MPI_Wtime();
        // std::cout<<"Y solve time = "<<(t_end_y-t_strt_y)<<std::endl;


        // cout << "Rank " << rank << " apply dirichlet" << endl;
         



        
        //start_time = std::chrono::high_resolution_clock::now();
        M.apply_dirichlet_conditions(rank, size);
        //end_time = std::chrono::high_resolution_clock::now();
        //apply_us += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
     
        
        //start_time = std::chrono::high_resolution_clock::now();
#pragma omp parallel for collapse(2)
        for (int j = 0; j < M.mesh.y_coordinates.size(); j++)
        {

            for (int i = 0; i < M.mesh.x_coordinates.size(); i++)
            {

                int index = i * M.thomas_i_jump + j * M.thomas_j_jump;
                //(*M.p_density_vectors)[n] /= M.thomas_denomz[0];
                for (int d = 0; d < M.thomas_denomz[0].size(); d++)
                {
                    M.p_density_vectors[index + d] /= M.thomas_denomz[0][d];
                }

                // should be an empty loop if mesh.z_coordinates.size() < 2
                for (int k = 1; k < M.mesh.z_coordinates.size(); k++)
                {

                    int index_inc = index + M.thomas_k_jump;
                    // axpy(&(*M.p_density_vectors)[n], M.thomas_constant1, (*M.p_density_vectors)[n - M.thomas_k_jump]);
                    for (int d = 0; d < M.thomas_k_jump; d++)
                    {
                        M.p_density_vectors[index_inc + d] += M.thomas_constant1[d] * M.p_density_vectors[index + d];
                    }
                    //(*M.p_density_vectors)[n] /= M.thomas_denomz[k];
                    for (int d = 0; d < M.thomas_denomz[k].size(); d++)
                    {
                        M.p_density_vectors[index_inc + d] /= M.thomas_denomz[k][d];
                    }

                    index = index_inc;
                }

                index = i * M.thomas_i_jump + j * M.thomas_j_jump + (M.thomas_k_jump * (M.mesh.z_coordinates.size() - 1));
                for (int k = M.mesh.z_coordinates.size() - 2; k >= 0; k--)
                {

                    int index_dec = index - M.thomas_k_jump;
                    // naxpy(&(*M.p_density_vectors)[n], M.thomas_cz[k], (*M.p_density_vectors)[n + M.thomas_k_jump]);
                    for (int d = 0; d < M.thomas_k_jump; d++)
                    {
                        M.p_density_vectors[index_dec + d] -= M.thomas_cz[k][d] * M.p_density_vectors[index + d];
                    }
                    index = index_dec;
                }
            }
        }
		/*
        end_time = std::chrono::high_resolution_clock::now();
        duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
        if (rank == 0)
            file << duration_us << ",";
		*/
    
        //start_time = std::chrono::high_resolution_clock::now();
        M.apply_dirichlet_conditions(rank, size);
        /*end_time = std::chrono::high_resolution_clock::now();
        apply_us += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
        if (rank == 0)
            file << apply_us/4.0 << std::endl;
		*/
        
        return;
}

//Jose
void diffusion_decay_solver__constant_coefficients_LOD_2D( Microenvironment& M, double dt )
{
	if( M.mesh.regular_mesh == false )
	{
		std::cout << "Error: This algorithm is written for regular Cartesian meshes. Try: something else." << std::endl << std::endl; 
		return; 
	}
	std::cout << "Diffusion decay solver constant coefficients LOD 2D is desactivated" << std::endl;
	// constants for the linear solver (Thomas algorithm) 
	/*
	if( !M.diffusion_solver_setup_done )
	{
		std::cout << std::endl << "Using method " << __FUNCTION__ << " (2D LOD with Thomas Algorithm) ... " << std::endl << std::endl;  
		
		M.thomas_denomx.resize( M.mesh.x_coordinates.size() , M.zero );
		M.thomas_cx.resize( M.mesh.x_coordinates.size() , M.zero );

		M.thomas_denomy.resize( M.mesh.y_coordinates.size() , M.zero );
		M.thomas_cy.resize( M.mesh.y_coordinates.size() , M.zero );
		
		// define constants and pre-computed quantities 

		M.thomas_i_jump = M.number_of_densities(); 
		M.thomas_j_jump = M.mesh.x_coordinates.size(); 

		M.thomas_constant1 =  M.diffusion_coefficients; //   dt*D/dx^2 
		M.thomas_constant1a = M.zero; // -dt*D/dx^2; 
		M.thomas_constant2 =  M.decay_rates; // (1/2)*dt*lambda 
		M.thomas_constant3 = M.one; // 1 + 2*constant1 + constant2; 
		M.thomas_constant3a = M.one; // 1 + constant1 + constant2; 
		
		M.thomas_constant1 *= dt; 
		M.thomas_constant1 /= M.mesh.dx; 
		M.thomas_constant1 /= M.mesh.dx; 

		M.thomas_constant1a = M.thomas_constant1; 
		M.thomas_constant1a *= -1.0; 

		M.thomas_constant2 *= dt; 
		M.thomas_constant2 *= 0.5; // for splitting via LOD

		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant2; 

		M.thomas_constant3a += M.thomas_constant1; 
		M.thomas_constant3a += M.thomas_constant2; 
		
		// Thomas solver coefficients 

		M.thomas_cx.assign( M.mesh.x_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomx.assign( M.mesh.x_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomx[0] = M.thomas_constant3a; 
		M.thomas_denomx[ M.mesh.x_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.x_coordinates.size() == 1 )
		{ M.thomas_denomx[0] = M.one; M.thomas_denomx[0] += M.thomas_constant2; } 

		M.thomas_cx[0] /= M.thomas_denomx[0]; 
		for( unsigned int i=1 ; i <= M.mesh.x_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomx[i] , M.thomas_constant1 , M.thomas_cx[i-1] ); 
			M.thomas_cx[i] /= M.thomas_denomx[i]; // the value at  size-1 is not actually used  
		}

		M.thomas_cy.assign( M.mesh.y_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomy.assign( M.mesh.y_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomy[0] = M.thomas_constant3a; 
		M.thomas_denomy[ M.mesh.y_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.y_coordinates.size() == 1 )
		{ M.thomas_denomy[0] = M.one; M.thomas_denomy[0] += M.thomas_constant2; } 

		M.thomas_cy[0] /= M.thomas_denomy[0]; 
		for( unsigned int i=1 ; i <= M.mesh.y_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomy[i] , M.thomas_constant1 , M.thomas_cy[i-1] ); 
			M.thomas_cy[i] /= M.thomas_denomy[i]; // the value at  size-1 is not actually used  
		}

		M.diffusion_solver_setup_done = true; 
	}

	// set the pointer
	
	M.apply_dirichlet_conditions();

	// x-diffusion 
	#pragma omp parallel for 
	for( unsigned int j=0; j < M.mesh.y_coordinates.size() ; j++ )
	{
		// Thomas solver, x-direction

		// remaining part of forward elimination, using pre-computed quantities 
		unsigned int n = M.voxel_index(0,j,0);
		(*M.p_density_vectors)[n] /= M.thomas_denomx[0]; 

		n += M.thomas_i_jump; 
		for( unsigned int i=1; i < M.mesh.x_coordinates.size() ; i++ )
		{
			axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_i_jump] ); 
			(*M.p_density_vectors)[n] /= M.thomas_denomx[i]; 
			n += M.thomas_i_jump; 
		}

		// back substitution 
		n = M.voxel_index( M.mesh.x_coordinates.size()-2 ,j,0); 

		for( int i = M.mesh.x_coordinates.size()-2 ; i >= 0 ; i-- )
		{
			naxpy( &(*M.p_density_vectors)[n] , M.thomas_cx[i] , (*M.p_density_vectors)[n+M.thomas_i_jump] ); 
			n -= M.thomas_i_jump; 
		}
	}

	// y-diffusion 

	M.apply_dirichlet_conditions();
	#pragma omp parallel for 
	for( unsigned int i=0; i < M.mesh.x_coordinates.size() ; i++ )
	{
		// Thomas solver, y-direction

		// remaining part of forward elimination, using pre-computed quantities 

		int n = M.voxel_index(i,0,0);
		(*M.p_density_vectors)[n] /= M.thomas_denomy[0]; 

		n += M.thomas_j_jump; 
		for( unsigned int j=1; j < M.mesh.y_coordinates.size() ; j++ )
		{
			axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_j_jump] ); 
			(*M.p_density_vectors)[n] /= M.thomas_denomy[j]; 
			n += M.thomas_j_jump; 
		}

		// back substitution 
		n = M.voxel_index( i,M.mesh.y_coordinates.size()-2, 0); 

		for( int j = M.mesh.y_coordinates.size()-2 ; j >= 0 ; j-- )
		{
			naxpy( &(*M.p_density_vectors)[n] , M.thomas_cy[j] , (*M.p_density_vectors)[n+M.thomas_j_jump] ); 
			n -= M.thomas_j_jump; 
		}
	}

	M.apply_dirichlet_conditions();
	
	// reset gradient vectors 
//	M.reset_all_gradient_vectors(); 
	
	return; */
}

void diffusion_decay_explicit_uniform_rates( Microenvironment& M, double dt )
{
    std::cout << "Diffusion decay explicit uniform rates is desactivated" << std::endl;
    /*
	using std::vector; 
	using std::cout; 
	using std::endl; 

	// static int n_jump_i = 1; 
	// static int n_jump_j = M.mesh.x_coordinates.size(); 
	// static int n_jump_k = M.mesh.x_coordinates.size() * M.mesh.y_coordinates.size(); 

	if( !M.diffusion_solver_setup_done )
	{	
		M.thomas_i_jump = 1; 
		M.thomas_j_jump = M.mesh.x_coordinates.size(); 
		M.thomas_k_jump = M.thomas_j_jump * M.mesh.y_coordinates.size(); 
	
		M.diffusion_solver_setup_done = true; 
	}
	
	if( M.mesh.uniform_mesh == false )
	{
		cout << "Error: This algorithm is written for uniform Cartesian meshes. Try: something else" << endl << endl; 
		return; 
	}

	// double buffering to reduce memory copy / allocation overhead 

	static vector< vector<double> >* pNew = &(M.temporary_density_vectors1);
	static vector< vector<double> >* pOld = &(M.temporary_density_vectors2);

	// swap the buffers 

	vector< vector<double> >* pTemp = pNew; 
	pNew = pOld; 
	pOld = pTemp; 
	M.p_density_vectors = pNew; 

	// static bool reaction_diffusion_shortcuts_are_set = false; 

	static vector<double> constant1 = (1.0 / ( M.mesh.dx * M.mesh.dx )) * M.diffusion_coefficients; 
	static vector<double> constant2 = dt * constant1; 
	static vector<double> constant3 = M.one + dt * M.decay_rates;

	static vector<double> constant4 = M.one - dt * M.decay_rates;

	#pragma omp parallel for
	for( unsigned int i=0; i < (*(M.p_density_vectors)).size() ; i++ )
	{
		unsigned int number_of_neighbors = M.mesh.connected_voxel_indices[i].size(); 

		double d1 = -1.0 * number_of_neighbors; 

		(*pNew)[i] = (*pOld)[i];  
		(*pNew)[i] *= constant4; 

		for( unsigned int j=0; j < number_of_neighbors ; j++ )
		{
			axpy( &(*pNew)[i], constant2, (*pOld)[  M.mesh.connected_voxel_indices[i][j] ] ); 
		}
		vector<double> temp = constant2; 
		temp *= d1; 
		axpy( &(*pNew)[i] , temp , (*pOld)[i] ); 
	}
	
	// reset gradient vectors 
//	M.reset_all_gradient_vectors(); 
*/
	return; 
}

void diffusion_decay_solver__constant_coefficients_LOD_1D( Microenvironment& M, double dt )
{
	if( M.mesh.regular_mesh == false )
	{
		std::cout << "Error: This algorithm is written for regular Cartesian meshes. Try: something else." << std::endl << std::endl; 
		return; 
	}
	
    std::cout << "Diffusion decay solver constant 1D is desactivated" << std::endl;

	// constants for the linear solver (Thomas algorithm) 
	/*
	if( !M.diffusion_solver_setup_done )
	{
		std::cout << std::endl << "Using method " << __FUNCTION__ << " (2D LOD with Thomas Algorithm) ... " << std::endl << std::endl;  
		
		M.thomas_denomx.resize( M.mesh.x_coordinates.size() , M.zero );
		M.thomas_cx.resize( M.mesh.x_coordinates.size() , M.zero );

		// define constants and pre-computed quantities 

		M.thomas_i_jump = 1; 
		M.thomas_j_jump = M.mesh.x_coordinates.size(); 

		M.thomas_constant1 =  M.diffusion_coefficients; //   dt*D/dx^2 
		M.thomas_constant1a = M.zero; // -dt*D/dx^2; 
		M.thomas_constant2 =  M.decay_rates; // (1/2)*dt*lambda 
		M.thomas_constant3 = M.one; // 1 + 2*constant1 + constant2; 
		M.thomas_constant3a = M.one; // 1 + constant1 + constant2; 
		
		M.thomas_constant1 *= dt; 
		M.thomas_constant1 /= M.mesh.dx; 
		M.thomas_constant1 /= M.mesh.dx; 

		M.thomas_constant1a = M.thomas_constant1; 
		M.thomas_constant1a *= -1.0; 

		M.thomas_constant2 *= dt; 
		M.thomas_constant2 *= 1; // no splitting via LOD

		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant2; 

		M.thomas_constant3a += M.thomas_constant1; 
		M.thomas_constant3a += M.thomas_constant2; 
		
		// Thomas solver coefficients 

		M.thomas_cx.assign( M.mesh.x_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomx.assign( M.mesh.x_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomx[0] = M.thomas_constant3a; 
		M.thomas_denomx[ M.mesh.x_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.x_coordinates.size() == 1 )
		{ M.thomas_denomx[0] = M.one; M.thomas_denomx[0] += M.thomas_constant2; } 

		M.thomas_cx[0] /= M.thomas_denomx[0]; 
		for( unsigned int i=1 ; i <= M.mesh.x_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomx[i] , M.thomas_constant1 , M.thomas_cx[i-1] ); 
			M.thomas_cx[i] /= M.thomas_denomx[i]; // the value at  size-1 is not actually used  
		}

		M.diffusion_solver_setup_done = true; 
	}

	// set the pointer
	
	M.apply_dirichlet_conditions();

	// x-diffusion 
	#pragma omp parallel for 
	for( unsigned int j=0; j < M.mesh.y_coordinates.size() ; j++ )
	{
		// Thomas solver, x-direction

		// remaining part of forward elimination, using pre-computed quantities 
		unsigned int n = M.voxel_index(0,j,0);
		(*M.p_density_vectors)[n] /= M.thomas_denomx[0]; 

		n += M.thomas_i_jump; 
		for( unsigned int i=1; i < M.mesh.x_coordinates.size() ; i++ )
		{
			axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_i_jump] ); 
			(*M.p_density_vectors)[n] /= M.thomas_denomx[i]; 
			n += M.thomas_i_jump; 
		}

		// back substitution 
		n = M.voxel_index( M.mesh.x_coordinates.size()-2 ,j,0); 

		for( int i = M.mesh.x_coordinates.size()-2 ; i >= 0 ; i-- )
		{
			naxpy( &(*M.p_density_vectors)[n] , M.thomas_cx[i] , (*M.p_density_vectors)[n+M.thomas_i_jump] ); 
			n -= M.thomas_i_jump; 
		}
	}

	M.apply_dirichlet_conditions();
	
	// reset gradient vectors 
//	M.reset_all_gradient_vectors(); 
	*/
	return; 
}


};
