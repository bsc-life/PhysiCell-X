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

/*--------------------------------------------------------------------------------------*/
/* Parallel version of the 3-D Thomas Solver, uses X-decomposition i.e. the X-direction */
/* of BioFVM going from left to right is divided into slices, each slice given to 1 MPI */
/* process. 
/*--------------------------------------------------------------------------------------*/

void diffusion_decay_solver__constant_coefficients_LOD_3D( Microenvironment& M, double dt, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	MPI_Request send_req, recv_req;
	
	double t_strt_set, t_end_set; 
	double t_strt_x, t_end_x;
	double t_strt_y,t_end_y;
	double t_strt_z,t_end_z; 
    
	
  if( M.mesh.uniform_mesh == false || M.mesh.Cartesian_mesh == false )
	{
		std::cout << "Error: This algorithm is written for uniform Cartesian meshes. Try: other solvers!" << std::endl << std::endl; 
		return; 
	}

/*----------------------------------------------*/
/* define constants and pre-computed quantities */
/*----------------------------------------------*/
	
	if( ! M.diffusion_solver_setup_done )
	{
		 if(world.rank == 0)
		  	std::cout << std::endl << "Using method " << __FUNCTION__ << " (implicit 3-D LOD with Thomas Algorithm) ... " 
		  	<< std::endl << std::endl;  
		
/*-------------------------------------------------------------------------------------------*/
/* x_coordinates are of size local_x_nodes,(see function resize() of class Cartesian_Mesh 	 */ 
/* in BioFVM_Cartesian.cpp. Each line of Voxels going from left to right forms a tridiagonal */ 
/* system of Equations. Now these lines are going to split in the X decomposition.				   */
/*-------------------------------------------------------------------------------------------*/
		
    M.thomas_denomx.resize( M.mesh.x_coordinates.size() , M.zero );           //sizeof(x_coordinates) = local_x_nodes, denomx is the main diagonal elements
		M.thomas_cx.resize( M.mesh.x_coordinates.size() , M.zero );               //Both b and c of tridiagonal matrix are equal, hence just one array needed
		
/*-----------------------------------------------------------------------*/
/* y_coordinates are of size of local_y_nodes. Each line of Voxels going */
/* from bottom to top forms a tridiagonal system of Equations.  				 */
/*-----------------------------------------------------------------------*/

		M.thomas_denomy.resize( M.mesh.y_coordinates.size() , M.zero );           
		M.thomas_cy.resize( M.mesh.y_coordinates.size() , M.zero );
        
/*-----------------------------------------------------------------------*/
/* z_coordinates are of size of local_z_nodes. Each line of Voxels going */         
/* from front to back forms a tridiagonal system of Equations. 	 				 */
/*-----------------------------------------------------------------------*/
		
		M.thomas_denomz.resize( M.mesh.z_coordinates.size() , M.zero );
		M.thomas_cz.resize( M.mesh.z_coordinates.size() , M.zero );

/*-------------------------------------------------------------*/
/* For X-decomposition thomas_i_jump - 1 can be in the previous*/
/* process and thomas_i_jump+1 can be in the next process      */
/* hence we can use thomas_j_jump and thomas_k_jump safely     */
/* but we CANNOT use thomas_i_jump safely                      */              
/*-------------------------------------------------------------*/

		M.thomas_i_jump = 1; 
		M.thomas_j_jump = M.mesh.x_coordinates.size(); 
		M.thomas_k_jump = M.thomas_j_jump * M.mesh.y_coordinates.size(); 
        
/*-------------------------------------------------------------*/
/* This part below of defining constants SHOULD typically      */
/* not change during parallelization.                          */
/*-------------------------------------------------------------*/

		M.thomas_constant1  = M.diffusion_coefficients;      // dt*D/dx^2 
		M.thomas_constant1a = M.zero;                        // -dt*D/dx^2; 
		M.thomas_constant2  = M.decay_rates;                 // (1/3)* dt*lambda 
		M.thomas_constant3  = M.one;                         // 1 + 2*constant1 + constant2; 
		M.thomas_constant3a = M.one;                         // 1 + constant1 + constant2; 		
			
		M.thomas_constant1 *= dt; 
		M.thomas_constant1 /= M.mesh.dx; 
		M.thomas_constant1 /= M.mesh.dx; 

		M.thomas_constant1a = M.thomas_constant1; 
		M.thomas_constant1a *= -1.0; 

		M.thomas_constant2 *= dt; 
		M.thomas_constant2 /= 3.0;                            // for the LOD splitting of the source, division by 3 is for 3-D 

		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant2; 

		M.thomas_constant3a += M.thomas_constant1; 
		M.thomas_constant3a += M.thomas_constant2; 

		// Thomas solver coefficients 
        
/*-----------------------------------------------------------------------*/
/* In 1-D X decomposition, x-lines are partitioned and hence assignments */
/* for the X-dimension WILL change. For y/z, they SHOULD not change.		 */
/*-----------------------------------------------------------------------*/

		M.thomas_cx.assign( M.mesh.x_coordinates.size() , M.thomas_constant1a );                  //Fill b and c elements with -D * dt/dx^2 
		M.thomas_denomx.assign( M.mesh.x_coordinates.size()  , M.thomas_constant3 );              //Fill diagonal elements with (1 + 1/3 * lambda * dt + 2*D*dt/dx^2)
		
		if(world.rank == 0)
            M.thomas_denomx[0] = M.thomas_constant3a;                                                 //First diagonal element is   (1 + 1/3 * lambda * dt + 1*D*dt/dx^2)
            
		if(world.rank == (world.size-1))
            M.thomas_denomx[ M.mesh.x_coordinates.size()-1 ] = M.thomas_constant3a;                   //Last diagonal element  is   (1 + 1/3 * lambda * dt + 1*D*dt/dx^2) 
		
		if(world.rank == 0)
       if( M.mesh.x_coordinates.size() == 1 )                                                    //This is an extreme case, won't exist, still if it does
        {                                                                                         //then this must be at rank 0
           M.thomas_denomx[0] = M.one; 
           M.thomas_denomx[0] += M.thomas_constant2; 
        } 
        
    if(world.rank == 0)
           M.thomas_cx[0] /= M.thomas_denomx[0];                                                     //The first c element of tridiagonal matrix is div by first diagonal el.

/*------------------------------------------------*/		
/* axpy(1st, 2nd, 3rd) => 1st = 1st + 2nd * 3rd 	*/
/* The value at  size-1 is not actually used  		*/
/* Since value of size-1 is not used, it means it */
/* is the value after the last Diagonal element   */
/* -----------------------------------------------*/
                                                                                                 
        for(int ser_ctr=0; ser_ctr<=world.size-1; ser_ctr++)
        {
            if(world.rank == ser_ctr)
            {
            		/* LATER CHECK IF CONDITIONS CAN BE MERGED - BECAUSE IT LOOKS LIKE A REPETITION */
            		
                if(world.rank == 0 && world.rank <= world.size-1)                   //If size=1, then this process does not send data
                {
                    
                    for( int i=1 ; i <= M.mesh.x_coordinates.size()-1 ; i++ )
                    { 
                        axpy( &M.thomas_denomx[i] , M.thomas_constant1 , M.thomas_cx[i-1] ); 
                        M.thomas_cx[i] /= M.thomas_denomx[i];                  // the value at  size-1 is not actually used  
                    }
                }
                else
                {
                    for( int i=1 ; i <= M.mesh.x_coordinates.size()-1 ; i++ )
                    { 
                        axpy( &M.thomas_denomx[i] , M.thomas_constant1 , M.thomas_cx[i-1] ); 
                        M.thomas_cx[i] /= M.thomas_denomx[i];                   // the value at  size-1 is not actually used  
                    }
                }
                    
                if(world.rank < (world.size-1))
                {
                    MPI_Isend(&(M.thomas_cx[M.mesh.x_coordinates.size()-1][0]), M.thomas_cx[M.mesh.x_coordinates.size()-1].size(), MPI_DOUBLE, ser_ctr+1, 1111, cart_topo.mpi_cart_comm, &send_req);
                    MPI_Wait(&send_req, MPI_STATUS_IGNORE);
                }
            }
                    
            if(world.rank == (ser_ctr+1) && (ser_ctr+1) <= (world.size-1))
            {
                    
                std::vector<double> temp_cx(M.thomas_cx[0].size());
                    
                MPI_Irecv(&temp_cx[0], temp_cx.size(), MPI_DOUBLE, ser_ctr, 1111, cart_topo.mpi_cart_comm, &recv_req);
                MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
                    
                   
                axpy( &M.thomas_denomx[0] , M.thomas_constant1 , temp_cx );                        //CHECK IF &temp_cz[0] is OK, axpy() in BioFVM_vector.cpp 
                M.thomas_cx[0] /= M.thomas_denomx[0];                                                  // the value at  size-1 is not actually used  
            }
            
            MPI_Barrier(cart_topo.mpi_cart_comm); 
        }
                                                                                                 
/*--------------------------------------------------------------------*/
/* In 1-D X decomposition, z and y-lines are contiguous and typically */
/* the assignments below for z, y should not be changed.              */
/* Both the first voxel i.e. index 0 and last voxel i.e. index=       */
/* y_coordinates.size()-1 are on the same process (same for Z)        */
/*--------------------------------------------------------------------*/
		
		M.thomas_cy.assign( M.mesh.y_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomy.assign( M.mesh.y_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomy[0] = M.thomas_constant3a; 
		M.thomas_denomy[ M.mesh.y_coordinates.size()-1 ] = M.thomas_constant3a;
		 
		if( M.mesh.y_coordinates.size() == 1 )
		{ 
        M.thomas_denomy[0] = M.one; 
        M.thomas_denomy[0] += M.thomas_constant2;        
    } 
    
		M.thomas_cy[0] /= M.thomas_denomy[0];
		 
		for( int i=1 ; i <= M.mesh.y_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomy[i] , M.thomas_constant1 , M.thomas_cy[i-1] ); 
			M.thomas_cy[i] /= M.thomas_denomy[i];                                         // the value at  size-1 is not actually used  
		}

/*----------------------------------*/
/* Setting Z-coordinate parameters  */
/*----------------------------------*/

		M.thomas_cz.assign( M.mesh.z_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomz.assign( M.mesh.z_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomz[0] = M.thomas_constant3a; 
		M.thomas_denomz[ M.mesh.z_coordinates.size()-1 ] = M.thomas_constant3a; 
		
		if( M.mesh.z_coordinates.size() == 1 )
		{ 
            M.thomas_denomz[0] = M.one; 
            M.thomas_denomz[0] += M.thomas_constant2; 
    }
     
		M.thomas_cz[0] /= M.thomas_denomz[0]; 
		
		for( int i=1 ; i <= M.mesh.z_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomz[i] , M.thomas_constant1 , M.thomas_cz[i-1] ); 
			M.thomas_cz[i] /= M.thomas_denomz[i];                                        // the value at  size-1 is not actually used  
		}	

		M.diffusion_solver_setup_done = true; 
	
	}

	// x-diffusion 
	
/*------------------------------------------------------------*/
/* CHECK LATER IF DIRICHLET CONDITIONS ARE BEING SET PROPERLY	*/
/*------------------------------------------------------------*/

	M.apply_dirichlet_conditions();																									
    
 /*-----------------------------------------------------------------------------------*/
 /*                        FORWARD ELIMINATION - x DIRECTION/DECOMPOSITION            */
 /*-----------------------------------------------------------------------------------*/
 
 /*----------------------------------------------------------------------------------------------------*/
 /* 2-D Y-Z planes are to be packed and sent to adjacent processes.                                    */
 /* My direction of traversing is go up up up i.e. y direction points then go in i.e. Z-direction      */
 /* Remember to visualize 3-D as 2-D plates kept after one another. Hence Z-direction data is farther  */
 /* apart than X/Y direction.                                                                          */
 /*----------------------------------------------------------------------------------------------------*/
 

 int y_size = M.mesh.y_coordinates.size();
 int z_size = M.mesh.z_coordinates.size();
 int p_size = (*M.p_density_vectors)[0].size();     //All p_density_vectors elements have same size, use anyone
 
 
 int snd_data_size = z_size * y_size * p_size;      //Number of data elements to be sent                                      
 int rcv_data_size = z_size * y_size * p_size;      //Number of data elements to be received
 
 
 std::vector<double> snd_data(snd_data_size);
 std::vector<double> rcv_data(rcv_data_size);

/*--------------------------------------------------------------------------------------------------*/ 
/* So row is along Z axis, column of each row is along Y-axis and each element has p_density_vector */
/*--------------------------------------------------------------------------------------------------*/ 
 
 std::vector<std::vector<std::vector<double>>> block3d(z_size, std::vector<std::vector<double>>( y_size, std::vector<double>(p_size)));
    
 //t_strt_x = MPI_Wtime();  

/*-----------------------------------------------------------------------------------------*/
/* The zeroth element on each process in the x-direction with rank >= 1 is to be treated	 */ 
/* separately because its previous element is on the process before it and so data must be */
/* caught in block3D and then calculations are to be done. 																 */
/* Data is sent using a 1-dimensional array, received in 1-dimensional array BUT then that */
/* data is transferred to a 3D array called 'block3d', its easier to use this block3d in 	 */
/* 'axpy' calculations. 																																	 */
/*-----------------------------------------------------------------------------------------*/	
	
for(int ser_ctr=0; ser_ctr <= world.size-1; ser_ctr++)
 {
     if(world.rank == ser_ctr)
     {
         if(world.rank == 0)
         {
            #pragma omp parallel for
            for(int k=0; k<= M.mesh.z_coordinates.size()-1; k++)
            {
                for(int j=0; j<= M.mesh.y_coordinates.size()-1; j++)
                {
                    int n = M.voxel_index(0,j,k);
                    (*M.p_density_vectors)[n] /= M.thomas_denomx[0];
                    
                    for( int i=1; i < M.mesh.x_coordinates.size() ; i++ )
                    {
                        int n  = M.voxel_index(i,j,k); 
                        //int n1 = M.voxel_index(i-1,j,k);            //Can remove this overhead of finding index now
                        axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_i_jump] ); 
                        (*M.p_density_vectors)[n] /= M.thomas_denomx[i]; 
                    }
                }
            }
         }
         else
         {
            #pragma omp parallel for 
            for(int k=0; k<= M.mesh.z_coordinates.size()-1; k++)
            {
                for(int j=0; j<= M.mesh.y_coordinates.size()-1; j++)
                {
                    int n = M.voxel_index(0,j,k);                               //Need to consider case separately for k=0, as k-1 would be -1 ! 
                    axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , block3d[k][j]); 
                    (*M.p_density_vectors)[n] /= M.thomas_denomx[0];
                    
                    for( int i=1; i < M.mesh.x_coordinates.size() ; i++ )
                    {
                        int n = M.voxel_index(i,j,k); 
                        //int n1 = M.voxel_index(i-1,j,k);
                        axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_i_jump] ); 
                        (*M.p_density_vectors)[n] /= M.thomas_denomx[i]; 
                    }
                }
            }
         }
         
         if(world.rank < (world.size-1))
         {
             int x_end = M.mesh.x_coordinates.size()-1; 
             int ctr=0;
             
             for(int k=0; k<= M.mesh.z_coordinates.size()-1; k++)
             {
                for(int j=0; j<= M.mesh.y_coordinates.size()-1; j++)
                {
                    int n = M.voxel_index(x_end,j,k); 
                    for(int ele=0; ele <= (*M.p_density_vectors)[n].size()-1; ele++)
                        snd_data[ctr++] = (*M.p_density_vectors)[n][ele]; 
                }
             }
             MPI_Isend(&snd_data[0], snd_data_size, MPI_DOUBLE, ser_ctr+1, 1111, cart_topo.mpi_cart_comm, &send_req);
             MPI_Wait(&send_req, MPI_STATUS_IGNORE); 
         }         
     }
     if(world.rank == (ser_ctr+1) && world.rank <= (world.size-1))
     {
             //Receive the data here and try to put in same format as vector of vectors in 'block3d'
             
             MPI_Irecv(&rcv_data[0], rcv_data_size, MPI_DOUBLE, ser_ctr, 1111,  cart_topo.mpi_cart_comm, &recv_req);
             MPI_Wait(&recv_req,MPI_STATUS_IGNORE);
              
             int ctr = 0;
             for(int m=0; m<z_size; m++)
             {
                 for(int n=0; n<y_size; n++)
                 {
                     for(int p=0; p<p_size; p++)
                     {
                         block3d[m][n][p]=rcv_data[ctr++];
                     }
                 }
             }
    }
    MPI_Barrier(cart_topo.mpi_cart_comm);
 }
 
 
 /*----------------------------------------------------------------------------------------*/
 /*  BACK SUBSITUTION for x-direction: while doing Forward elimination, we have gone like  */
 /* left to right, bottom to top then front to back - so we have reached the back-top-right*/ 						
 /* corner of the physical domain. Now we begin from the same point and going towards			 */
 /* front left bottom corner. This way there is a possibility that last used data may be 	 */
 /* found in the cache hierarchy. 																												 */
 /*----------------------------------------------------------------------------------------*/
 
 
 for(int ser_ctr = world.size-1; ser_ctr >= 0; ser_ctr--)
 {
     if(world.rank == ser_ctr)
     {
         if(world.rank == (world.size-1))
         {
            #pragma omp parallel for
            for(int k=M.mesh.z_coordinates.size()-1; k>=0 ; k--)
            {
                for(int j=M.mesh.y_coordinates.size()-1; j>=0 ; j--)
                {
                    for( int i=M.mesh.x_coordinates.size()-2; i>=0 ; i-- )
                    {
                        int n  = M.voxel_index(i,j,k); 
                        //int n2 = M.voxel_index(i+1,j,k);                                    //can remove overhead of finding index here
                        naxpy( &(*M.p_density_vectors)[n] , M.thomas_cx[i] , (*M.p_density_vectors)[n+M.thomas_i_jump] );
                    }
                }
            }
         }
         else
         {
            #pragma omp parallel for 
            for(int k=M.mesh.z_coordinates.size()-1; k>=0 ; k--)
            {
                for(int j=M.mesh.y_coordinates.size()-1; j>=0 ; j--)
                {
                    int n = M.voxel_index(M.mesh.x_coordinates.size()-1,j,k);                                
                    naxpy( &(*M.p_density_vectors)[n] , M.thomas_cx[M.mesh.x_coordinates.size()-1] , block3d[k][j]);
                    
                    for( int i=M.mesh.x_coordinates.size()-2; i >=0; i-- )
                    {
                        int n  = M.voxel_index(i,j,k); 
                        //int n2 = M.voxel_index(i+1,j,k);
                        naxpy( &(*M.p_density_vectors)[n] , M.thomas_cx[i] , (*M.p_density_vectors)[n+M.thomas_i_jump] );
                    }
                }
            }
         }
         if(world.rank > 0)
         {
             int x_start = 0; 
             int ctr=0;
             
             for(int k=0; k<= M.mesh.z_coordinates.size()-1; k++)
             {
                for(int j=0; j<= M.mesh.y_coordinates.size()-1; j++)
                {
                    int n = M.voxel_index(x_start,j,k); 
                    for(int ele=0; ele <= (*M.p_density_vectors)[n].size()-1; ele++)
                        snd_data[ctr++] = (*M.p_density_vectors)[n][ele]; 
                }
             }
             MPI_Isend(&snd_data[0], snd_data_size, MPI_DOUBLE, ser_ctr-1, 1111, cart_topo.mpi_cart_comm, &send_req);
             MPI_Wait(&send_req, MPI_STATUS_IGNORE); 
         }         
     }
     if(world.rank == (ser_ctr-1) && world.rank >= 0)
     {
             //Receive the data here and try to put in same format as vector of vectors
             MPI_Irecv(&rcv_data[0], rcv_data_size, MPI_DOUBLE, ser_ctr, 1111,  cart_topo.mpi_cart_comm, &recv_req);
             MPI_Wait(&recv_req,MPI_STATUS_IGNORE); 
             
             int ctr = 0;
             for(int m=0; m<z_size; m++)
             {
                 for(int n=0; n<y_size; n++)
                 {
                     for(int p=0; p<p_size; p++)
                     {
                         block3d[m][n][p]=rcv_data[ctr++];
                     }
                 }
             }
    }
    MPI_Barrier(cart_topo.mpi_cart_comm);
 }
  //t_end_x = MPI_Wtime();
  //std::cout<<"X solve time = "<<(t_end_x-t_strt_x)<<std::endl;
	
	/*-----------------------------------------------------*/
	/* y-diffusion : contiguous data i.e. no decomposition */
	/*-----------------------------------------------------*/

	M.apply_dirichlet_conditions();
	
	//t_strt_y = MPI_Wtime();
	#pragma omp parallel for 
	for( int k=0; k < M.mesh.z_coordinates.size() ; k++ )
	{
		for( int i=0; i < M.mesh.x_coordinates.size() ; i++ )
		{
			/*----------------------------------------------------------------------*/
			/* Thomas solver, y-direction : Forward Elimination Phase 							*/
			/* remaining part of forward elimination, using pre-computed quantities */
			/*----------------------------------------------------------------------*/

			int n = M.voxel_index(i,0,k);
			(*M.p_density_vectors)[n] /= M.thomas_denomy[0]; 

			for( int j=1; j < M.mesh.y_coordinates.size() ; j++ )
			{
				n = M.voxel_index(i,j,k); 
				axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_j_jump] ); 
				(*M.p_density_vectors)[n] /= M.thomas_denomy[j]; 
			}

			/*----------------------------------------------------------------------*/
			/* Thomas solver, y-direction : Backward Substitution  Phase 						*/
			/* n = voxel_index( mesh.x_coordinates.size()-2 ,j,k);									*/
			/*----------------------------------------------------------------------*/  

			for( int j = M.mesh.y_coordinates.size()-2 ; j >= 0 ; j-- )
			{
				n = M.voxel_index(i,j,k); 
				naxpy( &(*M.p_density_vectors)[n] , M.thomas_cy[j] , (*M.p_density_vectors)[n+M.thomas_j_jump] ); 
			}
  	}
 }
 //t_end_y = MPI_Wtime();
 //std::cout<<"Y solve time = "<<(t_end_y-t_strt_y)<<std::endl;

	/*-----------------------------------------------------*/
	/* z-diffusion : contiguous data i.e. no decomposition */
	/*-----------------------------------------------------*/

 M.apply_dirichlet_conditions();
 
    /*------------------------------------------------------------------------------*/
    /* NOTES: What if the Z-direction was decomposed ? Notes for a split z-direction*/
    /* NOTES REMAIN THE SAME FOR X-DECOMPOSITION																		*/
    /* PROCESSING OF 0TH ELEMENT IN ARRAY                                           */
    /* The processing on the 0th element of array will be done                      */
    /* only on rank = 0, or basically the processes containing                      */
    /* the front boundary of the domain.                                            */
    /* MPI divides domain into thick slices in Z-direction                          */
    /* OpenMP divides the thick slices in vertical direction i.e. horizontal lines. */
    /*------------------------------------------------------------------------------*/
 
    /*------------------------------------------------------------------------------*/
    /* FORWARD SUBSITUTION AND DATA SENDING OF BACK PLANE                           */
    /* Let all the processing finish at rank 0, now we need to send  the 'back'     */
    /* square of p_density_vectors[] to rank 1. Do the same at rank 2 so on...      */
    /* Possibly the thomas_k_jump will also change                                  */
    /* Instead of using thomas_k_jump, we can use prev=voxel_index(i,j,k-1)         */
    /* and then use it in (*M.p_density_vectors)[prev]                              */
    /* Also for rank > 0, it needs to start at k=0                                  */
    /* The values of p_density_vectors sent from the previous neighbour will be     */
    /* stored in temp_p_density_vectors -- only used for k=1 else just use          */
    /* p_density_vectors. Vector of vectors is not stored contiguously              */ 
    /* So we need to pack them into a single array and then send it to other end    */
    /*------------------------------------------------------------------------------*/
    
    /*---------------------------------------------------------------------------------*/
    /* BACK SUBSTITUTION AND DATA SENDING OF FRONT PLANE                               */
    /* This starts from the last process i.e. rank = size-1 then proceeds backwards    */
    /* i.e. from back boundary to front boundary. On the last process k starts from    */
    /* z_coordinates.size()-2 but on other processes, it starts from                   */
    /* z_coordinates.size()-1.                                                         */
    /*---------------------------------------------------------------------------------*/
 
 		/*------------------COMMENT OUT FROM HERE-----------------------------------------*/
 		/* So the code is like M.apply_dirichlet_conditions() then code for Z-diffusion   */
 		/* Then M.apply_dirichlet_conditions();                                           */
 		/*--------------------------------------------------------------------------------*/
 
		/*----------------------------------------------------------------------------------------------------------*/
		/*Declaring a matrix with z_size rows, y_size columns and each element of matrix is a vector of size p_size.*/
		/*First will store received data from rcv_data into this 3-d block because we need to pass base address     */
		/*of vectors in the routine axpy()                                                                          */
		/*----------------------------------------------------------------------------------------------------------*/

 
 		//t_strt_z = MPI_Wtime();
 
   #pragma omp parallel for 
   for( int j=0; j < M.mesh.y_coordinates.size() ; j++ )
  	{
   		for( int i=0; i < M.mesh.x_coordinates.size() ; i++ )
   		{
       /*----------------------------------------------------------------------*/
			 /* Thomas solver, Z-direction : Forward Elimination Phase 							 */
			 /* remaining part of forward elimination, using pre-computed quantities */
			 /*----------------------------------------------------------------------*/ 
     		
     		int n = M.voxel_index(i,j,0);
 				(*M.p_density_vectors)[n] /= M.thomas_denomz[0]; 
     
     		//should be an empty loop if mesh.z_coordinates.size() < 2
     	  
 				for( int k=1; k < M.mesh.z_coordinates.size() ; k++ )
 				{
 					n = M.voxel_index(i,j,k); 
 					axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_k_jump] ); 
 					(*M.p_density_vectors)[n] /= M.thomas_denomz[k]; 
 				}
 				
 			/*----------------------------------------------------------------------*/
			/* Thomas solver, z-direction : Backward Substitution  Phase 						*/
			/*----------------------------------------------------------------------*/
 				
 				for( int k = M.mesh.z_coordinates.size()-2 ; k >= 0 ; k-- )
 				{
 					n = M.voxel_index(i,j,k); 
 					naxpy( &(*M.p_density_vectors)[n] , M.thomas_cz[k] , (*M.p_density_vectors)[n+M.thomas_k_jump] ); 
 				}
   		}
   }
 
 	//t_end_z = MPI_Wtime();
 	//std::cout<<"Z solve time = "<<(t_end_z-t_strt_z)<<std::endl;
	
 M.apply_dirichlet_conditions();
	
	//reset gradient vectors : 				Commented out in original code 
	//M.reset_all_gradient_vectors(); Commented out in original code

	return; 
}


void diffusion_decay_solver__constant_coefficients_LOD_2D( Microenvironment& M, double dt )
{
	if( M.mesh.regular_mesh == false )
	{
		std::cout << "Error: This algorithm is written for regular Cartesian meshes. Try: something else." << std::endl << std::endl; 
		return; 
	}
	
	// constants for the linear solver (Thomas algorithm) 
	
	if( !M.diffusion_solver_setup_done )
	{
		std::cout << std::endl << "Using method " << __FUNCTION__ << " (2D LOD with Thomas Algorithm) ... " << std::endl << std::endl;  
		
		M.thomas_denomx.resize( M.mesh.x_coordinates.size() , M.zero );
		M.thomas_cx.resize( M.mesh.x_coordinates.size() , M.zero );

		M.thomas_denomy.resize( M.mesh.y_coordinates.size() , M.zero );
		M.thomas_cy.resize( M.mesh.y_coordinates.size() , M.zero );
		
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
	
	return; 
}

void diffusion_decay_explicit_uniform_rates( Microenvironment& M, double dt )
{
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

	return; 
}

};
