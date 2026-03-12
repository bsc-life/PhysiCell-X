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

#ifndef __BioFVM_mesh_h__
#define __BioFVM_mesh_h__

#include <iostream>
#include <vector> 

#include "BioFVM_matlab.h"
#include "../DistPhy/DistPhy_Environment.h"
#include "../DistPhy/DistPhy_Cartesian.h" 

using namespace DistPhy::mpi; 

namespace BioFVM{

 /*! \brief Voxels are the basic spatial container for densities, which are networked into meshes. 
  * 
  * Voxels are the basic spatial container for a finite volume method. Voxels are connected to 
  * other voxels into a General_Mesh, here most likely a Cartesian_Mesh. Voxel boundaries are Voxel_Faces. 
  * 
  * A Microenvironment Domain will include a network of Voxels (and Voxel_Faces), 
  * with a vector<double> of densities for each Voxel, along with rate constants, etc. 
  * The Domain may also include a vector<double> of flux coefficients for each Voxel_Face. 
 */
 
class Voxel
{

 private:
	friend std::ostream& operator<<(std::ostream& os, const Voxel& v); 
	/*!< outputs the Voxel to an open ostream 
	 * \param os -- the stream 
	 * \param mv -- the voxel you use this friendly friend operator on
	 * Example: Voxel v; 
	 *          cout << v << endl; 
	*/ 

 public:
	Voxel(); 
	int mesh_index; /*!< voxel's index in a General_Mesh :---> In parallel this is local voxel index (Gaurav Saxena)*/ 

	/*-----------------------------------------------------------------------------------*/
	/* Adding an index global_mesh_index to signify a global voxel index (Gaurav Saxena) */
  /*-----------------------------------------------------------------------------------*/
    
    int global_mesh_index; 
	
	double volume; /*!< voxel's volume (cubic spatial units) */ 
	std::vector<double> center; /*!< center of volume */
	bool is_Dirichlet;
	void stream_output_with_units( std::ostream& os , std::string units ) const;
};

class Voxel_Face
{
 private:
	friend std::ostream& operator<<(std::ostream& os , const Voxel_Face& vf ); 
	
 public:
	Voxel_Face(); 
	int mesh_index; 
	
	double surface_area; 
	std::vector<double> center; 
	std::vector<double> outward_normal; 
	std::vector<double> inward_normal; 
	
	void stream_output_with_units( std::ostream& os , std::string units ) const;
};

class General_Mesh
{
 private: 
	friend std::ostream& operator<<(std::ostream& os, const General_Mesh& mesh);  
	
	// this stores the indexing of the voxel faces (connect voxel i to voxel j, face stored at k)
	// only for use in a future release
	// std::unordered_map< int,std::unordered_map<int,int> > voxel_face_index_mapping; 
	
 public:
	General_Mesh();  
	
	// [xmin ymin zmin xmax ymax zmax ]
	std::vector<double> bounding_box; 
	
	std::vector<double> local_bounding_box; //Sub-domain boundaries
	
	std::vector<Voxel> voxels; 
	std::vector<Voxel_Face> voxel_faces; 
	// each voxel[k] has a list of connected voxels -- helpful for some numerical methods 
	std::vector< std::vector<int> > connected_voxel_indices;
    
  /*-----------------------------------------------------------------------------*/
  /* Need to specify global indices connected_voxel_global_indices (Gaurav Saxena) */
  /*-------------------------------------------------------------------------------*/
    std::vector< std::vector<int> > connected_voxel_global_indices; //(added by Gaurav Saxena)
	
	int nearest_voxel_index( std::vector<double>& position );

	   
	bool is_position_valid(double x, double y, double z);
	
	/*=======================================================================================*/
	/* Parallel prototype of two functions which determine if cell has crossed sub-domain		 */
	/*=======================================================================================*/
	
	bool has_it_crossed_to_left_subdomain (double x, mpi_Environment &world);
	bool has_it_crossed_to_right_subdomain(double x, mpi_Environment &world);
	
	/*=====================================================================================*/
	/* Parallel prototype of a function which handles cell positions at the sub-domain		 */
	/* boundaries. Remember it can change the x-coordinates in the "position" vector.			 */
	/*=====================================================================================*/
	
    void correct_position_within_subdomain(std::vector<double> &pos, mpi_Environment &world, mpi_Cartesian &cart_topo);
  
    /* In parallel scenario, return the local_bounding_box vector */
  
    std::vector<double> get_subdomain_limits();

	/* the following help manage the voxel faces */ 
	// returns the index of the voxel face connecting from voxels[i] to voxels[j] 
	int voxel_face_index( int i, int j ); 
	
	// returns the Voxel_Face connecting voxels[i] to voxels[j] 
	Voxel_Face& voxel_face(int i, int j );   
	// returns the normal vector from voxels[i] to voxels[j] 
	std::vector<double>& outward_normal( int i, int j ); 
	
	/*! This creates a Voxel_Face from voxels[i] to voxels[j], and another from voxels[j] to 
	    voxels[i], both with surface area SA. It also auto-updates connected_voxel_indices[i] 
	    and connected_voxel_indices[j]. */ 
	void connect_voxels(int i,int j, double SA);   
	
	void connect_voxels_faces_only(int i,int j, double SA); 
	void connect_voxels_indices_only(int i,int j, double SA);
    
    /*-----------------------------------------------------------------------------*/
    /* For parallel implementation, need to maintain global indices adjacency list */
    /*-----------------------------------------------------------------------------*/
    
    void connect_voxels_global_indices_only(int i,int j, double SA); //---> to maintain a list of connected global indices (added Gaurav Saxena)
	
	/*! This removes all connections between voxels[i] and voxels[j], and deletes the associated 
	    Voxel_Face(s). */
	void disconnect_voxels(int i, int j); 
	void clear_voxel_face_index_mapping( void );  
	
	bool Cartesian_mesh; 
	bool uniform_mesh; 
	bool regular_mesh;
	bool use_voxel_faces; 
	
	std::string units; 
	
	void display_information( std::ostream& os); 
	
	void write_to_matlab( std::string filename ); 
	
	/*==========================================*/
	/* Parallel prototype of the function above */	
	/*==========================================*/
	void write_to_matlab( std::string filename, mpi_Environment &world, mpi_Cartesian &cart_topo);

	void read_from_matlab( std::string filename ); 
};

class Cartesian_Mesh : public General_Mesh
{
 private:
 
 public:
	std::vector<double> x_coordinates; 
	std::vector<double> y_coordinates;
	std::vector<double> z_coordinates; 	
	std::vector< std::vector<int> > moore_connected_voxel_indices; // Keeps the list of voxels in the Moore nighborhood
	
	/*============================================================================*/
	/* Adding 2 moore lists that will contain global indexes of voxels 						*/
	/* one for left and one for right sub-domain boundary.												*/
	/*============================================================================*/
	
	std::vector< std::vector<int> > moore_connected_voxel_global_indices_left;
	std::vector< std::vector<int> > moore_connected_voxel_global_indices_right;
	 
	void create_moore_neighborhood(void);
	
	/*============================================================================*/
	/* Parallel version of create_mmore_neighbourhood() function above. 					*/
	/* Builds a local list of voxel indices and 2 lists of global indices.				*/
	/*============================================================================*/
	
	void create_moore_neighborhood(mpi_Environment &world, mpi_Cartesian &cart_topo);
	
	unsigned int voxel_index( unsigned int i, unsigned int j, unsigned int k ); 
	std::vector<unsigned int> cartesian_indices( unsigned int n ); 
	
	double dx;
	double dy;
	double dz; 
	
	double dV; 	
	double dS;

	double dS_xy;
	double dS_yz; 
	double dS_xz;

	int x_size;
	int y_size;
	int z_size;

	int n_substrates;

	double local_x_start;
	double local_y_start;
	double local_z_start;
	
	Cartesian_Mesh(); // done 
	
	Cartesian_Mesh( int , int , int );  
	
	void create_voxel_faces( void ); 

	void resize( int,int,int ); 
	void resize( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , int x_nodes, int y_nodes, int z_nodes ); 
	void resize( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx, double dy, double dz );
    
    /*======================================*/
    /* Parallel prototype of function above */
    /*======================================*/
    
    void resize( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx_new, double dy_new , double dz_new, mpi_Environment &world, mpi_Cartesian &cart_topo);
    
	void resize_uniform( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx ); 
	
	int nearest_voxel_index( std::vector<double>& position ); 
	
	int nearest_lcl_voxel_index( std::vector<double>& position );
	
	/*-----------------------------------------------------------------------------------*/
  /* Parallel prototype of a new function which returns a process specific local index */
  /* in which the Basic_Agent resides. 																								 */
  /*-----------------------------------------------------------------------------------*/
	int nearest_voxel_local_index( std::vector<double>& position, mpi_Environment &world, mpi_Cartesian &cart_topo );
	  
	int nearest_voxel_face_index( std::vector<double>& position );  
	std::vector<unsigned int> nearest_cartesian_indices( std::vector<double>& position ); 
	Voxel& nearest_voxel( std::vector<double>& position ); 
	
	void display_information( std::ostream& os ); 
	
	void read_from_matlab( std::string filename ); 
};

class Voronoi_Mesh : public General_Mesh
{
 private:
 
 public:
	void display_information( std::ostream& os); 
};

};

#endif
