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

#include "BioFVM_vector.h" 
#include "BioFVM_mesh.h" 
#include "../DistPhy/DistPhy_Environment.h"
#include "../DistPhy/DistPhy_Cartesian.h"
#include <algorithm>

using namespace DistPhy::mpi; 

namespace BioFVM{
	
/* voxels */
const int mesh_min_x_index=0;
const int mesh_min_y_index=1;
const int mesh_min_z_index=2;
const int mesh_max_x_index=3;
const int mesh_max_y_index=4;
const int mesh_max_z_index=5;	
Voxel::Voxel()
{
	mesh_index = 0; 
	volume = 10*10*10; 
	center.assign( 3 , 0.0 ); 
	is_Dirichlet = false;
}

std::ostream& operator<<(std::ostream& os, const Voxel& v)  
{
	static std::string tabbing = "\t\t\t\t"; 
	static std::string tabbing2 = "\t\t\t\t\t"; 
	os	<< tabbing << "<voxel ID=\"" << v.mesh_index << "\">"  << std::endl
		<< tabbing2 << "<center " << v.center << " />" << std::endl  
		<< tabbing2 << "<volume>" << v.volume << "</volume>" << std::endl  
		<< tabbing  << "</voxel>"; 

 return os; 
}

void Voxel::stream_output_with_units( std::ostream& os , std::string units ) const 
{
	static std::string tabbing = "\t\t\t\t"; 
	static std::string tabbing2 = "\t\t\t\t\t"; 
	os	<< tabbing << "<voxel ID=\"" << mesh_index << "\">"  << std::endl
		<< tabbing2 << "<center " << center << " units=\"" << units << "\" />" << std::endl 
		<< tabbing2 << "<volume units=\"cubic " << units << "\">" << volume << "</volume>" << std::endl
		<< tabbing  << "</voxel>"; 
	return; 
}

/* voxel faces */ 

Voxel_Face::Voxel_Face()
{
	mesh_index = 0;  
	
	surface_area = 10*10; 
	center.assign( 3 , 0.0  ); 
	outward_normal.assign( 3 , 0.0 ); 
	inward_normal.assign( 3 , 0.0 ); 
}

std::ostream& operator<<(std::ostream& os, const Voxel_Face& vf)  
{
	static std::string tabbing = "\t\t\t\t"; 
	static std::string tabbing2 = "\t\t\t\t\t"; 
	os	<< tabbing << "<voxel_face ID=\"" << vf.mesh_index << "\">"  << std::endl
		<< tabbing2 << "<center " << vf.center << " />" << std::endl  
		<< tabbing2 << "<outward_normal " << vf.outward_normal << " />" << std::endl  
		<< tabbing2 << "<inward_normal " << vf.inward_normal << " />" << std::endl  
		<< tabbing2 << "<surface_area>" << vf.surface_area << "</surface_area>" << std::endl  
		<< tabbing  << "</voxel_face>"; 

 return os; 
}

void Voxel_Face::stream_output_with_units( std::ostream& os , std::string units ) const
{
	static std::string tabbing = "\t\t\t\t"; 
	static std::string tabbing2 = "\t\t\t\t\t"; 
	os	<< tabbing << "<voxel_face ID=\"" << mesh_index << "\">"  << std::endl
		<< tabbing2 << "<center units=\"" << units << "\" " << center << " />" << std::endl  
		<< tabbing2 << "<outward_normal units=\"" << units << "\" " << outward_normal << " />" << std::endl  
		<< tabbing2 << "<inward_normal units=\"" << units << "\" " << inward_normal << " />" << std::endl  
		<< tabbing2 << "<surface_area units=\"square " << units << "\">" << surface_area << "</surface_area>" << std::endl  
		<< tabbing  << "</voxel_face>"; 
}

/* general meshes */ 

General_Mesh::General_Mesh()
{
	// x1, x2, y1, y2, z1, z2 
	bounding_box.assign(6,0.0); 
	bounding_box[mesh_min_x_index] = -0.5; 
	bounding_box[mesh_min_y_index] = -0.5; 
	bounding_box[mesh_min_z_index] = -0.5; 
	bounding_box[mesh_max_x_index] = 0.5; 
	bounding_box[mesh_max_y_index] = 0.5; 
	bounding_box[mesh_max_z_index] = 0.5; 
	
	/*----------------------------------------------------------------------*/
	/* Space allocated for a vector of doubles called 'local_bounding_box', */
	/* added by Gaurav Saxena, to contain local sub-domain boundaries				*/
	/*----------------------------------------------------------------------*/
	
	local_bounding_box.assign(6,0.0); 
	
	voxels.resize( 1 );  
	voxel_faces.resize( 0 ); 
	
	connected_voxel_indices.resize( 1 ); 
	connected_voxel_indices[0].clear(); 

	// voxel_face_index_mapping.clear();
	
	Cartesian_mesh = false; 
	uniform_mesh = false; 
	regular_mesh = false; 
	use_voxel_faces = true; 	
	units = "none"; 
}
 
std::ostream& operator<<(std::ostream& os, const General_Mesh& mesh)  
{
	std::boolalpha( os ); 
	static std::string tabbing = "\t"; 
	static std::string tabbing2 = "\t\t"; 
	static std::string tabbing3 = "\t\t\t"; 
	os	<< tabbing << "<mesh type=\"general\" uniform=\"" << mesh.uniform_mesh << "\" regular=\"" << mesh.regular_mesh << "\" units=\"" << mesh.units << "\">"  << std::endl
		<< tabbing2 << "<voxels>" << std::endl;
	for( unsigned int i=0; i < mesh.voxels.size() ; i++ )
	{ os << mesh.voxels[i] << std::endl; } 
	os 	<< tabbing2 << "</voxels>" << std::endl 
		<< tabbing2 << "<voxel_faces>" << std::endl;
	for( unsigned int i=0; i < mesh.voxel_faces.size() ; i++ )
	{ os << mesh.voxel_faces[i] << std::endl; } 
	os 	<< tabbing2 << "</voxel_faces>" << std::endl; 
	
	os	<< tabbing2 << "<voxel_connections>" << std::endl;
	for( unsigned int i=0 ; i < mesh.connected_voxel_indices.size() ; i++ )
	{
		os << tabbing3 << "<connected_voxel_indices ID=\"" << i << "\">" << std::endl; 
		for( unsigned int j=0; j < mesh.connected_voxel_indices[i].size() ; j++ )
		{
			os 	<< tabbing3 << "\t<index>" << (mesh.connected_voxel_indices[i])[j] << "</index>" << std::endl; 
		}
		os << tabbing3 << "</connected_voxel_indices>" << std::endl; 
	}
	
	os << tabbing2 << "</voxel_connections>" << std::endl; 
	
	
	os	<< tabbing  << "</mesh>"; 
		
		// later: output information on connected faces 
		
		// later: 

 return os; 
}

/*--------------------------------------------------------------------------------------*/
/* Two parallel functions to determine if the cell has crossed to the left sub-domain or*/
/* the right sub-domain. 																																*/
/*--------------------------------------------------------------------------------------*/

bool General_Mesh::has_it_crossed_to_left_subdomain (double x, mpi_Environment &world)
{
	if(world.rank != 0)
		if(x < local_bounding_box[mesh_min_x_index])
			return true; 
		else
			return false;
	else
		return false; 		
}

bool General_Mesh::has_it_crossed_to_right_subdomain (double x, mpi_Environment &world)
{
	if(world.rank != (world.size-1))
		if(x > local_bounding_box[mesh_max_x_index])
			return true; 
		else
			return false;
	else
		return false; 		
}

std::vector<double> General_Mesh::get_subdomain_limits()
{
	return local_bounding_box; 
}

bool General_Mesh::is_position_valid(double x, double y, double z)
{
	if(x< bounding_box[mesh_min_x_index] || x>bounding_box[mesh_max_x_index])
		return false;
	if(y< bounding_box[mesh_min_y_index] || y>bounding_box[mesh_max_y_index])
		return false;
	if(z< bounding_box[mesh_min_z_index] || z>bounding_box[mesh_max_z_index])
		return false;
	return true;
}

/*------------------------------------------------------------------------*/
/* A parallel function which determines if the position of the cell within*/
/* a sub-domain is valid or not. Can be used anytime we have a cell at		*/
/* the boundary of a sub-domain (i) if at left boundary then shift right  */
/* by +eps and (ii) if at right boundary then shift left by -eps.					*/
/* No need to examine the left and right boundary of process 0 and process*/
/* n-1, respectively, as if a cell is outside these, then it is marked 		*/
/* as a non-functional cell.																							*/
/*------------------------------------------------------------------------*/

void General_Mesh::correct_position_within_subdomain(std::vector<double> &pos, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	/*---------------------------------------------------*/
	/* Only for 1-dimensional pure x-decomposition 			 */
	/*---------------------------------------------------*/
	
	double eps2 = 2 * 1e-5; //1e-16 <--- this is the previous value, possibly causing errors like : x - eps2 = x
	
	/*----------------------------------------------------------------------------*/
	/* None of this is needed if there exists only 1 process i.e. world.size == 1 */
	/* If world.size == 1, then this code should not be executed otherwise it 		*/
	/* will pull back the cells from outside the domain into the domain.					*/
	/* IMPORTANT: This function pulls back the parent after cell division IF it 	*/
	/* has crossed the sub-domain boundary (in order to preserve center of mass)  */ 
	/* This function should NOT be called from the update_position() function.		*/
	/* TBD: If both daughter and parent cross sub-domain boundary, then the parent*/
	/* can be pulled back by 2 * eps units instead of just 1 * eps units - just 	*/
	/* to make the cell x-coordinate different (Later ask Miguel/Arnau to give a	*/
	/* better value or think of a strategy to remove it. 													*/
	/*----------------------------------------------------------------------------*/
	
	if(world.size > 1)
	{
		if(world.rank == 0)
		{
			if(pos[0] >= local_bounding_box[mesh_max_x_index])					//If cell position is right of right sub-boundary
			 	 pos[0] =  local_bounding_box[mesh_max_x_index] - eps2; 		//pull cell a bit left
		}
	
		if(world.rank == (world.size-1))
		{
			if(pos[0] <= local_bounding_box[mesh_min_x_index])					//If cell position is left of left sub-boundary
			   pos[0] =  local_bounding_box[mesh_min_x_index] + eps2; 			//pull cell a bit right
		}
	
		if(world.rank > 0 && world.rank < (world.size-1))							//Do both the above for processes 1 to (n-1)
		{
			if(pos[0] <= local_bounding_box[mesh_min_x_index])
			 	 pos[0] =  local_bounding_box[mesh_min_x_index] + eps2;
			 
			if(pos[0] >= local_bounding_box[mesh_max_x_index])
			 	 pos[0] =  local_bounding_box[mesh_max_x_index] - eps2;				
		}	
	}
}


void General_Mesh::connect_voxels_faces_only(int i,int j, double SA)
{
	// check to see if the voxels are connected -- implement later! 
	// Jose: to do this, (i,j) ids required within Voxel_Face class

	// create a new Voxel_Face connecting i to j

	Voxel_Face VF1; 
	int k = voxel_faces.size(); 
	VF1.mesh_index = k; 
	VF1.surface_area = SA; 
	VF1.outward_normal = voxels[j].center - voxels[i].center ; 
	normalize( &(VF1.outward_normal) );
	VF1.inward_normal = VF1.outward_normal; 
	VF1.inward_normal *= -1.0; 
	
	// convention: face is oriented from lower index to higher index 
	if( j < i )
	{ VF1.outward_normal *= -1.0; VF1.inward_normal *= -1.0; }

	// add it to the vector of voxel faces 
	
	voxel_faces.push_back( VF1 ); 

	return; 
}

void General_Mesh::connect_voxels_indices_only(int i,int j, double SA)
{
	// add j to the list of connected voxels for voxel i 
	
	//connected_voxel_indices[i].push_back( j ); 
	//connected_voxel_indices[j].push_back( i ); 
	if (std::find(connected_voxel_indices[i].begin(), connected_voxel_indices[i].end(), j) == connected_voxel_indices[i].end())
    {
        connected_voxel_indices[i].push_back(j);
    }

    // Check if voxel i is already in the list of connected voxels for voxel j
    if (std::find(connected_voxel_indices[j].begin(), connected_voxel_indices[j].end(), i) == connected_voxel_indices[j].end())
    {
        connected_voxel_indices[j].push_back(i);
    }
	return;
}

/*--------------------------------------------------------------------------------------*/
/* connect_voxels_global_indices_only() is a function added for parallel implementation */
/* When it connects two voxels, it uss the global indexes (added Gaurav Saxena)         */
/*---------------------------------------------------------------------------------------*/

void General_Mesh::connect_voxels_global_indices_only(int i,int j, double SA) // done 
{
	
	if (std::find(connected_voxel_global_indices[i].begin(), connected_voxel_global_indices[i].end(), j) == connected_voxel_global_indices[i].end())
    {
        connected_voxel_global_indices[i].push_back(j);
    }

    // Check if voxel i is already in the list of connected voxels for voxel j
    if (std::find(connected_voxel_global_indices[j].begin(), connected_voxel_global_indices[j].end(), i) == connected_voxel_global_indices[j].end())
    {
        connected_voxel_global_indices[j].push_back(i);
    }
	return;
}

void General_Mesh::connect_voxels(int i,int j, double SA)
{
	// check to see if the voxels are connected -- implement later!
	// Jose: to do this, (i,j) ids required within Voxel_Face class


	// create a new Voxel_Face connecting i to j

	Voxel_Face VF1; 
	int k = voxel_faces.size(); 
	VF1.mesh_index = k; 
	VF1.surface_area = SA; 
	VF1.outward_normal = voxels[j].center - voxels[i].center ; 
	normalize( &(VF1.outward_normal) );
	VF1.inward_normal = VF1.outward_normal; 
	VF1.inward_normal *= -1.0; 
	
	// convention: face is oriented from lower index to higher index 
	if( j < i )
	{ VF1.outward_normal *= -1.0; VF1.inward_normal *= -1.0; }

	// add it to the vector of voxel faces 
	
	voxel_faces.push_back( VF1 ); 

	return; 
}

void General_Mesh::display_information( std::ostream& os )
{
	os << std::endl << "Mesh information: " << std::endl 
	<< "type: general mesh" << std::endl 
	<< "Domain: " 
	<< "[" << bounding_box[0] << "," << bounding_box[3] << "] " <<  units << " x " 
	<< "[" << bounding_box[1] << "," << bounding_box[4] << "] " <<  units << " x " 
	<< "[" << bounding_box[2] << "," << bounding_box[5] << "] " <<  units << std::endl
	<< "   voxels: " << voxels.size() << std::endl
	<< "   voxel faces: " << voxel_faces.size() << std::endl
	<< "   volume: "; 

	double total_volume = 0.0; 
	for( unsigned int i=0; i < voxels.size(); i++ )
	{ total_volume += voxels[i].volume; }
	os << total_volume << " cubic " << units << std::endl; 
 
	return; 
}

void General_Mesh::write_to_matlab( std::string filename )
{ 
	unsigned int number_of_data_entries = voxels.size();
	unsigned int size_of_each_datum = 3 + 1; // x,y,z, volume 

	FILE* fp = write_matlab_header( size_of_each_datum, number_of_data_entries,  filename, "mesh" );  

	// storing data as cols 
	for( unsigned int i=0; i < number_of_data_entries ; i++ )
	{
		fwrite( (char*) &( voxels[i].center[0] ) , sizeof(double) , 1 , fp ); 
		fwrite( (char*) &( voxels[i].center[1] ) , sizeof(double) , 1 , fp ); 
		fwrite( (char*) &( voxels[i].center[2] ) , sizeof(double) , 1 , fp ); 
		fwrite( (char*) &( voxels[i].volume ) , sizeof(double) , 1 , fp ); 
	}

	fclose( fp ); 
}

/*----------------------------------------------------------------------------------------------*/
/* Parallel equivalent of the function above - remember now BioFVM has actually 3 MATLAB file 	*/
/* writing functions/snippets: (1) Writing Mesh (2) Writing Microenvironment (3) Writing agents */
/*----------------------------------------------------------------------------------------------*/

void General_Mesh::write_to_matlab( std::string filename, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	
    MPI_File fh;
    MPI_Offset file_size, offset; 
    MPI_Datatype etype, filetype; 
    double *buffer;                         //Will contain center[0/1/2] and volume in a contiguous buffer
    char char_filename[filename.size()+1];
    int elements_to_write; 
    
    /*----------------------------------------------------------------------------------------*/
    /* Now total data entries is now the sum of all entries on all processes                  */
    /* size of datum remains the same                                                         */
    /*----------------------------------------------------------------------------------------*/

    int number_of_data_entries = voxels.size();
	int size_of_each_datum = 3 + 1 ;		//density values are NOT to be taken into account here (see function above) 

    //Possibly we do not need to return anything over here, we can write a separate file at Master
    //All processes call this function - because William Groppe says MPI_File_open is collective operation
    
    /*----------------------------------------------------------------------------------------------------*/
    /* IMPORTANT: A problem will arise when writing the basic_agents header. Each process contains a 			*/
    /* different number of basic_agents and hence multiplying world.size * number_of_data_entries 	 			*/
    /* inside write_matlab4_header() would be INCORRECT. 																						 			*/
    /* SOLUTION: Instead of SENDING number_of_data_entries to write_matlab_header(), send world.size 			*/
    /* * number_of_data_entries to this function (for Mesh and Microenvironment MATLAB functions)		 			*/
    /* When dealing with basic_agents, do an MPI_Reduce() at the root process into number_of_data_entries */
    /* REMOVE the multiplication of world.size * number_of_data_entries from write_matlab4_header()				*/
    /*----------------------------------------------------------------------------------------------------*/

    write_matlab_header( size_of_each_datum, world.size * number_of_data_entries,  filename, "mesh", world, cart_topo );
    
    /*-----------------------------------------------------------*/
    /* IMPORTANT: POSSIBLY TRY TO REMOVE FOR PERFORMANCE REASONS */
    /*-----------------------------------------------------------*/
    
    MPI_Barrier(cart_topo.mpi_cart_comm); 
       
		// storing data as cols 
    buffer = new double[number_of_data_entries * size_of_each_datum];
    
    //std::cout<<"CX	"<<"CY	"<<"CZ	"<<"Vol	"<<"Density	\n"; 
    
    int n = 0;
		for( int i=0; i < number_of_data_entries ; i++ )
		{
      buffer[n++] = voxels[i].center[0];
      //std::cout<<mesh.voxels[i].center[0]<<"	";
      buffer[n++] = voxels[i].center[1];
      //std::cout<<mesh.voxels[i].center[1]<<"	";
      buffer[n++] = voxels[i].center[2];
      //std::cout<<mesh.voxels[i].center[2]<<"	";
      buffer[n++] = voxels[i].volume;
      //std::cout<<mesh.voxels[i].volume<<"	";
            
    	//fwrite( (char*) &( mesh.voxels[i].center[0] ) , sizeof(double) , 1 , fp ); 
			//fwrite( (char*) &( mesh.voxels[i].center[1] ) , sizeof(double) , 1 , fp );   
			//fwrite( (char*) &( mesh.voxels[i].center[2] ) , sizeof(double) , 1 , fp ); 
			//fwrite( (char*) &( mesh.voxels[i].volume ) , sizeof(double) , 1 , fp ); 
		}
	
		strcpy(char_filename, filename.c_str());
    
		MPI_File_open(cart_topo.mpi_cart_comm, char_filename, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);      //This file is already created while writing Matlab header
    MPI_File_get_size(fh,&file_size);
    
    offset = file_size + world.rank * sizeof(double) * number_of_data_entries * size_of_each_datum;
    etype = MPI_DOUBLE;
    filetype = MPI_DOUBLE; 
    elements_to_write = number_of_data_entries * size_of_each_datum; 
    
    MPI_File_set_view(fh, offset, etype, filetype, "native", MPI_INFO_NULL); 
    MPI_File_write(fh, buffer, elements_to_write, MPI_DOUBLE, MPI_STATUS_IGNORE);

     
		MPI_File_close(&fh);
    delete [] buffer;
    
		return;
}
 

void General_Mesh::read_from_matlab( std::string filename )
{
	unsigned int size_of_each_datum; 
	unsigned int number_of_data_entries; 
	FILE* fp = read_matlab_header( &size_of_each_datum, &number_of_data_entries,  filename ); 

	voxel_faces.resize( 0 ); 
	
	connected_voxel_indices.resize( 1 ); 
	connected_voxel_indices[0].clear(); 
	
	Cartesian_mesh = false; 
	uniform_mesh = false; 
	regular_mesh = false; 
	use_voxel_faces = false; 

	// resize the internal data structure 

	voxels.resize( number_of_data_entries );
	connected_voxel_indices.resize( voxels.size() ); 
	
	// read in the data
	// assumes each column has: x,y,z, dV
	
	bounding_box[0] = 9e99; 
	bounding_box[1] = 9e99; 
	bounding_box[2] = 9e99; 

	bounding_box[3] = -9e99; 
	bounding_box[4] = -9e99; 
	bounding_box[5] = -9e99; 
 
        size_t result;
	for( unsigned int i=0; i < number_of_data_entries ; i++ )
	{
		result = fread( (char*) &( voxels[i].center[0] ) , sizeof(double) , 1 , fp ); 
		result = fread( (char*) &( voxels[i].center[1] ) , sizeof(double) , 1 , fp ); 
		result = fread( (char*) &( voxels[i].center[2] ) , sizeof(double) , 1 , fp ); 
		result = fread( (char*) &( voxels[i].volume ) , sizeof(double) , 1 , fp ); 
		
		// estimate the bounding box; 
		if( voxels[i].center[0] < bounding_box[0] )
		{ bounding_box[0] = voxels[i].center[0]; }
		if( voxels[i].center[0] > bounding_box[3] )
		{ bounding_box[3] = voxels[i].center[0]; }

		if( voxels[i].center[1] < bounding_box[1] )
		{ bounding_box[1] = voxels[i].center[1]; }
		if( voxels[i].center[1] > bounding_box[4] )
		{ bounding_box[4] = voxels[i].center[1]; }

		if( voxels[i].center[2] < bounding_box[2] )
		{ bounding_box[2] = voxels[i].center[2]; }
		if( voxels[i].center[2] > bounding_box[5] )
		{ bounding_box[5] = voxels[i].center[2]; }
	} 
	
	std::cout << "Warning: General_Mesh::read_from_matlab is incomplete. No connection information read." << std::endl; 

	fclose( fp) ; 
	return; 
}

/* Cartesian meshes */

Cartesian_Mesh::Cartesian_Mesh()
{
	Cartesian_mesh = true; 
	uniform_mesh = true; 
	regular_mesh = true; 
	use_voxel_faces = false; 
	
	x_coordinates.assign( 1 , 0.0 ); 
	y_coordinates.assign( 1 , 0.0 ); 
	z_coordinates.assign( 1 , 0.0 );
	
	/*------------------------------------------------------------------------------------*/ 
	/* By calculating dx, dy and dz like this, I think they are assuming that there is 		*/
	/* only 1 voxel present initially and it spans the complete domain. 									*/
	/*------------------------------------------------------------------------------------*/ 

	
	dx = bounding_box[3] - bounding_box[0]; 
	dy = bounding_box[4] - bounding_box[1]; 
	dz = bounding_box[5] - bounding_box[2]; 
	
	static double tolerance = 1e-12; 
	
	if( fabs( dx - dy ) > tolerance ||
		fabs( dy - dz ) > tolerance ||
		fabs( dx - dz ) > tolerance )
	{ uniform_mesh = false; }
	
	dV = dx*dy*dz; 
	dS = dx*dy; 
	dS_xy = dx*dy; 
	dS_yz = dy*dz; 
	dS_xz = dx*dz; 
	
	Voxel template_voxel;
	template_voxel.volume = dV; 
 
	voxels.assign( x_coordinates.size() * y_coordinates.size() * z_coordinates.size() , template_voxel ); 
	voxels[0].center[0] = x_coordinates[0]; 
	voxels[0].center[1] = y_coordinates[0]; 
	voxels[0].center[2] = z_coordinates[0]; 
}


void Cartesian_Mesh::create_voxel_faces( void )
{
	// make connections 
	
	connected_voxel_indices.resize( voxels.size() ); 
	
	int i_jump = y_coordinates.size() * z_coordinates.size(); 
	int j_jump = y_coordinates.size(); 
	int k_jump = 1; 
		
	// x-aligned connections 
	for( unsigned int k=0 ; k < z_coordinates.size() ; k++ )
	{
		for( unsigned int j=0 ; j < y_coordinates.size() ; j++ )
		{
			for( unsigned int i=0 ; i < x_coordinates.size()-1 ; i++ )
			{
				int n = voxel_index(i,j,k); 
				connect_voxels_faces_only(n,n+i_jump, dS_yz ); 
			}
		}
	}
	// y-aligned connections 
	for( unsigned int k=0 ; k < z_coordinates.size() ; k++ )
	{
		for( unsigned int i=0 ; i < x_coordinates.size() ; i++ )
		{
			for( unsigned int j=0 ; j < y_coordinates.size()-1 ; j++ )
			{
				int n = voxel_index(i,j,k); 
				connect_voxels_faces_only(n,n+j_jump, dS_xz ); 
			}
		}
	}	
	// z-aligned connections 
	for( unsigned int j=0 ; j < y_coordinates.size() ; j++ )
	{
		for( unsigned int i=0 ; i < x_coordinates.size() ; i++ )
		{
			for( unsigned int k=0 ; k < z_coordinates.size()-1 ; k++ )
			{
				int n = voxel_index(i,j,k); 
				connect_voxels_faces_only(n,n+k_jump, dS_xy ); 
			}
		}
	}	

	return; 
}

Cartesian_Mesh::Cartesian_Mesh( int xnodes, int ynodes, int znodes )
{
	x_coordinates.assign( xnodes , 0.0 ); 
	y_coordinates.assign( ynodes , 0.0 ); 
	z_coordinates.assign( znodes , 0.0 ); 
	
	dx = 1;
	dy = 1; 
	dz = 1; 
	
	dV = dx*dx*dz; 
	dS = dx*dy; 
	
	dS_xy = dS; 
	dS_yz = dS; 
	dS_xz = dS; 

	uniform_mesh = true; 
	regular_mesh = true; 
	use_voxel_faces = false; 
 
	for( unsigned int i=0; i < x_coordinates.size() ; i++ )
	{ x_coordinates[i] = i*dx; }
	for( unsigned int i=0; i < y_coordinates.size() ; i++ )
	{ y_coordinates[i] = i*dy; }
	for( unsigned int i=0; i < z_coordinates.size() ; i++ )
	{ z_coordinates[i] = i*dz; }	
	
	bounding_box[0] = x_coordinates[0]-dx/2.0; 
	bounding_box[3] = x_coordinates[x_coordinates.size()-1]+dx/2.0; 
	bounding_box[1] = y_coordinates[0]-dy/2.0; 
	bounding_box[4] = y_coordinates[y_coordinates.size()-1]+dy/2.0; 
	bounding_box[2] = z_coordinates[0]-dz/2.0; 
	bounding_box[5] = z_coordinates[z_coordinates.size()-1]+dz/2.0; 
	
	Voxel template_voxel;
	template_voxel.volume = dV; 

	units = "none"; 
	
	voxels.assign( x_coordinates.size() * y_coordinates.size() * z_coordinates.size() , template_voxel ); 

	// initializing and connecting voxels 
 
	int n=0; 
	for( unsigned int i=0 ; i < x_coordinates.size() ; i++ )
	{
		for( unsigned int j=0 ; j < y_coordinates.size() ; j++ )
		{
			for( unsigned int k=0 ; k < z_coordinates.size() ; k++ )
			{
				voxels[n].center[0] = x_coordinates[i]; 
				voxels[n].center[1] = y_coordinates[j]; 
				voxels[n].center[2] = z_coordinates[k]; 
				voxels[n].mesh_index = n; 
				voxels[n].volume = dV; 

				n++; 
			}
		}
	}
	
	// make connections 
	
	connected_voxel_indices.resize( voxels.size() ); 
	
	int i_jump = y_coordinates.size()*z_coordinates.size(); 
	int j_jump = z_coordinates.size(); 
	int k_jump = 1; 
		
	// x-aligned connections 
	for( unsigned int k=0 ; k < z_coordinates.size() ; k++ )
	{
		for( unsigned int j=0 ; j < y_coordinates.size() ; j++ )
		{
			for( unsigned int i=0 ; i < x_coordinates.size()-1 ; i++ )
			{
				int n = voxel_index(i,j,k); 
				connect_voxels_indices_only(n,n+i_jump, dS_yz ); 
			}
		}
	}
	// y-aligned connections 
	for( unsigned int k=0 ; k < z_coordinates.size() ; k++ )
	{
		for( unsigned int i=0 ; i < x_coordinates.size() ; i++ )
		{
			for( unsigned int j=0 ; j < y_coordinates.size()-1 ; j++ )
			{
				int n = voxel_index(i,j,k); 
				connect_voxels_indices_only(n,n+j_jump, dS_xz ); 
			}
		}
	}	
	// z-aligned connections 
	for( unsigned int j=0 ; j < y_coordinates.size() ; j++ )
	{
		for( unsigned int i=0 ; i < x_coordinates.size() ; i++ )
		{
			for( unsigned int k=0 ; k < z_coordinates.size()-1 ; k++ )
			{
				int n = voxel_index(i,j,k); 
				connect_voxels_indices_only(n,n+k_jump, dS_xy ); 
			}
		}
	}	
	
	if( use_voxel_faces )
	{ create_voxel_faces(); }
}	 

void Cartesian_Mesh::create_moore_neighborhood()
{
	//BioFVM-B optimized generation of creating moore neighbours
	//Loops can be parallelized and collapsed but it has shown downgrades in performance
	moore_connected_voxel_indices.resize(voxels.size());
	for (int i = 0; i < x_coordinates.size(); i++)
	{
		for (int j = 0; j < y_coordinates.size(); j++)
		{
			for (int k = 0; k < z_coordinates.size(); k++)
			{
				int center_inex = voxel_index(i, j, k);

				int size = 0;
				for (int ii = -1; ii <= 1; ii++)
					for (int jj = -1; jj <= 1; jj++)
						for (int kk = -1; kk <= 1; kk++)
							if (i + ii >= 0 && i + ii < x_coordinates.size() &&
								j + jj >= 0 && j + jj < y_coordinates.size() &&
								k + kk >= 0 && k + kk < z_coordinates.size() &&
								!(ii == 0 && jj == 0 && kk == 0))
							{
								++size;
							}
				moore_connected_voxel_indices[center_inex].resize(size, 0);
				int neighbor = 0;
				for (int ii = -1; ii <= 1; ii++)
					for (int jj = -1; jj <= 1; jj++)
						for (int kk = -1; kk <= 1; kk++)
							if (i + ii >= 0 && i + ii < x_coordinates.size() &&
								j + jj >= 0 && j + jj < y_coordinates.size() &&
								k + kk >= 0 && k + kk < z_coordinates.size() &&
								!(ii == 0 && jj == 0 && kk == 0))
							{
								int neighbor_index = voxel_index(i + ii, j + jj, k + kk);
								moore_connected_voxel_indices[center_inex][neighbor] = neighbor_index;
								++neighbor;
							}
			}
		}
	}
}

/*-------------------------------------------------------------*/
/* Parallel version of create_moore_neighborhood 							 */
/*-------------------------------------------------------------*/

void Cartesian_Mesh::create_moore_neighborhood(mpi_Environment &world, mpi_Cartesian &cart_topo)
{
	/*--------------------------------------------------------------------------------------*/
	/* I need to divide moore neighborhood construction into 3 parts: 											*/
	/* (1) Non-subdomain boundary voxels: create a list of local 26 neighbouring voxels 		*/
	/* (2) Left sub-domain boundary voxels: create a list of global 9 voxels, residing 			*/
	/* in neighbouring process only, in moore_connected_global_voxel_indices_left  					*/
	/* (3) Right sub-domain boundary voxels: create a list of global 9 voxels only 					*/
	/* residing in neighbouring process only, in moore_connected_global_voxel_indices_right */
	/*--------------------------------------------------------------------------------------*/
	
	/* Create local list of voxels NOT ONLY for non-boundary voxels BUT ALSO for boundary voxels 			*/
	/* because we want the local neighbouring voxels for boundary voxels to be contained in this list */
	/* and we will make 2 more separate lists */
	
	/* For optimization loop-order can be changed from outer to inner : z then y then x */
	
	moore_connected_voxel_indices.resize( voxels.size() );
	
	for( unsigned int j=0 ; j < y_coordinates.size() ; j++ )
	{
		for( unsigned int i=0 ; i < x_coordinates.size() ; i++ )
		{
			for( unsigned int k=0 ; k < z_coordinates.size() ; k++ )
			{
				int center_inex = voxel_index(i,j,k); 
				for(int ii=-1;ii<=1;ii++)
					for(int jj=-1;jj<=1;jj++)
						for(int kk=-1;kk<=1;kk++)
							if(i+ii >= 0  && i+ii < x_coordinates.size() &&
								 j+jj >= 0  && j+jj < y_coordinates.size() &&
								 k+kk >= 0  && k+kk < z_coordinates.size() &&
								 !(ii==0  && jj==0 && kk==0)
								)
								 {
									int neighbor_index = voxel_index(i+ii,j+jj,k+kk);
									moore_connected_voxel_indices[center_inex].push_back( neighbor_index );
								 }
			}
		}
	}
	
	/* Allocate space to _left/_right moore lists = number of voxels in yz-plane */
	
	moore_connected_voxel_global_indices_left.resize (y_coordinates.size() * z_coordinates.size());
	moore_connected_voxel_global_indices_right.resize(y_coordinates.size() * z_coordinates.size());
	
	
	int x_nodes = (bounding_box[3]-bounding_box[0])/dx;
	int y_nodes = (bounding_box[4]-bounding_box[1])/dy;
	int z_nodes = (bounding_box[5]-bounding_box[2])/dz;
	
	int global_z_jump = 1;
	int global_y_jump = z_nodes; 
	int global_x_jump = z_nodes * y_nodes;
	
	/* Now build the moore list of global voxel indices for the left sub-domain boundary voxels */
	
	if(world.rank > 0)
	{
		
		int vxl_indx_ctr = 0;
		
		for(int j=0; j<y_coordinates.size(); j++)
		{
			for(int k=0; k<z_coordinates.size(); k++)
			{
				int center_inex_local  = voxel_index(0, j, k); 
				int center_inex_global = voxels[center_inex_local].global_mesh_index;
				
				/* First tackle Inner region (IC of a sub-domain !) 																*/
				/* Such a voxel will have 9 global neighbouring voxels on neighbouring left process */
				
				if(j > 0 && j < y_coordinates.size()-1 && k > 0 && k < z_coordinates.size()-1)
				{
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump); 								 // LEFT-1
				
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_y_jump); // LEFT-1,UP+1
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_y_jump); // LEFT-1,DOWN-1
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_z_jump); // LEFT-1,IN+1	
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_z_jump); // LEFT-1,OUT-1
				
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_y_jump + global_z_jump); //LEFT-1,UP+1,IN+1
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_y_jump - global_z_jump); //LEFT-1,UP+1,OUT-1
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_y_jump + global_z_jump); //LEFT-1,DOWN-1,IN+1
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_y_jump - global_z_jump); //LEFT-1,DOWN-1,OUT-1
				}
				
				/* Now tackle 4 extreme corner voxels : SW, SE, NW, NE */
				/* Such voxels will have 4 neighbouring voxels in neighbouring left process */
				
				if(j == 0 && k == 0)
				{
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_y_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_z_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_y_jump + global_z_jump);
				}
				
				if(j == y_coordinates.size()-1 && k == 0)
				{
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump); // LEFT-1
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_y_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_z_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_y_jump + global_z_jump);
				}
				
				if(j == 0 && k == z_coordinates.size()-1)
				{
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_y_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_z_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_y_jump - global_z_jump);
				}
				
				if(j == y_coordinates.size()-1 && k == z_coordinates.size()-1)
				{
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_y_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_z_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_y_jump - global_z_jump);
				}
				
				/* Now tackle 2 vertical strips and 2 horizontal strips not containing the 4 extreme corner voxels  */
				/* Such voxels have 6 neighbouring voxels in neighbouring left process 															*/
				
				/* Front vertical strip */
				if(k == 0 && j > 0 && j < y_coordinates.size()-1)
				{
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_y_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_y_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_z_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_z_jump + global_y_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_z_jump - global_y_jump);
				}
				
				/* Back Vertical strip */
				if(k == z_coordinates.size()-1 && j > 0 && j < y_coordinates.size()-1)
				{
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_y_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_y_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_z_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_z_jump + global_y_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_z_jump - global_y_jump);
				}
				
				/* Bottom horizontal strip */
				
				if(j == 0 && k > 0 && k < z_coordinates.size()-1)
				{
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_z_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_z_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_y_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_y_jump + global_z_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_y_jump - global_z_jump);
				}
				
				/* Top horizontal strip */
				
				if(j == y_coordinates.size()-1 && k > 0 && k < z_coordinates.size()-1)
				{
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump + global_z_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_z_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_y_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_y_jump + global_z_jump);
					moore_connected_voxel_global_indices_left[vxl_indx_ctr].push_back(center_inex_global-global_x_jump - global_y_jump - global_z_jump);						
				}
				
				vxl_indx_ctr = vxl_indx_ctr + 1;
			}
		}
	}
	
	/* Now build the moore list of global voxel indices for the right sub-domain boundary voxels */

	if(world.rank < (world.size-1))
	{
		int vxl_indx_ctr = 0;
		
		for(int j=0; j<y_coordinates.size(); j++)
		{
			for(int k=0; k<z_coordinates.size(); k++)
			{
				int center_inex_local  = voxel_index(x_coordinates.size()-1, j, k); 
				int center_inex_global = voxels[center_inex_local].global_mesh_index;
				
				/* First tackle Inner region (IC of a sub-domain !) 																*/
				/* Such a voxel will have 9 global neighbouring voxels on neighbouring right process */
				
				if(j > 0 && j < y_coordinates.size()-1 && k > 0 && k < z_coordinates.size()-1)
				{
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump); 								 // RIGHT+1
				
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_y_jump); // RIGHT+1,UP+1
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_y_jump); // RIGHT+1,DOWN-1
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_z_jump); // RIGHT+1,IN+1	
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_z_jump); // RIGHT+1,OUT-1
				
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_y_jump + global_z_jump); //RIGHT+1,UP+1,IN+1
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_y_jump - global_z_jump); //RIGHT+1,UP+1,OUT-1
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_y_jump + global_z_jump); //RIGHT+1,DOWN-1,IN+1
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_y_jump - global_z_jump); //RIGHT+1,DOWN-1,OUT-1
				}
				
				/* Now tackle 4 extreme corner voxels : SW, SE, NW, NE */
				/* Such voxels will have 4 neighbouring voxels in neighbouring right process */
				
				if(j == 0 && k == 0)
				{
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_y_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_z_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_y_jump + global_z_jump);
				}
				
				if(j == y_coordinates.size()-1 && k == 0)
				{
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump); 
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_y_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_z_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_y_jump + global_z_jump);
				}
				
				if(j == 0 && k == z_coordinates.size()-1)
				{
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_y_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_z_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_y_jump - global_z_jump);
				}
				
				if(j == y_coordinates.size()-1 && k == z_coordinates.size()-1)
				{
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_y_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_z_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_y_jump - global_z_jump);
				}
				
				/* Now tackle 2 vertical strips and 2 horizontal strips not containing the 4 extreme corner voxels  */
				/* Such voxels have 6 neighbouring voxels in neighbouring right process 															*/
				
				/* Front vertical strip */
				
				if(k == 0 && j > 0 && j < y_coordinates.size()-1)
				{
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_y_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_y_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_z_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_z_jump + global_y_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_z_jump - global_y_jump);
				}
				
				/* Back Vertical strip */
				
				if(k == z_coordinates.size()-1 && j > 0 && j < y_coordinates.size()-1)
				{
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_y_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_y_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_z_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_z_jump + global_y_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_z_jump - global_y_jump);
				}
				
				/* Bottom horizontal strip */
				
				if(j == 0 && k > 0 && k < z_coordinates.size()-1)
				{
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_z_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_z_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_y_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_y_jump + global_z_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_y_jump - global_z_jump);
				}
				
				/* Top horizontal strip */
				
				if(j == y_coordinates.size()-1 && k > 0 && k < z_coordinates.size()-1)
				{
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump + global_z_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_z_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_y_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_y_jump + global_z_jump);
					moore_connected_voxel_global_indices_right[vxl_indx_ctr].push_back(center_inex_global+global_x_jump - global_y_jump - global_z_jump);						
				}
				
				vxl_indx_ctr = vxl_indx_ctr + 1;
			}
		}	
	}
}

unsigned int Cartesian_Mesh::voxel_index( unsigned int i, unsigned int j, unsigned int k )
{
    
 /*-------------------------------------------------------------*/   
 /* Returns the loval index of the voxel, given the coordinates */ 
 /* size of x/y/z_coordinates is same as local_x/y/z_nodes      */
 /*-------------------------------------------------------------*/
    
 return ( i * y_coordinates.size() * z_coordinates.size() + j * z_coordinates.size() + k) ; //New indexation---> Jose 
}

std::vector<unsigned int> Cartesian_Mesh::cartesian_indices( unsigned int n )
{
	std::vector<unsigned int> out(3, -1 ); 

	//Jose: new indexation
	// figure out i; 
	unsigned int YZ = y_coordinates.size() * z_coordinates.size();
	out[0] = (unsigned int) floor( n/YZ ); 
 
	// figure out j; 
	out[1] = (unsigned int) floor(   (n - out[0]*YZ) / z_coordinates.size() );
 
	// figure out k; 
	out[2] = n - z_coordinates.size()*(   out[1] + y_coordinates.size()*out[0] ); 

	return out; 
}

//Resize in single node execution
void Cartesian_Mesh::resize( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , int x_nodes, int y_nodes, int z_nodes )
{
	x_coordinates.assign( x_nodes , 0.0 ); 
	y_coordinates.assign( y_nodes , 0.0 ); 
	z_coordinates.assign( z_nodes , 0.0 ); 

	dx = ( x_end - x_start )/( (double) x_nodes ); 
	if( x_nodes < 2 )
	{ dx = 1; }
	dy = ( y_end - y_start )/( (double) y_nodes ); 
	if( y_nodes < 2 )
	{ dy = 1; }
	dz = ( z_end - z_start )/( (double) z_nodes  ); 
	if( z_nodes < 2 )
	{ dz = 1; }

	uniform_mesh = true; 
	regular_mesh = true; 
	static double tol = 1e-16; 
 
	if( fabs( dx - dy ) > tol && x_nodes > 1 && y_nodes > 1 )
	{ uniform_mesh = false; }
	if( fabs( dy - dz ) > tol && y_nodes > 1 && z_nodes > 1 )
	{ uniform_mesh = false; }
	if( fabs( dx - dz ) > tol && x_nodes > 1 && z_nodes > 1 )
	{ uniform_mesh = false; }

	for( unsigned int i=0; i < x_coordinates.size() ; i++ )
	{ x_coordinates[i] = x_start + (i+0.5)*dx; }
	for( unsigned int i=0; i < y_coordinates.size() ; i++ )
	{ y_coordinates[i] = y_start + (i+0.5)*dy; }
	for( unsigned int i=0; i < z_coordinates.size() ; i++ )
	{ z_coordinates[i] = z_start + (i+0.5)*dz; }

	bounding_box[0] = x_start; 
	bounding_box[3] = x_end; 
	bounding_box[1] = y_start; 
	bounding_box[4] = y_end; 
	bounding_box[2] = z_start; 
	bounding_box[5] = z_end; 

	dV = dx*dy*dz; 
	dS = dx*dy; 
	
	dS_xy = dx*dy; 
	dS_yz = dy*dz; 
	dS_xz = dx*dz; 

	Voxel template_voxel;
	template_voxel.volume = dV; 

	voxels.assign( x_coordinates.size() * y_coordinates.size() * z_coordinates.size() , template_voxel ); 

	int n=0;
	for( unsigned int i=0 ; i < x_coordinates.size() ; i++ ) 
	{
		for( unsigned int j=0 ; j < y_coordinates.size() ; j++ )
		{
			for( unsigned int k=0 ; k < z_coordinates.size() ; k++ )
			{
				voxels[n].center[0] = x_coordinates[i]; 
				voxels[n].center[1] = y_coordinates[j]; 
				voxels[n].center[2] = z_coordinates[k]; 
				voxels[n].mesh_index = n; 
				voxels[n].volume = dV; 

				n++; 
			}
		}
	}
	
	// make connections 
	
	connected_voxel_indices.resize( voxels.size() ); 
	voxel_faces.clear(); 
	
	for( unsigned int i=0; i < connected_voxel_indices.size() ; i++ )
	{ connected_voxel_indices[i].clear(); }
	int i_jump = y_coordinates.size() * z_coordinates.size(); 
	int j_jump = z_coordinates.size(); 
	int k_jump = 1; 
	
	// x-aligned connections 
	int count = 0; 
	for( unsigned int k=0 ; k < z_coordinates.size() ; k++ )
	{
		for( unsigned int j=0 ; j < y_coordinates.size() ; j++ )
		{
			for( unsigned int i=0 ; i < x_coordinates.size()-1 ; i++ )
			{
				int n = voxel_index(i,j,k); 
				connect_voxels_indices_only(n,n+i_jump, dS_yz ); 
				count++; 
			}
		}
	}
	// y-aligned connections 
	for( unsigned int k=0 ; k < z_coordinates.size() ; k++ )
	{
		for( unsigned int i=0 ; i < x_coordinates.size() ; i++ )
		{
			for( unsigned int j=0 ; j < y_coordinates.size()-1 ; j++ )
			{
				int n = voxel_index(i,j,k); 
				connect_voxels_indices_only(n,n+j_jump, dS_xz ); 
			}
		}
	}	
	// z-aligned connections 
	for( unsigned int j=0 ; j < y_coordinates.size() ; j++ )
	{
		for( unsigned int i=0 ; i < x_coordinates.size() ; i++ )
		{
			for( unsigned int k=0 ; k < z_coordinates.size()-1 ; k++ )
			{
				int n = voxel_index(i,j,k); 
				connect_voxels_indices_only(n,n+k_jump, dS_xy ); 
			}
		}
	}
	
	if( use_voxel_faces )
	{ create_voxel_faces(); }
	
	create_moore_neighborhood();
	return; 
}

void Cartesian_Mesh::resize( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx_new, double dy_new , double dz_new )
{
	dx = dx_new;
	dy = dy_new; 
	dz = dz_new; 

	double eps = 1e-16; 
	int x_nodes = (int) ceil( eps + (x_end-x_start)/dx ); 
	int y_nodes = (int) ceil( eps + (y_end-y_start)/dy ); 
	int z_nodes = (int) ceil( eps + (z_end-z_start)/dz ); 

	x_coordinates.assign( x_nodes , 0.0 ); 
	y_coordinates.assign( y_nodes , 0.0 ); 
	z_coordinates.assign( z_nodes , 0.0 ); 

	uniform_mesh = true; 
	regular_mesh = true; 
	double tol = 1e-16; 
	if( fabs( dx - dy ) > tol || fabs( dy - dz ) > tol || fabs( dx - dz ) > tol )
	{ uniform_mesh = false; }

	for( unsigned int i=0; i < x_coordinates.size() ; i++ )
	{ x_coordinates[i] = x_start + (i+0.5)*dx; }
	for( unsigned int i=0; i < y_coordinates.size() ; i++ )
	{ y_coordinates[i] = y_start + (i+0.5)*dy; }
	for( unsigned int i=0; i < z_coordinates.size() ; i++ )
	{ z_coordinates[i] = z_start + (i+0.5)*dz; }

	bounding_box[0] = x_start; 
	bounding_box[3] = x_end; 
	bounding_box[1] = y_start; 
	bounding_box[4] = y_end; 
	bounding_box[2] = z_start; 
	bounding_box[5] = z_end; 

	dV = dx*dy*dz; 
	dS = dx*dy; 

	dS_xy = dx*dy; 
	dS_yz = dy*dz; 
	dS_xz = dx*dz; 
	
	Voxel template_voxel;
	template_voxel.volume = dV; 

	voxels.assign( x_coordinates.size() * y_coordinates.size() * z_coordinates.size() , template_voxel ); 
	
	int n=0; 
	for( unsigned int i=0 ; i < x_coordinates.size() ; i++ )
	{
		for( unsigned int j=0 ; j < y_coordinates.size() ; j++ )
		{
			for( unsigned int k=0 ; k < z_coordinates.size() ; k++ )
			{
				voxels[n].center[0] = x_coordinates[i]; 
				voxels[n].center[1] = y_coordinates[j]; 
				voxels[n].center[2] = z_coordinates[k]; 
				voxels[n].mesh_index = n; 
				voxels[n].volume = dV; 

				n++; 
			}
		}
	}
	
	// make connections 
	
	connected_voxel_indices.resize( voxels.size() ); 
	voxel_faces.clear(); 
	
	for( unsigned int i=0; i < connected_voxel_indices.size() ; i++ )
	{ connected_voxel_indices[i].clear(); }
	
	int i_jump = y_coordinates.size() * z_coordinates.size(); 
	int j_jump = z_coordinates.size(); 
	int k_jump = 1; 
	
	// x-aligned connections 
	int count = 0; 
	for( unsigned int k=0 ; k < z_coordinates.size() ; k++ )
	{
		for( unsigned int j=0 ; j < y_coordinates.size() ; j++ )
		{
			for( unsigned int i=0 ; i < x_coordinates.size()-1 ; i++ )
			{
				int n = voxel_index(i,j,k); 
				connect_voxels_indices_only(n,n+i_jump, dS_yz ); 
				count++; 
			}
		}
	}

	// y-aligned connections 
	for( unsigned int k=0 ; k < z_coordinates.size() ; k++ )
	{
		for( unsigned int i=0 ; i < x_coordinates.size() ; i++ )
		{
			for( unsigned int j=0 ; j < y_coordinates.size()-1 ; j++ )
			{
				int n = voxel_index(i,j,k); 
				connect_voxels_indices_only(n,n+j_jump, dS_xz ); 
			}
		}
	}	

	// z-aligned connections 
	for( unsigned int j=0 ; j < y_coordinates.size() ; j++ )
	{
		for( unsigned int i=0 ; i < x_coordinates.size() ; i++ )
		{
			for( unsigned int k=0 ; k < z_coordinates.size()-1 ; k++ )
			{
				int n = voxel_index(i,j,k); 
				connect_voxels_indices_only(n,n+k_jump, dS_xy ); 
			}
		}
	}
	
	if( use_voxel_faces )
	{ create_voxel_faces(); }
	
	create_moore_neighborhood();
	return; 
}

/*=============================================*/
/* Parallel version of resize() function above */
/*=============================================*/

void Cartesian_Mesh::resize( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx_new, double dy_new , double dz_new, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
    int dims[3], coords[3]; 
    
	/*------------------------------------------------------------------------------------------------*/
	/*		local_x/y/z_start give local starting coordinates at each process	                      */
	/*------------------------------------------------------------------------------------------------*/
    
    double local_x_start;
    double local_y_start;
    double local_z_start;
    
  /*------------------------------------------------------------------------------------------------*/
	/*		To find global_mesh_index, we declare following	                                          */
	/*------------------------------------------------------------------------------------------------*/
    
    int x_index;
    int y_index;
    int z_index;
    int local_start_of_global_index; 
    
    /*-------------------------------------------*/
    /* Just to reduce typing _new everytime .... */
    /*-------------------------------------------*/
    
    dx = dx_new;
	dy = dy_new; 
	dz = dz_new; 

	double eps = 1e-16; 
    
  /*--------------------------------------------*/
	/*		Global Nodes	                      */
	/*--------------------------------------------*/
    
	int x_nodes = (int) ceil( eps + (x_end-x_start)/dx ); 
	int y_nodes = (int) ceil( eps + (y_end-y_start)/dy ); 
	int z_nodes = (int) ceil( eps + (z_end-z_start)/dz ); 
    
    /*-----------------------------------------------------------------------------------------------------------------*/
    /* Need to store dims[index] = cart_topo.dims[index] so we can use dims array easily, need to declare this array above */
    /*-----------------------------------------------------------------------------------------------------------------*/
    
    dims[0] = cart_topo.mpi_dims[0]; 
    dims[1] = cart_topo.mpi_dims[1];
    dims[2] = cart_topo.mpi_dims[2];
    
  /*--------------------------------------------*/
	/*		Local Nodes on MPI Processes	      */
	/*--------------------------------------------*/
    
    int local_x_nodes = x_nodes/dims[1];
    
   /*----------------------------------------------------------------*/ 
   /* Here I need to insert a perfect divisibility check and message */
   /*----------------------------------------------------------------*/ 
   
   //TODO: residuals to be allocated in first ranks
   // if (x_nodes % dims[1] > mpi_rank) ++local_x_nodes;
   if(x_nodes % dims[1] != 0)
   {
   	if(world.rank == 0)
   		{
   			std::cout<<"Error: Total voxels in X-Direction are NOT perfectly divisible by the total MPI processes"<<std::endl;
   			std::cout<<"Aborting the program"<<std::endl;
   		}
   	MPI_Abort(cart_topo.mpi_cart_comm, -1); 	
   }
    
    int local_y_nodes = y_nodes/dims[0];
    int local_z_nodes = z_nodes/dims[2];

	/*------------------------------------------------------------------------------------------*/
	/*		Assign space to coordinate arrays according to 'local' nodes and NOT 'global' nodes */
	/*------------------------------------------------------------------------------------------*/
    
    x_coordinates.assign( local_x_nodes , 0.0 ); 
	y_coordinates.assign( local_y_nodes , 0.0 ); 
	z_coordinates.assign( local_z_nodes , 0.0 ); 

	uniform_mesh = true; 
	regular_mesh = true; 
	double tol = 1e-16; 
    
	if( fabs( dx - dy ) > tol || fabs( dy - dz ) > tol || fabs( dx - dz ) > tol )
	{ 
        uniform_mesh = false; 
    }
    
    /*------------------------------------------------------------------------------------------*/
    /* Need to store coordinates of processes from car_topo to the coords array for convenience */
    /*------------------------------------------------------------------------------------------*/
    
    coords[0] = cart_topo.mpi_coords[0]; 
    coords[1] = cart_topo.mpi_coords[1]; 
    coords[2] = cart_topo.mpi_coords[2]; 
    
    /*-----------------------------------------------------------------------------------------------------------------*/
    /* First find local x coordinate start of each process = global_x_start + process_y_coodinate * x_local_nodes * dx */
    /*-----------------------------------------------------------------------------------------------------------------*/
    
    local_x_start = x_start + (coords[1] * local_x_nodes * dx);
    
	for( unsigned int i=0; i < x_coordinates.size() ; i++ )
	{ 
        x_coordinates[i] = local_x_start + (i+0.5)*dx;         
    }
    
    /*----------------------------------------------------------------------------------------------------------------------*/
    /* First find local y coordinate start of each process = global_y_start + process_x_coodinate * y_local_nodes * dy      */
    /* BioFVM y-coordinate is increasing from bottom to top while MPI process x-coordinate is decreasing from bottom to top */
    /*----------------------------------------------------------------------------------------------------------------------*/
    
    local_y_start = y_start + ((dims[0] - coords[0] - 1) * local_y_nodes * dy);
    
	for( unsigned int i=0; i < y_coordinates.size() ; i++ )
	{ 
        y_coordinates[i] = local_y_start + (i+0.5)*dy; 
    }
    
	/*-----------------------------------------------------------------------------------------------------------------*/
  /* First find local z coordinate start of each process = global_z_start + process_z_coodinate * z_local_nodes * dz */
  /*-----------------------------------------------------------------------------------------------------------------*/
	
    local_z_start = z_start + (coords[2] * local_z_nodes * dz);
    
	for( unsigned int i=0; i < z_coordinates.size() ; i++ )
	{ 
        z_coordinates[i] = local_z_start + (i+0.5)*dz;         
    }

  /*--------------------------*/
	/* Global Domain dimensions */
	/*--------------------------*/
    
	bounding_box[0] = x_start; 
	bounding_box[3] = x_end; 
	
	bounding_box[1] = y_start; 
	bounding_box[4] = y_end; 
	
	bounding_box[2] = z_start; 
	bounding_box[5] = z_end;
	
	/*--------------------------------------*/
	/* Local sub-domain bounding box limits */
	/*--------------------------------------*/
	
	local_bounding_box[0] = local_x_start;
	local_bounding_box[3] = local_x_start + local_x_nodes * dx;
	
	local_bounding_box[1] = local_y_start;
	local_bounding_box[4] = local_y_start + local_y_nodes * dy;
	
	local_bounding_box[2] = local_z_start; 	
	local_bounding_box[5] = local_z_start + local_z_nodes * dz;

  /*----------------------------------------*/
	/* Voxel volume and various surface areas */
  /*----------------------------------------*/
    
  dV = dx*dy*dz; 
	dS = dx*dy; 

	dS_xy = dx*dy; 
	dS_yz = dy*dz; 
	dS_xz = dx*dz; 
	
  /*-------------------------------------------------------------------------------------------------------------*/
	/* Creates a Voxel defined in BioFVM_mesh.cpp to a default Voxel having index 0, center (0,0,0), volume = 1000 */
  /* The template Voxel now has a new volume i.e. template_voxel.volume = dV; see below                          */
  /* Added a new field global_mesh_index to Voxel class in BioFVM_Mesh.h as we need both local voxel index       */
  /* and global voxel index. 
	/*-------------------------------------------------------------------------------------------------------------*/
    
	Voxel template_voxel;
	template_voxel.volume = dV; 

    /*------------------------------------------------*/
    /* voxels is std::vector<Voxel> voxels;           */
    /* Size of x/y/z_coordinates is local_x/y/z_nodes */
    /* Next statement assigns voxels on local process */
    /*------------------------------------------------*/
    
    voxels.assign( x_coordinates.size() * y_coordinates.size() * z_coordinates.size() , template_voxel );
    
    
    /*--------------------------------------------------------------*/
    /* Next Loop assigns x,y, and z coordinates to voxel centers    */
    /* 2 indexes for every voxel are maintained: local & global     */
    /* The starting global index of the first voxel on each process */
    /* must be calculated to maintain global_mesh_index             */
    /*--------------------------------------------------------------*/
    //New layout BioFVM-B	
	//TODO: local_x_nodes could be not equal through all sub-domains
    //local_start_of_global_index = (coords[1] * z_nodes * y_nodes * local_x_nodes) +       //Imagine 3rd plate 'beginning' point (leftmost bottom point)
    //                              (dims[0]-coords[0]-1) * z_nodes * local_y_nodes +       //Imagine going up in 3rd plate
    //                              (coords[2] * local_z_nodes) ;                           //Imagine going right in 3rd plate
    local_start_of_global_index =  world.rank * z_nodes * y_nodes * local_x_nodes;                                     
	int n = 0;
	//#pragma omp parallel for collapse(3)
	for (int i = 0; i < x_coordinates.size(); i++)
	{
		for (int j = 0; j < y_coordinates.size(); j++)
		{
			for (int k = 0; k < z_coordinates.size(); k++)
			{
				int z_index = k;
				int y_index = j * z_nodes;  
				int x_index = i * y_nodes * z_nodes;
				voxels[n].center[0] = x_coordinates[i]; 
				voxels[n].center[1] = y_coordinates[j]; 
				voxels[n].center[2] = z_coordinates[k]; 
				voxels[n].mesh_index = n;                                                             //This now becomes the local index
				voxels[n].global_mesh_index = local_start_of_global_index + z_index + y_index + x_index;    //This is now the global index of the Voxel in the global mesh.
				voxels[n].volume = dV; 
				n++; 
			}
		}
	}
	
	/*--------------------------*/
	/* Make Connections next ...*/
  /*--------------------------*/
	
	connected_voxel_indices.resize( voxels.size() );             //This uses local indices
    
    /*-------------------------------------------------------------------------------------*/
    /* Another array that gives connection indices according to global numbering of voxels */
    /*-------------------------------------------------------------------------------------*/
    
    connected_voxel_global_indices.resize(voxels.size());         //--->This field in class Voxels has been added by Gaurav Saxena
    
	voxel_faces.clear(); 
	
	for( unsigned int i=0; i < connected_voxel_indices.size() ; i++ )
	{ 
        connected_voxel_indices[i].clear(); 
        connected_voxel_global_indices[i].clear();                //I forgot this while parallelizing BioFVM code, possibly won't matter but for completeness...
  }
	
	/*----------------------------------------------------------------------------------*
	/* Need to define two kinds of jumps between immediate neighbours in same direction */
  /* (1) Local jump and (2) Global jump.
  /*----------------------------------------------------------------------------------*/
	
	int i_jump = local_y_nodes * local_z_nodes; //Local x-jump
	int j_jump = local_z_nodes; //Local y-jump
	int k_jump = 1; //Local z-jump

	int i_global_jump = z_nodes*y_nodes; //Global x-jump
	int j_global_jump = z_nodes; //Global y-jump
	int k_global_jump = k_jump; //Global z-jump
	
	
	/*-------------------------------------------------------------------------------*/
	/* In a parallel environment, making connections will be divided into two parts: */ 
  /* (1) Inner part - left, right, front, back, top, bottom neighbours all exist   */
  /* (2) Boundary parts: immediate neighbours my exist across processes or may not */
  /* exist, if the process touches the boundary                                    */
  /*-------------------------------------------------------------------------------*/
	
	int x_size = x_coordinates.size();
	int y_size = y_coordinates.size();
	int z_size = z_coordinates.size();

	//Optimized generation of if connected voxels with 1 dimensional displacement
	for (int i = 0; i < x_size; ++i) {
		for (int j = 0; j < y_size; ++j){
			for (int k = 0; k < z_size; ++k) {
				int n = voxel_index(i, j, k);
				//X-neigbours
				if ( i > 0 and i < x_size - 1){
					connected_voxel_indices[n].push_back(n + i_jump);
					connected_voxel_indices[n].push_back(n - i_jump);
					connected_voxel_global_indices[n].push_back(voxels[n + i_jump].global_mesh_index);
					connected_voxel_global_indices[n].push_back(voxels[n - i_jump].global_mesh_index);
				}
				//Left boundary of each process
				if (i == 0) {
					if (voxels[n].center[0] - dx / 2 > x_start){
						// First connect this to right neighbour then right neighbour to this.
						connected_voxel_indices[n].push_back(n + i_jump);
						connected_voxel_global_indices[n].push_back(voxels[n + i_jump].global_mesh_index);
						connected_voxel_global_indices[n].push_back(voxels[n].global_mesh_index - i_global_jump);}
					else {
						connected_voxel_indices[n].push_back(n + i_jump);
						connected_voxel_global_indices[n].push_back(voxels[n + i_jump].global_mesh_index);}
				}
				//Right boundary of each process
				if (i == x_size -1){
					if (voxels[n].center[0] + dx / 2 < x_end){ // i.e. it is not a process aligned with right physical boundary
						connected_voxel_global_indices[n].push_back(voxels[n].global_mesh_index + i_global_jump);} // But there is a neighbour on next process, so use global index
				}
				//Y-aligned parts
				if ( j > 0 and j < y_size - 1){
					connected_voxel_indices[n].push_back(n + j_jump);
					connected_voxel_indices[n].push_back(n - j_jump);
					connected_voxel_global_indices[n].push_back(voxels[n + j_jump].global_mesh_index);
					connected_voxel_global_indices[n].push_back(voxels[n - j_jump].global_mesh_index);
				}
				if (j == 0 ){
					if (voxels[n].center[1] - dy / 2 > y_start) // i.e. it is not a process aligned with bottom physical boundary
					{
						connected_voxel_indices[n].push_back(n + j_jump);
						connected_voxel_global_indices[n].push_back(voxels[n + j_jump].global_mesh_index);
						connected_voxel_global_indices[n].push_back(voxels[n].global_mesh_index - j_global_jump); // But there is a neighbour on previous process, so use global index
					}
					else // It is the process that is aligned with left physical boundary
					{
						connected_voxel_indices[n].push_back(n + j_jump);
						connected_voxel_global_indices[n].push_back(voxels[n + j_jump].global_mesh_index);
					}
				}
				if (j == y_size -1) {
						if (voxels[n].center[1] + dy / 2 < y_end)            // i.e. it is not a process aligned with right physical boundary
					{
						connected_voxel_global_indices[n].push_back(voxels[n].global_mesh_index + j_global_jump); // But there is a neighbour on next process, so use global index
					}
				}
				//Z-aligned
				if ( k > 0 and k < z_size - 1){
					connected_voxel_indices[n].push_back(n + k_jump);
					connected_voxel_indices[n].push_back(n - k_jump);
					connected_voxel_global_indices[n].push_back(voxels[n + k_jump].global_mesh_index);
					connected_voxel_global_indices[n].push_back(voxels[n - k_jump].global_mesh_index);
				}
				if (k == 0) {
					if (voxels[n].center[2] - dz / 2 > z_start) // i.e. it is not a process aligned with bottom physical boundary
					{
						connected_voxel_indices[n].push_back(n + k_jump);
						connected_voxel_global_indices[n].push_back(voxels[n + k_jump].global_mesh_index);
						connected_voxel_global_indices[n].push_back(voxels[n].global_mesh_index - k_global_jump); // But there is a neighbour on previous process, so use global index
					}
					else // It is the process that is aligned with left physical boundary
					{
						connected_voxel_indices[n].push_back(n + k_jump);
						connected_voxel_global_indices[n].push_back(voxels[n + k_jump].global_mesh_index);
					}
				}
				if (k == z_size -1){
					if (voxels[n].center[2] + dz / 2 < z_end)            // i.e. it is not a process aligned with right physical boundary
					{
						connected_voxel_global_indices[n].push_back(voxels[n].global_mesh_index + k_global_jump); // But there is a neighbour on next process, so use global index
					}
				}
			}
		}
	}
	
	/*----------------------------------------------------------------------------------------*/
	/* Checking through DDT step-by-step execution, use_voxel_faces is false for this example */
    /* Leave parallelization for later                                                        */
    /*----------------------------------------------------------------------------------------*/
	
	if( use_voxel_faces )
	{ 
        create_voxel_faces(); 
    }
	
	/*------------------------------------------------------------------------------*/
	/* Creating moore neighbourhood is now parallelized in the sense that: 					*/
	/* Inner sub-domain region maintains 26 local voxel indices and sub-domain			*/
	/* boundary voxels maintain (26-9=17) local voxel indices												*/
	/* The remaining 9 neighbours (not within the process) are maintained as global	*/
	/* voxels indices list in 2 different lists i.e. _left/_right										*/
	/* Now call the parallelized version 																						*/
  /*------------------------------------------------------------------------------*/
    
	create_moore_neighborhood(world, cart_topo);
    
	return; 
}

void Cartesian_Mesh::resize( int x_nodes, int y_nodes, int z_nodes )
{ return resize( 0-.5, x_nodes-1+.5 , 0-.5 , y_nodes-1+.5 , 0-.5 , z_nodes - 1+.5 , x_nodes, y_nodes, z_nodes ); } 

void Cartesian_Mesh::resize_uniform( double x_start, double x_end, double y_start, double y_end, double z_start, double z_end , double dx_new )
{ return resize( x_start, x_end, y_start, y_end, z_start, z_end , dx_new, dx_new , dx_new ); }

int Cartesian_Mesh::nearest_voxel_index( std::vector<double>& position )
{
	unsigned int i = (unsigned int) floor( (position[0]-local_bounding_box[0])/dx ); 
	unsigned int j = (unsigned int) floor( (position[1]-bounding_box[1])/dy ); 
	unsigned int k = (unsigned int) floor( (position[2]-bounding_box[2])/dz ); 

	//  add some bounds checking -- truncate to inside the computational domain   

	if( i >= x_coordinates.size() ){ i = x_coordinates.size()-1; }
	if( i < 0 ){ i = 0; }

	if( j >= y_coordinates.size() ){ j = y_coordinates.size()-1; }
	if( j < 0 ){ j = 0; }

	if( k >= z_coordinates.size() ){ k = z_coordinates.size()-1; }
	if( k < 0 ){ k = 0; }

	return ( i*y_coordinates.size() + j )*z_coordinates.size() + k; 
}

int Cartesian_Mesh::nearest_lcl_voxel_index( std::vector<double>& position )
{
	unsigned int i = (unsigned int) floor( (position[0]-local_bounding_box[0])/dx ); 
	unsigned int j = (unsigned int) floor( (position[1]-local_bounding_box[1])/dy ); 
	unsigned int k = (unsigned int) floor( (position[2]-local_bounding_box[2])/dz ); 

	//  add some bounds checking -- truncate to inside the computational domain   

	if( i >= x_coordinates.size() ){ i = x_coordinates.size()-1; }
	if( i < 0 ){ i = 0; }

	if( j >= y_coordinates.size() ){ j = y_coordinates.size()-1; }
	if( j < 0 ){ j = 0; }

	if( k >= z_coordinates.size() ){ k = z_coordinates.size()-1; }
	if( k < 0 ){ k = 0; }

	return ( i*y_coordinates.size() + j )*z_coordinates.size() + k; 
}

/*-----------------------------------------------------------------------------------*/
/* Parallel new function that returns a process specific local index 								 */
/* in which the Basic_Agent resides. 																								 */
/*-----------------------------------------------------------------------------------*/
//Jose: voy por aqui
int Cartesian_Mesh::nearest_voxel_local_index( std::vector<double>& position, mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	/*----------------------------------------------------*/
    /* Routine should return the local index of the voxel */
    /* of the process having rank world.rank that contains*/
    /* the Basic_Agent. The local index is needed because */
    /* voxels[global_index] is not a valid position       */
    /* voxels[local_index] = some_global_index is ok      */  
    /*----------------------------------------------------*/
    
    /*----------------------------------------------------*/
    /* Coordinates of Voxel containing Basic_Agent, 			*/
    /* dx, dy and dz are members of Cartesian_Mesh class  */
    /*----------------------------------------------------*/
    
    int x_vox = (int) floor( (position[0]-bounding_box[0])/dx ); 
	int y_vox = (int) floor( (position[1]-bounding_box[1])/dy ); 
	int z_vox = (int) floor( (position[2]-bounding_box[2])/dz );
    
    /*----------------------------------------------------*/
    /* Global Voxels in each directions                   */
    /*----------------------------------------------------*/
    
    int global_num_x_voxels = (bounding_box[3]-bounding_box[0])/dx; 
    int global_num_y_voxels = (bounding_box[4]-bounding_box[1])/dy;
    int global_num_z_voxels = (bounding_box[5]-bounding_box[2])/dz;
    
    /*----------------------------------------------------*/
    /* Local Voxels in each directions                    */
    /*----------------------------------------------------*/
    
	//TODO. adjust local_x_voxels with residuals from mpi_dim division
    int local_num_x_voxels = (bounding_box[3]-bounding_box[0])/(cart_topo.mpi_dims[1] * dx); //Not true when not perfectly divisible
    int local_num_y_voxels = (bounding_box[4]-bounding_box[1])/(cart_topo.mpi_dims[0] * dy);
    int local_num_z_voxels = (bounding_box[5]-bounding_box[2])/(cart_topo.mpi_dims[2] * dz);
    
    /*---------------------------------------------------------------*/
    /* bounds checking - truncate to inside the computational domain */                   
    /*---------------------------------------------------------------*/

	if( x_vox >= global_num_x_voxels ){x_vox = global_num_x_voxels-1;}
    if( x_vox < 0 ){x_vox = 0;}

	if( y_vox >= global_num_y_voxels ){y_vox = global_num_y_voxels-1;}
    if( y_vox < 0 ){y_vox = 0;}

	if( z_vox >= global_num_z_voxels ){z_vox = global_num_z_voxels-1;}
    if( z_vox < 0 ){ z_vox = 0;}
    
    /*---------------------------------------------------------------*/
    /* Find process coordinates using mpi_Rank and mpi_Dims					 */
    /* Now there is no need as coordinates are in cart_topo object   */                   
    /*---------------------------------------------------------------*/
    
// 		 int prod12 = mpi_Dims[1] * mpi_Dims[2]; 
//     int proc_x_coord = floor(mpi_Rank/prod12);
//     int proc_y_coord = floor((mpi_Rank - proc_x_coord * prod12)/mpi_Dims[2]); 
//     int proc_z_coord = mpi_Rank - proc_x_coord * prod12 - proc_y_coord * mpi_Dims[2];
		
		int proc_x_coord = cart_topo.mpi_coords[0];
		int proc_y_coord = cart_topo.mpi_coords[1];
		int proc_z_coord = cart_topo.mpi_coords[2]; 
    
    /*---------------------------------------------------------------*/
    /* Calculate the X/Y/Z coordinate of the first voxel             */                   
    /* of the process (given its process coordinates as above)       */
    /* Remember X/Y Mesh direction and MPI are different             */
    /*---------------------------------------------------------------*/
    
    int proc_start_vox_x_coord = proc_y_coord * local_num_x_voxels;         
    int proc_start_vox_y_coord = (cart_topo.mpi_dims[0] - 1 - proc_x_coord) * local_num_y_voxels; 
    int proc_start_vox_z_coord = proc_z_coord * local_num_z_voxels;
    
    /*---------------------------------------------------------------*/
    /* Calculate the difference between x/y/z coord of the Voxel     */                   
    /* that contains the Basic_Agent and the first starting Voxel    */
    /* of that process. Clearly, this diff >= 0.                     */
    /*---------------------------------------------------------------*/
    
    int diff_x_coord = x_vox - proc_start_vox_x_coord; 
    int diff_y_coord = y_vox - proc_start_vox_y_coord;
    int diff_z_coord = z_vox - proc_start_vox_z_coord;
    
    /*---------------------------------------------------------------*/
    /* Now calculate how many voxels are between the starting voxel  */
    /* and the voxel that contains the Basic_Agent.                  */
    /*---------------------------------------------------------------*/
    
    //int process_local_index_of_voxel_containing_basic_agent = diff_z_coord * local_num_x_voxels * local_num_y_voxels + \
                                                              diff_y_coord * local_num_x_voxels + \
                                                              diff_x_coord;
	int process_local_index_of_voxel_containing_basic_agent = diff_x_coord * local_num_y_voxels * local_num_z_voxels + \
                                                              diff_y_coord * local_num_z_voxels + \
                                                              diff_z_coord;
                                                              
    // std::cout<<"position[0]="<<position[0]<<" Local Voxel Index="<< \
		// process_local_index_of_voxel_containing_basic_agent<<std::endl;  
    
    return (process_local_index_of_voxel_containing_basic_agent); 
}



std::vector<unsigned int> Cartesian_Mesh::nearest_cartesian_indices( std::vector<double>& position )
{
	std::vector<unsigned int> out; 
	out.assign(3, 0 ); 
	out[0] = (unsigned int) floor( (position[0]-bounding_box[0])/dx ); 
	out[1] = (unsigned int) floor( (position[1]-bounding_box[1])/dy ); 
	out[2] = (unsigned int) floor( (position[2]-bounding_box[2])/dz ); 

	//  add some bounds checking -- truncate to inside the computational domain  

	if( out[0] >= x_coordinates.size() ){ out[0] = x_coordinates.size()-1; }
	if( out[0] < 0 ){ out[0] = 0; }

	if( out[1] >= y_coordinates.size() ){ out[1] = y_coordinates.size()-1; }
	if( out[1] < 0 ){ out[1] = 0; }

	if( out[2] >= z_coordinates.size() ){ out[2] = z_coordinates.size()-1; }
	if( out[2] < 0 ){ out[2] = 0; }

	return out; 
}

Voxel& Cartesian_Mesh::nearest_voxel( std::vector<double>& position )
{ return voxels[ nearest_voxel_index( position ) ]; }

void Cartesian_Mesh::display_information( std::ostream& os )
{
	os << std::endl << "Mesh information: " << std::endl;
	if( uniform_mesh ) 
	{ os << "type: uniform Cartesian" << std::endl; }
	else
	{
		if( regular_mesh )
		{ os << "type: regular Cartesian" << std::endl; }
		else
		{ os << "type: general Cartesian" << std::endl; }
	}
	os << "Domain: " 
	<< "[" << bounding_box[0] << "," << bounding_box[3] << "] " <<  units << " x " 
	<< "[" << bounding_box[1] << "," << bounding_box[4] << "] " <<  units << " x " 
	<< "[" << bounding_box[2] << "," << bounding_box[5] << "] " <<  units << std::endl
	<< "   resolution: dx = " << dx << " " << units; 
	if( !uniform_mesh )
	{
		os	<< ", dy = " << dy << " " << units 
			<< ", dz = " << dz << " " << units ; 
	}
	os << std::endl 
	<< "   voxels (per-process): " << voxels.size() << std::endl //--->These are voxels per-process now (Gaurav Saxena)
	<< "   voxel faces: " << voxel_faces.size() << std::endl
	<< "   volume: " << ( bounding_box[3]-bounding_box[0] )*( bounding_box[4]-bounding_box[1] )*( bounding_box[5]-bounding_box[2] ) 
		<< " cubic " << units << std::endl; 	

	return; 
}

void Cartesian_Mesh::read_from_matlab( std::string filename )
{
	unsigned int size_of_each_datum; 
	unsigned int number_of_data_entries; 
	FILE* fp = read_matlab_header( &size_of_each_datum, &number_of_data_entries,  filename ); 

	voxel_faces.resize( 0 ); 
	
	connected_voxel_indices.resize( 1 ); 
	connected_voxel_indices[0].clear(); 
		
	Cartesian_mesh = true; 
	uniform_mesh = false; 
	regular_mesh = true; 
	use_voxel_faces = false; 

	// resize the internal data structure 

	voxels.resize( number_of_data_entries );
	connected_voxel_indices.resize( voxels.size() ); 
	
	x_coordinates.resize( number_of_data_entries );
	y_coordinates.resize( number_of_data_entries );
	z_coordinates.resize( number_of_data_entries );
	
	// read in the data
	// assumes each column has: x,y,z, dV
	
	bounding_box[0] = 9e99; 
	bounding_box[1] = 9e99; 
	bounding_box[2] = 9e99; 

	bounding_box[3] = -9e99; 
	bounding_box[4] = -9e99; 
	bounding_box[5] = -9e99; 
 
    size_t result;
	for( unsigned int i=0; i < number_of_data_entries ; i++ )
	{
		result = fread( (char*) &( voxels[i].center[0] ) , sizeof(double) , 1 , fp ); 
		result = fread( (char*) &( voxels[i].center[1] ) , sizeof(double) , 1 , fp ); 
		result = fread( (char*) &( voxels[i].center[2] ) , sizeof(double) , 1 , fp ); 
		result = fread( (char*) &( voxels[i].volume ) , sizeof(double) , 1 , fp ); 
		
		// estimate the bounding box; 
		if( voxels[i].center[0] < bounding_box[0] )
		{ bounding_box[0] = voxels[i].center[0]; }
		if( voxels[i].center[0] > bounding_box[3] )
		{ bounding_box[3] = voxels[i].center[0]; }

		if( voxels[i].center[1] < bounding_box[1] )
		{ bounding_box[1] = voxels[i].center[1]; }
		if( voxels[i].center[1] > bounding_box[4] )
		{ bounding_box[4] = voxels[i].center[1]; }

		if( voxels[i].center[2] < bounding_box[2] )
		{ bounding_box[2] = voxels[i].center[2]; }
		if( voxels[i].center[2] > bounding_box[5] )
		{ bounding_box[5] = voxels[i].center[2]; }
	} 
	
	// figure out dx, dy, dz 

	double xmin = bounding_box[0]; // voxels[0].center[0]; 
	double ymin = bounding_box[1]; // voxels[0].center[1]; 
	double zmin = bounding_box[2]; // voxels[0].center[2]; 

	// int n = voxels.size(); 
	double xmax = bounding_box[3]; // voxels[n-1].center[0]; 
	double ymax = bounding_box[4]; // voxels[n-1].center[1]; 
	double zmax = bounding_box[5]; // voxels[n-1].center[2]; 

	// figure out number of x nodes  
	int xnodes = 0; 
	while( fabs( voxels[xnodes].center[0] - xmax ) > 1e-15 )
	{ xnodes++; }
	xnodes++; 

	// figure out number of y nodes 
	int ynodes = 0; 

	while( fabs( voxels[ynodes*xnodes].center[1] - ymax ) > 1e-15 )
	{ ynodes += 1; }
	ynodes++;

	// figure out number of z nodes 

	int znodes = 0; 

	while( fabs( voxels[ynodes*xnodes*znodes].center[2] - zmax ) > 1e-15 )
	{ znodes += 1; }
	znodes++;

	// figure out differentials

	dx = ( xmax - xmin ) / ( (double) xnodes - 1.0  ); 
	dy = ( ymax - ymin ) / ( (double) ynodes - 1.0  ); 
	dz = ( zmax - zmin ) / ( (double) znodes - 1.0  ); 
	
	dV = dx*dy*dz; 
	dS = dx*dy;

	dS_xy = dx*dy;
	dS_yz = dy*dz; 
	dS_xz = dx*dz;

	uniform_mesh = true; 
	double tol = 1e-16; 
	if( fabs( dx - dy ) > tol || fabs( dy - dz ) > tol || fabs( dx - dz ) > tol )
	{ uniform_mesh = false; }
	
	// correct the bounding box 
	
	double half_step = dx * 0.5; 
	
	bounding_box[0] -= half_step; 
	bounding_box[3] += half_step;
	
	half_step = dy * 0.5; 
	bounding_box[1] -= half_step; 
	bounding_box[4] += half_step;
	
	half_step = dz * 0.5; 
	bounding_box[2] -= half_step; 
	bounding_box[5] += half_step;
 
	// write out the x,y,z coordinates; 
	x_coordinates.resize( xnodes ); 
	y_coordinates.resize( ynodes ); 
	z_coordinates.resize( znodes ); 

	for( unsigned int i=0; i < x_coordinates.size() ; i++ )
	{ x_coordinates[i] = xmin + i*dx ;   }

	for( unsigned int i=0; i < y_coordinates.size() ; i++ )
	{ y_coordinates[i] = ymin + i*dy ; }

	for( unsigned int i=0; i < z_coordinates.size() ; i++ )
	{ z_coordinates[i] = zmin + i*dz ; }

	dV = dx*dy*dz; 
	
	units = "none";  
	
	// still need to figure out connected indices later 
	std::cout << "Warning: Cartesian_Mesh::read_from_matlab is incomplete. No connection information read." << std::endl; 

	fclose( fp) ; 
	return; 
}

void Voronoi_Mesh::display_information( std::ostream& os )
{
	os << std::endl << "Mesh information: " << std::endl 
	<< "type: Voronoi (not implemented!)" << std::endl
	<< "Domain: " 
	<< "[" << bounding_box[0] << "," << bounding_box[3] << "] " <<  units << " x " 
	<< "[" << bounding_box[1] << "," << bounding_box[4] << "] " <<  units << " x " 
	<< "[" << bounding_box[2] << "," << bounding_box[5] << "] " <<  units << std::endl
	<< "   voxels: " << voxels.size() << std::endl
	<< "   voxel faces: " << voxel_faces.size() << std::endl 
	<< "   volume: "; 

	double total_volume = 0.0; 
	for( unsigned int i=0; i < voxels.size(); i++ )
	{ total_volume += voxels[i].volume; }
	os << total_volume << " cubic " << units << std::endl; 

	return; 
}

};
