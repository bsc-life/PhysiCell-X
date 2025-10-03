/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2025, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./PhysiCell_SVG.h"

bool Write_SVG_start( std::ostream& os, double width, double height )
{
 os << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" << std::endl 
    << "<!-- Created with PhysiCell (http://PhysiCell.MathCancer.org/) -->" << std::endl; 

 os << "<svg " << std::endl
    << " xmlns:dc=\"http://purl.org/dc/elements/1.1/\" " << std::endl
    << " xmlns:cc=\"http://creativecommons.org/ns#\" " << std::endl
    << " xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" " << std::endl
    << " xmlns:svg=\"http://www.w3.org/2000/svg\" " << std::endl
    << " xmlns=\"http://www.w3.org/2000/svg\" " << std::endl
    << " version=\"1.1\" " << std::endl
    << " width=\"" << width << "\" " << std::endl
    << " height=\"" << height << "\" " << std::endl
    << " id=\"svg2\">" << std::endl;
	
	return true; 
}

/*----------------------------------------------------*/
/* Parallel version of the function Write_SVG_start() */
/*----------------------------------------------------*/

bool Write_SVG_start( std::string & file_str, double width, double height, mpi_Environment &world, mpi_Cartesian &cart_topo )
{
 
 int count = 16;
 std::string str[count]; 
 
 str[0] = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"; 				
 str[1] = "<!-- Created with PhysiCell (http://PhysiCell.MathCancer.org/) -->\n"; 

 str[2] = "<svg \n"; //5+1=6 bytes
 str[3] = " xmlns:dc=\"http://purl.org/dc/elements/1.1/\" \n";	
 str[4] = " xmlns:cc=\"http://creativecommons.org/ns#\" \n"; 		
 str[5] = " xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" \n";	
 str[6] = " xmlns:svg=\"http://www.w3.org/2000/svg\" \n";	
 str[7] = " xmlns=\"http://www.w3.org/2000/svg\" \n";			
 str[8] = " version=\"1.1\" \n";	
 str[9] = " width=\""; 
 str[10]= std::to_string(width); 
 str[11]= "\" \n";	
 str[12]= " height=\"";
 str[13]= std::to_string(height);
 str[14]= "\" \n"; 
 str[15]= " id=\"svg2\">\n";	
 
 for(int i=0; i<count; i++)
 	file_str.append(str[i]); 
 
	
	return true; 
}

bool Write_SVG_end( std::ostream& os )
{
 os << "</svg>" << std::endl;
 return true; 
}

/*----------------------------------------------------*/
/* Parallel version of the function Write_SVG_end()  */
/*----------------------------------------------------*/

bool Write_SVG_end( std::string &file_str, mpi_Environment &world, mpi_Cartesian &cart_topo )
{
 std::string str = "</svg>\n";
 file_str.append(str); 
 return true;
}

bool Write_SVG_text( std::ostream& os, const char* str , double position_x, double position_y, double font_size , const char* color , const char* font)
{
 os << "  <text x=\"" << position_x << "\" y=\""  << position_y << "\"" << std::endl
    << "   font-family=\"" << font << "\" font-size=\"" << font_size << "\" fill=\"" << color << "\" >" << std::endl
    << "   " << str << std::endl << "  </text>" << std::endl; 
  return true; 
}

/*----------------------------------------------------*/
/* Parallel version of the function Write_SVG_text()  */
/*----------------------------------------------------*/

bool Write_SVG_text( std::string &file_str, const char* str , double position_x, double position_y, double font_size , const char* color , const char* font, mpi_Environment &world, mpi_Cartesian &cart_topo)
{
  int count = 16;
  std::string str_arr[count];
  
  str_arr[0] = "  <text x=\"";
  str_arr[1] = std::to_string(position_x);
  str_arr[2] = "\" y=\"";
  str_arr[3] = std::to_string(position_y);
  str_arr[4] = "\"\n";
  str_arr[5] = "   font-family=\"";
  str_arr[6] = font;
  str_arr[7] = "\" font-size=\""; 
  str_arr[8] = std::to_string(font_size);
  str_arr[9] = "\" fill=\"";
  str_arr[10]= color;
  str_arr[11]= "\" >\n";
  str_arr[12]= "   ";
  str_arr[13]= str;
  str_arr[14]= "\n";
  str_arr[15]= "  </text>\n";
  
  for(int i=0; i<count; i++)
  	file_str.append(str_arr[i]);
  
   
  return true; 
}


bool Write_SVG_circle( std::ostream& os, double center_x, double center_y, double radius, double stroke_size, 
                       std::string stroke_color , std::string fill_color )
{
 os << "  <circle cx=\"" << center_x << "\" cy=\"" << center_y << "\" r=\"" << radius << "\" stroke-width=\"" << stroke_size 
    << "\" stroke=\"" << stroke_color << "\" fill=\"" << fill_color << "\"/>" << std::endl; 
 return true; 
}

/*----------------------------------------*/
/* Parallel version of Write_SVG_circle() */
/*----------------------------------------*/


bool Write_SVG_circle( std::string& file_str, double center_x, double center_y, double radius, double stroke_size, std::string stroke_color , std::string fill_color, mpi_Environment &world, mpi_Cartesian &cart_topo )
{
		int count = 13; 
		std::string str_arr[count];
		
		str_arr[0] = "  <circle cx=\"";
		str_arr[1] = std::to_string(center_x);
		str_arr[2] = "\" cy=\"";
		str_arr[3] = std::to_string(center_y);
		str_arr[4] = "\" r=\"";
		str_arr[5] = std::to_string(radius);
		str_arr[6] = "\" stroke-width=\"";
		str_arr[7] = std::to_string(stroke_size);
		str_arr[8] =  "\" stroke=\"";
		str_arr[9] = stroke_color;
		str_arr[10]= "\" fill=\"";
		str_arr[11]= fill_color;
		str_arr[12]= "\"/>\n"; 
		
		for(int i=0; i<count; i++)
			file_str.append(str_arr[i]); 
 return true; 
}



bool Write_SVG_rect( std::ostream& os , double UL_corner_x, double UL_corner_y, double width, double height, 
                     double stroke_size, std::string stroke_color , std::string fill_color )
{
 os << "  <rect x=\"" << UL_corner_x << "\" y=\"" << UL_corner_y << "\" width=\"" << width << "\" height=\"" 
    << height << "\" stroke-width=\"" << stroke_size 
    << "\" stroke=\"" << stroke_color << "\" fill=\"" << fill_color << "\"/>" << std::endl; 
 return true; 
}

/*----------------------------------------------------*/
/* Parallel version of the function Write_SVG_rect()  */
/*----------------------------------------------------*/

bool Write_SVG_rect( std::string & file_str , double UL_corner_x, double UL_corner_y, double width, double height, 
                     double stroke_size, std::string stroke_color , std::string fill_color, mpi_Environment &world, mpi_Cartesian &cart_topo  )
{
  int count = 15;
  std::string str[count]; 
  
  str[0] = "  <rect x=\"";
  str[1] = std::to_string(UL_corner_x);
  str[2] = "\" y=\"";
  str[3] = std::to_string(UL_corner_y);
  str[4] = "\" width=\"";
  str[5] = std::to_string(width);
  str[6] = "\" height=\""; 
  str[7] = std::to_string(height);
  str[8] = "\" stroke-width=\"";
  str[9] = std::to_string(stroke_size); 
  str[10]= "\" stroke=\"";
  str[11]= stroke_color;
  str[12]= "\" fill=\"";
  str[13]= fill_color;
  str[14]= "\"/>\n";
  
  for(int i=0; i<count; i++)
  	file_str.append(str[i]); 
   
 return true; 
}

bool Write_SVG_line( std::ostream& os , double start_x, double start_y, double end_x , double end_y, double thickness, 
                    std::string stroke_color )
{
 os << "  <line x1=\"" << start_x << "\" y1=\"" << start_y << "\" x2=\"" << end_x << "\" y2=\"" << end_y << "\" "
    << "stroke=\"" << stroke_color << "\" stroke-width=\"" << thickness << "\"/>" << std::endl; 
 return true; 
}

void Write_SVG_text(std::ostream& os, const char* str , double position_x, double position_y, double font_size , const char* color , const char* font, double rotation)
{
    double text_width = font_size * strlen(str) / 2.0;  // estimate the width of the text
    double text_height = font_size / 2.0;  // estimate the height of the text

    double center_x = position_x + text_width / 2.0;
    double center_y = position_y + text_height / 2.0;

    os << "<text x=\"" << position_x << "\" y=\"" << position_y << "\" font-size=\"" << font_size << "\" fill=\"" << color << "\" font-family=\"" << font << "\" transform=\"rotate(" << rotation << " " << center_x << " " << center_y << ")\">" << str << "</text>\n";
}

void Write_SVG_text(std::string & file_str, const char* str , double position_x, double position_y, double font_size , const char* color , const char* font, double rotation)
{
    double text_width = font_size * strlen(str) / 2.0;  // estimate the width of the text
    double text_height = font_size / 2.0;  // estimate the height of the text

    double center_x = position_x + text_width / 2.0;
    double center_y = position_y + text_height / 2.0;

    int count = 19;

    std::string aux[count];

    aux[0] = "<text x=\"";
    aux[1] = std::to_string(position_x);
    aux[2] = "\" y=\"";
    aux[3] = std::to_string(position_y);
    aux[4] = "\" font-size=\"";
    aux[5] = std::to_string(font_size);
    aux[6] = "\" fill=\"";
    aux[7] = color;
    aux[8] = "\" font-family=\"";
    aux[9] = font;
    aux[10] = "\" transform=\"rotate(";
    aux[11] = std::to_string(rotation);
    aux[12] = " ";
    aux[13] = std::to_string(center_x);
    aux[14] = " ";
    aux[15] = std::to_string(center_y);
    aux[16] = ")\">";
    aux[17] = str;
    aux[18] = "</text>\n";

    for(int i=0; i<count; i++)
  	   file_str.append(aux[i]); 

   }

