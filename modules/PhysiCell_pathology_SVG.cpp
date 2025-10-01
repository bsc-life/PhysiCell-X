
#include "./PhysiCell_pathology.h"

#include "./PhysiCell_SVG.h"

namespace PhysiCell{





void SVG_plot_mpi( std::string filename , Microenvironment& M, double z_slice , double time, std::vector<std::string> (*cell_coloring_function)(Cell*), std::string (*substrate_coloring_function)(double, double, double), mpi_Environment &world, mpi_Cartesian &cart_topo )
{
	MPI_File fh;
	MPI_Offset offset; 
	MPI_Datatype etype, filetype;
	etype = MPI_CHAR;
	filetype = MPI_CHAR;
	char char_filename[filename.size()+1];
	int elements_to_write;
	int mpi_error;
	int global_total_cell_count; 

	std::string file_str=""; 											//This string would contain EVERYTHING that a process has to write  
	int *sum 				= new int[world.size]; 							//Root process will receive individual file_str lengths
	int *cum_sum 		= new int[world.size];							//Array containing cumulative sum
	int local_chars = 0;								//C strings take 1 more 
	char *data 			= new char[1]; 			//C character array for MPI-IO
	int char_offset_of_process;													//Dont want to catch cum_sum in MPI_Offset type
 
	
	double X_lower = M.mesh.bounding_box[0];
	double X_upper = M.mesh.bounding_box[3];
 
	double Y_lower = M.mesh.bounding_box[1]; 
	double Y_upper = M.mesh.bounding_box[4]; 

	double plot_width = X_upper - X_lower; 
	double plot_height = Y_upper - Y_lower; 

	double font_size = 0.025 * plot_height; 				// PhysiCell_SVG_options.font_size; 
	double top_margin = font_size*(.2+1+.2+.9+.5 ); 


	/*-----------------------------------------------------------------*/
	/* Since we need to display a message if the file cannot be opened */
	/* I need to open the file here but defer the writing into the file*/
	/* till we have the complete strings to be written into the file 	 */
	/*-----------------------------------------------------------------*/


	/* open the file, BUT defer the "header" writing till ALL "file_str" strings are complete */
	
	
	strcpy(char_filename, filename.c_str()); 
	mpi_error = MPI_File_open(cart_topo.mpi_cart_comm, char_filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	
	MPI_Barrier(cart_topo.mpi_cart_comm);  // Gaurav Saxena is adding this to synchronize 
	
	if( mpi_error)
	{ 
		if(world.rank == 0)
		{
			std::cout << std::endl << "Error: Failed to open " << filename << " for SVG writing." << std::endl << std::endl; 

			std::cout << std::endl << "Error: We're not writing data like we expect. " << std::endl
								<< "Check to make sure your save directory exists. " << std::endl << std::endl
								<< "I'm going to exit with a crash code of -1 now until " << std::endl 
								<< "you fix your directory. Sorry!" << std::endl << std::endl; 
			exit(-1);
		} 
	} 
	
	/*----------------------------------------------------------------------------------*/
	/* Now all the writing functions TILL the "for" loop are called by the root process */
	/*----------------------------------------------------------------------------------*/
	
	if(world.rank == 0) {
		if(PhysiCell_settings.enable_substrate_plot == true && (*substrate_coloring_function) != NULL){

			double legend_padding = 200.0; // I have to add a margin on the left to visualize the bar plot and the values

			Write_SVG_start( file_str, plot_width + legend_padding, plot_height + top_margin, world, cart_topo );

		// draw the background
			Write_SVG_rect( file_str , 0 , 0 , plot_width + legend_padding, plot_height + top_margin , 0.002 * plot_height , "white", "white",world, cart_topo );

		}
		else{

			Write_SVG_start( file_str, plot_width , plot_height + top_margin, world, cart_topo );

			// draw the background
			Write_SVG_rect( file_str , 0 , 0 , plot_width, plot_height + top_margin , 0.002 * plot_height , "white", "white", world, cart_topo );

		}
	}

	// Write the simulation time to the top of the plot
 
	char* szString; 
	szString = new char [1024]; 
 
	/* Need to total the cells on each process and store global sum in global_total_cell_count */
	 
	int total_cell_count = all_cells->size();
	MPI_Reduce(&total_cell_count, &global_total_cell_count, 1, MPI_INT, MPI_SUM, 0, cart_topo.mpi_cart_comm ); 
 
	double temp_time = time; 

	std::string time_label = formatted_minutes_to_DDHHMM( temp_time ); 
 
	sprintf( szString , "Current time: %s, z = %3.2f %s", time_label.c_str(), 
		z_slice , PhysiCell_SVG_options.simulation_space_units.c_str() ); 
	
	if(world.rank == 0)
		Write_SVG_text( file_str, szString, font_size*0.5,  font_size*(.2+1), 
										font_size, PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str(), world, cart_topo );
		
	sprintf( szString , "%u agents" , global_total_cell_count ); 
	
	if(world.rank == 0)
		Write_SVG_text( file_str, szString, font_size*0.5,  font_size*(.2+1+.2+.9), 
										0.95*font_size, PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str(), world, cart_topo );
	
	delete [] szString; 

	// Add an outer "g" for coordinate transforms 
	
	if(world.rank == 0)
	{
		int count = 4;
		std::string str1[count];
		
	  str1[0] = " <g id=\"tissue\" \n";
	  str1[1] = "    transform=\"translate(0,";
	  str1[2] = std::to_string(plot_height+top_margin);
	  str1[3] = ") scale(1,-1)\">\n";
	  
	  for(int i=0; i<count; i++)
	  	file_str.append(str1[i]);
	} 
	   
	// Prepare to do mesh-based plot (later)
	
	double dx_stroma = M.mesh.dx; 
	double dy_stroma = M.mesh.dy; 
	
	if(world.rank == 0)
	{
		std::string str2;
		str2 = "  <g id=\"ECM\">\n"; 
		file_str.append(str2);
	}

  
	int ratio = 1; 
	double voxel_size = dx_stroma / (double) ratio ; 
  
	double half_voxel_size = voxel_size / 2.0; 
	double normalizer = 78.539816339744831 / (voxel_size*voxel_size*voxel_size); 
 
	//Include substrates in the SVG
	double max_conc_local, max_conc_global;
	double min_conc_local, min_conc_global;
	// color in the background ECM
	if(PhysiCell_settings.enable_substrate_plot == true && (*substrate_coloring_function) != NULL)
	{
		double dz_stroma = M.mesh.dz;

		std::string sub = PhysiCell_settings.substrate_to_monitor;
		int sub_index = M.find_density_index(sub); // check the substrate does actually exist
		if(sub_index == -1){
			if (world.rank == 0)
				std::cout << "ERROR SAMPLING THE SUBSTRATE: COULD NOT FIND THE SUBSTRATE " << sub << std::endl; //if not print error message
		}
		else
		{
			if(PhysiCell_settings.limits_substrate_plot){
			 max_conc_local = PhysiCell_settings.max_concentration;
			 min_conc_local = PhysiCell_settings.min_concentration;
			}
			else{
			 max_conc_local = M.density_vector(5)[sub_index];
			 min_conc_local = M.density_vector(5)[sub_index];	 // so here I am sampling the concentration to set a min and a mx
			//look for the max and min concentration among all the substrates
			for (int n = 0; n < M.number_of_voxels(); n++)
			{
				double concentration = M.density_vector(n)[sub_index];
				if (concentration > max_conc_local)
					max_conc_local = concentration;
				if (concentration < min_conc_local)
					min_conc_local = concentration;
			}
			//MPI allgather the max and the min from all to all
			MPI_Allreduce(&max_conc_local, &max_conc_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			MPI_Allreduce(&min_conc_local, &min_conc_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
			};

			//check that max conc is not zero otherwise it is a big problem!
			if(max_conc_global == 0){

				max_conc_global = 1.0;

			};

			for (int n = 0; n < M.number_of_voxels(); n++)
			{
				auto current_voxel = M.voxels(n);
				int z_center = current_voxel.center[2];
				double z_displ = z_center -  dz_stroma/2; 

				double z_compare = z_displ;

				if (default_microenvironment_options.simulate_2D == true){
					z_compare = z_center;
				};

				if (z_slice == z_compare){			//this is to make sure the substrate is sampled in the voxel visualized (so basically the slice)
					int x_center = current_voxel.center[0];
					int y_center = current_voxel.center[1];
					
					double x_displ = x_center -  dx_stroma/2;
					double y_displ = (y_center - dy_stroma) +  dy_stroma/2;

					double concentration = M.density_vector(n)[sub_index];

					std::string output = substrate_coloring_function(concentration, max_conc_global, min_conc_global );
					Write_SVG_rect( file_str , x_displ - X_lower , y_displ - Y_lower, dx_stroma, dy_stroma , 0 , "none", output, world, cart_topo );
				}
			}
		}
	}

	if(world.rank == world.size -1)
	{
		std::string str3 = "  </g>\n";
		file_str.append(str3); 
 	}

    //print header and ecm and clear all
	data = new char[file_str.size()];
    std::memcpy(data, file_str.data(), file_str.size());
	
	local_chars = file_str.size();
	MPI_Allgather(&local_chars, 1, MPI_INT, sum, 1, MPI_INT, cart_topo.mpi_cart_comm);
		
	
	for(int i=0; i<world.size; i++)
		cum_sum[i] = 0;
		
	for(int i=1; i<world.size; i++)
		cum_sum[i] = sum[i-1] + cum_sum[i-1]; 
	
	
	char_offset_of_process = 0;
	MPI_Scatter(cum_sum, 1, MPI_INT, &char_offset_of_process, 1, MPI_INT, 0, cart_topo.mpi_cart_comm);

	// compute MPI_Offset safely
	offset = static_cast<MPI_Offset>(char_offset_of_process);
	//std::cout << "\t\tRank " << world.rank << " offset :" << offset << " to write " << offset + local_chars << std::endl; 
    MPI_File_set_view(fh, offset, etype, filetype, "native", MPI_INFO_NULL); 
    MPI_File_write(fh, data, local_chars, MPI_CHAR, MPI_STATUS_IGNORE);

    file_str.clear();

	


	// plot intersecting cells
	
	if(world.rank == 0)
	{ 
		std::string str4 = "  <g id=\"cells\">\n";
		file_str.append(str4);
	}
	
    

	std::string str5[4];
	
	for( int i=0 ; i < total_cell_count ; i++ )
	{
		Cell* pC = (*all_cells)[i]; // global_cell_list[i]; 
  
		static std::vector<std::string> Colors; 
		if( fabs( (pC->position)[2] - z_slice ) < pC->phenotype.geometry.radius )
		{
			double r = pC->phenotype.geometry.radius ; 
			double rn = pC->phenotype.geometry.nuclear_radius ; 
			double z = fabs( (pC->position)[2] - z_slice) ; 
   
			Colors = cell_coloring_function( pC ); 

			str5[0] = "   <g id=\"cell";
			str5[1] = std::to_string(pC->ID); 
			str5[2] = "\">\n";
			
			for(int ctr=0; ctr<3; ctr++)
				file_str.append(str5[ctr]); 
  
			// figure out how much of the cell intersects with z = 0 
   
			double plot_radius = sqrt( r*r - z*z ); 

			Write_SVG_circle( file_str, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower, 
				plot_radius , 0.5, Colors[1], Colors[0], world, cart_topo ); 

			// plot the nucleus if it, too intersects z = 0;
			if( fabs(z) < rn && PhysiCell_SVG_options.plot_nuclei == true )
			{   
				plot_radius = sqrt( rn*rn - z*z ); 
			 	Write_SVG_circle( file_str, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower, 
					plot_radius, 0.5, Colors[3],Colors[2], world, cart_topo); 
			}					  
			str5[3] = "   </g>\n";
			
			file_str.append(str5[3]);
		}
	}
	
	/* After the for loop, all writing is to be done by rank = world.size-1 */
	
	if(world.rank == world.size-1)
	{
		std::string str6 = "  </g>\n"; 
		file_str.append(str6);
	}
    //print here
	// end of the <g ID="tissue">
	
	if(world.rank == world.size-1)
	{
		std::string str7 = " </g>\n"; 
		file_str.append(str7); 
 	}
 	
	// draw a scale bar
 
	double bar_margin = 0.025 * plot_height; 
	double bar_height = 0.01 * plot_height; 
	double bar_width = PhysiCell_SVG_options.length_bar; 
	double bar_stroke_width = 0.001 * plot_height; 
	
	std::string bar_units = PhysiCell_SVG_options.simulation_space_units; 
	
	// Convert from micron to mm
	
	double temp = bar_width;  

	if( temp > 999 && std::strstr( bar_units.c_str() , PhysiCell_SVG_options.mu.c_str() )   )
	{
		temp /= 1000;
		bar_units = "mm";
	}
	// Convert from mm to cm 
	
	if( temp > 9 && std::strcmp( bar_units.c_str() , "mm" ) == 0 )
	{
		temp /= 10; 
		bar_units = "cm";
	}
	
	szString = new char [1024];
	sprintf( szString , "%u %s" , (int) round( temp ) , bar_units.c_str() );
  
  if(world.rank == world.size-1)
  {
		Write_SVG_rect( file_str , plot_width - bar_margin - bar_width  , plot_height + top_margin - bar_margin - bar_height , 
										bar_width , bar_height , 0.002 * plot_height , "rgb(255,255,255)", "rgb(0,0,0)", world, cart_topo );
		Write_SVG_text( file_str, szString , plot_width - bar_margin - bar_width + 0.25*font_size , 
										plot_height + top_margin - bar_margin - bar_height - 0.25*font_size , 
										font_size , PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str(), world, cart_topo ); 
	}
	
	delete [] szString; 

	// Plot runtime 
	
	szString = new char [1024]; 
	RUNTIME_TOC(); 
	std::string formatted_stopwatch_value = format_stopwatch_value( runtime_stopwatch_value() );
	
	if(world.rank == world.size-1)
		Write_SVG_text( file_str, formatted_stopwatch_value.c_str() , bar_margin , top_margin + plot_height - bar_margin , 0.75 * font_size , 
										PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str(), world, cart_topo );
	delete [] szString; 

	//Draw legend for susbtrates


	// Draw a box around the plot window
	
	if(world.rank == world.size-1)
		Write_SVG_rect( file_str , 0 , top_margin, plot_width, plot_height , 0.002 * plot_height , "rgb(0,0,0)", "none", world, cart_topo );
	
	//
	if(world.rank == world.size-1) {
		if (PhysiCell_settings.enable_substrate_plot == true && (*substrate_coloring_function) != NULL) {

			// add legend for the substrate

			double conc_interval = (max_conc_global - min_conc_global) / 10; // setting the interval for the values in the legend.

			szString = new char [1024];
			double upper_left_x = plot_width + 25.0;
			double sub_rect_height = (plot_height - 25.0) / 10.0;
			for(int i = 0; i <= 9; i++){ //creating 10 rectangoles for the bar, each one with a different shade of color.

				double concentration_sample = min_conc_global + (conc_interval * (9-i)); // the color depends on the concentration, starting from the min concentration to the max (which was sampled before)

				std::string output = substrate_coloring_function(concentration_sample, max_conc_global, min_conc_global);

				double upper_left_y = sub_rect_height * i; // here I set the position of each rectangole

				Write_SVG_rect(file_str, upper_left_x, top_margin + upper_left_y, 25.0, sub_rect_height, 0.002 * plot_height , "none", output, world, cart_topo); //drawing each piece of the barplot

				if(i%2 != 0){ // of course I am not printing each value of the barplot, otherwise is too crowded, so just one each 2

					sprintf( szString , " %.2g", concentration_sample);
					Write_SVG_rect(file_str, upper_left_x + 25, top_margin + upper_left_y + sub_rect_height - (0.001 * plot_height), 3, 0.002 * plot_height, 0 , "rgb(0,0,0)", "rgb(0,0,0)", world, cart_topo);
					Write_SVG_text( file_str , szString, upper_left_x + 28, top_margin + upper_left_y + sub_rect_height, font_size ,
									PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str(), world, cart_topo ); // misterious values set with a trial and error approach due to OCD. But now the legend is coherent at pixel level
				}
			}

			sprintf( szString , "%.2g", max_conc_global);

			Write_SVG_rect(file_str, upper_left_x + 25, top_margin - (0.001 * plot_height), 3, 0.002 * plot_height, 0 , "rgb(0,0,0)", "rgb(0,0,0)", world, cart_topo);
			Write_SVG_text( file_str , szString, upper_left_x + 28, top_margin, font_size , 
				PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str(),world, cart_topo ); // misterious values set with a trial and error approach due to OCD. But now the legend is coherent at pixel level

			delete [] szString;

			// add a label to the right of the colorbar defined by above Write_SVG_rect calls
			Write_SVG_text(file_str, PhysiCell_settings.substrate_to_monitor.c_str(), upper_left_x + 35, top_margin + plot_height / 2, font_size,
							PhysiCell_SVG_options.font_color.c_str(), PhysiCell_SVG_options.font.c_str(), 90.0);
		}
	}
	
	// Close the svg tag
	
	if(world.rank == world.size-1)
		Write_SVG_end(file_str, world, cart_topo); 
		

	//last rank tell to first rank where it ended
	int offset_ecm;
	if (world.rank == world.size -1) {
		offset_ecm = offset + local_chars;
		MPI_Send(&offset_ecm, 1, MPI_INT, 0, 2, cart_topo.mpi_cart_comm);
	}
	else if (world.rank == 0) {
		MPI_Recv(&offset_ecm, 1, MPI_INT, world.rank -1, 2, cart_topo.mpi_cart_comm, MPI_STATUS_IGNORE);
	} 

	local_chars = file_str.size();
	data = new char[file_str.size()];
	if(local_chars > 0) memcpy(data, file_str.data(), local_chars);
	
	MPI_Gather(&local_chars, 1, MPI_INT, sum, 1, MPI_INT, 0, cart_topo.mpi_cart_comm);
		
	if(world.rank == 0)
	{
		cum_sum[0] = offset_ecm;
			
		for(int i=1; i<world.size; i++)
			cum_sum[i] = sum[i-1] + cum_sum[i-1]; 
	}
	
	MPI_Scatter(cum_sum, 1, MPI_INT, &char_offset_of_process, 1, MPI_INT, 0, cart_topo.mpi_cart_comm);
	
	offset = char_offset_of_process * sizeof(char); 
    //td::cout << "\t\tRank " << world.rank << " offset :" << offset << " to write " << offset + local_chars << std::endl; 

	MPI_File_set_view(fh, offset, etype, filetype, "native", MPI_INFO_NULL); 
	MPI_File_write(fh, data, local_chars, MPI_CHAR, MPI_STATUS_IGNORE);
  
	MPI_Barrier(cart_topo.mpi_cart_comm);  // Gaurav Saxena is adding this to synchronize 
    
	MPI_File_close(&fh);
  
	delete [] data;
	if(world.rank == 0) { delete [] sum; delete [] cum_sum; }
	file_str.clear(); 

	return; 
}

}