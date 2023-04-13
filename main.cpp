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
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <vector>
#include <string>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h" 

// custom user modules 

#include "./custom_modules/cancer_immune_3D.h" 
	
using namespace BioFVM;
using namespace PhysiCell;

int main( int argc, char* argv[] )
{
	// load and parse settings file(s)
	
	bool XML_status = false; 
	char copy_command [1024]; 
	if( argc > 1 )
	{
		XML_status = load_PhysiCell_config_file( argv[1] ); 
		sprintf( copy_command , "cp %s %s" , argv[1] , PhysiCell_settings.folder.c_str() ); 
	}
	else
	{
		XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml" );
		sprintf( copy_command , "cp ./config/PhysiCell_settings.xml %s" , PhysiCell_settings.folder.c_str() ); 
	}
	if( !XML_status )
	{ exit(-1); }
	
	// copy config file to output directry 
	system( copy_command ); 
	
	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);
	
	// PNRG setup 
	SeedRandom(); 
	
	// time setup 
	std::string time_units = "min"; 

	/* Microenvironment setup */ 
	
	setup_microenvironment(); 

	/* PhysiCell setup */ 
 	
	// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 30; 
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );
	
	create_cell_types();
	setup_tissue();
	
	/* Users typically start modifying here. START USERMODS */ 

	double il13_immune_activation_time = 
		parameters.doubles("il13_immune_activation_time"); // 60 * 24 * 14; // activate immune response at 14 days 

	double her2_immune_activation_time =
		parameters.doubles("her2_immune_activation_time");

	double egfr_immune_activation_time = parameters.doubles("egfr_immune_activation_time");

	double number_of_il13_immune_cells = parameters.ints("number_of_il13_immune_cells");
	double number_of_her2_immune_cells = parameters.ints("number_of_her2_immune_cells");

	double number_of_egfr_immune_cells = parameters.ints("number_of_egfr_immune_cells");

	
	/* Users typically stop modifying here. END USERMODS */ 
	
	// set MultiCellDS save options 

	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true ); 
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true );
	
	// save a simulation snapshot 

	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	// save a quick SVG cross section through z = 0, after setting its 
	// length bar to 200 microns 

	PhysiCell_SVG_options.length_bar = 200; 

	// for simplicity, set a pathology coloring function 
	
	std::vector<std::string> (*cell_coloring_function)(Cell*) = cancer_immune_coloring_function;
	
	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
	
	display_citations(); 
	
	// set the performance timers 

	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	
	std::ofstream report_file;
	if( PhysiCell_settings.enable_legacy_saves == true )
	{	
		sprintf( filename , "%s/simulation_report.txt" , PhysiCell_settings.folder.c_str() ); 
		
		report_file.open(filename); 	// create the data log file 
		report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;
	}
	
	// main loop 
	std::vector<std::vector<std::string>> cells_count;
	try 
	{	
		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt )
		{
			static bool immune_cells_il13_introduced = false; 
			static bool immune_cells_her2_introduced = false;
			static bool immune_cells_egfr_introduced = false;
			
			if( PhysiCell_globals.current_time > il13_immune_activation_time - 0.01*diffusion_dt && immune_cells_il13_introduced == false && number_of_il13_immune_cells != 0)
			{
				std::cout << "IL13 Therapy activated!" << std::endl << std::endl; 
				immune_cells_il13_introduced = true; 
				
				PhysiCell_settings.full_save_interval = 
					parameters.doubles("save_interval_after_therapy_start"); 
				PhysiCell_settings.SVG_save_interval = 
					parameters.doubles("SVG_interval_after_therapy_start"); 
				
				PhysiCell_globals.next_full_save_time = PhysiCell_globals.current_time; 
				PhysiCell_globals.next_SVG_save_time = PhysiCell_globals.current_time; 
				
				introduce_il13_immune_cells();
				
			} 
			if (PhysiCell_globals.current_time > egfr_immune_activation_time - 0.01 * diffusion_dt && immune_cells_egfr_introduced == false && number_of_egfr_immune_cells != 0)
			{
				std::cout << "EGFR Therapy activated!" << std::endl << std::endl;
				immune_cells_egfr_introduced = true;

				PhysiCell_settings.full_save_interval =
					parameters.doubles("save_interval_after_therapy_start");
				PhysiCell_settings.SVG_save_interval =
					parameters.doubles("SVG_interval_after_therapy_start");

				PhysiCell_globals.next_full_save_time = PhysiCell_globals.current_time;
				PhysiCell_globals.next_SVG_save_time = PhysiCell_globals.current_time;

				introduce_egfr_immune_cells();

			}

			if (PhysiCell_globals.current_time > her2_immune_activation_time - 0.01 * diffusion_dt && immune_cells_her2_introduced == false && number_of_her2_immune_cells != 0)
			{
				std::cout << "HER2 Therapy activated!" << std::endl << std::endl;
				immune_cells_her2_introduced = true;

				PhysiCell_settings.full_save_interval =
					parameters.doubles("save_interval_after_therapy_start");
				PhysiCell_settings.SVG_save_interval =
					parameters.doubles("SVG_interval_after_therapy_start");

				PhysiCell_globals.next_full_save_time = PhysiCell_globals.current_time;
				PhysiCell_globals.next_SVG_save_time = PhysiCell_globals.current_time;

				introduce_her2_immune_cells();
			}

			// save data if it's time. 
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
			{
				display_simulation_status( std::cout ); 
				if( PhysiCell_settings.enable_legacy_saves == true )
				{	
					log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);
				}
				
				if( PhysiCell_settings.enable_full_saves == true )
				{	
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
					
					save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
				}
				
				PhysiCell_globals.full_output_index++; 
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
			}
			
			// save SVG plot if it's time
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 0.01 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{	
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index ); 
					//for (int i = 0; i < (*all_cells).size(); i++)
					//{
						//Cell* pC = (*all_cells)[i];
						//std::cout << "current phase: " << pC->phenotype.cycle.current_phase().code << std::endl;
					//}
					//current phase=14;
					//std::cout << "swelling: " << PhysiCell_constants::necrotic_swelling << std::endl;
					//std::cout << "necrotic_lysed: " << PhysiCell_constants::necrotic_lysed << std::endl;
					//std::cout << "necrotic: " << PhysiCell_constants::necrotic << std::endl;
					//necrotic_swelling=101; necrotic_lysed=102; necrotic=103;
					int cancer_cells_counts = 0;
					
					std::vector<std::string> current_cells_counts;
					int immune_egfr_counts = 0;
					int immune_il13_counts = 0;
					int immune_her2_counts = 0;
					int il13_cancer_counts = 0;

					for (int i = 0; i < (*all_cells).size(); i++)
					{
						
						Cell* pC = (*all_cells)[i];
						static int index_il13 = pC->custom_data.find_variable_index("IL13_antigen_level");

						if (pC->type == 0)
						{
							cancer_cells_counts += 1;
							if (pC->custom_data[index_il13] >= 30)
							{
								il13_cancer_counts += 1;
							}
						}
						else if (pC->type == 1) {
							immune_egfr_counts += 1;
							
						}
						else if (pC->type == 2) {
							immune_il13_counts += 1;
							
						}
						else {
							immune_her2_counts += 1;
							
						}
						
					}
					current_cells_counts.push_back(std::to_string(cancer_cells_counts));
					current_cells_counts.push_back(std::to_string(il13_cancer_counts));
					current_cells_counts.push_back(std::to_string(immune_egfr_counts));
					current_cells_counts.push_back(std::to_string(immune_il13_counts));
					current_cells_counts.push_back(std::to_string(immune_her2_counts));
					cells_count.push_back(current_cells_counts);
					/*std::cout << "number of cancer cells: " << current_cells_counts[0] << std::endl;
					std::cout << "number of EGFR immune cells: " << current_cells_counts[1] << std::endl;
					std::cout << "number of IL13 immune cells: " << current_cells_counts[2] << std::endl;
					std::cout << "number of HER2 immune cells: " << current_cells_counts[3] << std::endl;
					*/
					SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
					
					PhysiCell_globals.SVG_output_index++; 
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
					
				}
			}
			
			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt );
			// if( default_microenvironment_options.calculate_gradients )
			// { microenvironment.compute_all_gradient_vectors(); }
			
			// run PhysiCell 
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );
			//std::cout << PhysiCell_globals.current_time << std::endl;
			PhysiCell_globals.current_time += diffusion_dt;
		}
		
		if( PhysiCell_settings.enable_legacy_saves == true )
		{			
			log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
			report_file.close();
		}
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}
	
	// save a final simulation snapshot 
	
	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
	//std::cout << cells_count.size() << std::endl;
	std::fstream myfile;
	myfile.open("cells_count.txt", std::ios_base::out);
	
	for (int i = 0; i < cells_count.size(); i++)
	{
		myfile << cells_count[i][0] <<", "<<cells_count[i][1]<<", "<<cells_count[i][2]<<", "<<cells_count[i][3] <<", "<<cells_count[i][4]<< std::endl;
	}
	myfile.close();
	return 0; 
}