#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "fstream"
#include <iomanip>

#include "../bmi/bmi.hxx"
#include "../include/all.hxx"
#include "../include/bmi_lgar.hxx"


std::string GetForcingFile(std::string config_file);

void ReadForcingData(std::string config_file, std::vector<std::string>& time, std::vector<double>& precip, std::vector<double>& pet);


#define SUCCESS 0
int main(int argc, char *argv[])
{
  
  BmiLGAR lgar_model;

  bool is_IO_supress = false;

  if (argc != 2) {
    printf("Usage: ./build/xlgar CONFIGURATION_FILE \n");
    printf("Run the LASAM (Lumped Arid/semi-aric Model through its BMI with a configuration file.\n");
    printf("Outputs are written to files `variables_data.csv and layers_data.csv`.\n");
    return SUCCESS;
  }

  time_t result = time(NULL);
  clock_t start_time, end_time;
  double elapsed;
  start_time = clock();
 
  lgar_model.Initialize(argv[1]);
  
  
  std::string var_name_precip = "precipitation_rate";
  std::string var_name_pet = "potential_evapotranspiration_rate";
  std::string var_name_wf = "soil_moisture_wetting_fronts";
  std::string var_name_thickness_wf = "soil_thickness_wetting_fronts";

  int num_output_var = 9;
  std::vector<std::string> output_var_names(num_output_var);
  std::vector<double> output_var_data(num_output_var);
  
  output_var_names[0] = "precipitation";
  output_var_names[1] = "potential_evapotranspiration";
  output_var_names[2] = "actual_evapotranspiration";
  output_var_names[3] = "surface_runoff"; // direct surface runoff
  output_var_names[4] = "giuh_runoff";
  output_var_names[5] = "soil_storage";
  output_var_names[6] = "total_discharge";
  output_var_names[7] = "infiltration";
  output_var_names[8] = "percolation";

  
  int nsteps = 7500; //2700; //7500.0;//57;
  std::vector<std::string> time;
  std::vector<double> precipitation;
  std::vector<double> PET;
  
  ReadForcingData(argv[1], time, precipitation, PET);

  if (verbosity.compare("high") == 0 && !is_IO_supress) {
    std::cout<<"Variables are written to file           : \'data_variables.csv\' \n";
    std::cout<<"Wetting fronts state is written to file : \'data_layers.csv\' \n"; 
  }

  FILE *outdata_fptr = NULL;
  FILE *outlayer_fptr = NULL;
  
  if (!is_IO_supress) {
    outdata_fptr = fopen("data_variables.csv", "w");      // write output variables (e.g. infiltration, storage etc.) to this file pointer
    outlayer_fptr = fopen("data_layers.csv", "w");   // write output layers to this file pointer

    // write heading (variable names)
    fprintf(outdata_fptr,"Time,");
    for (int j = 0; j < num_output_var; j++) {
      fprintf(outdata_fptr,"%s",output_var_names[j].c_str());
      if (j == num_output_var-1)
	fprintf(outdata_fptr,"\n");
      else
      fprintf(outdata_fptr,",");
    }

  }
  
  for (int i = 0; i < nsteps; i++) {
    
    if (verbosity.compare("none") != 0) {
      std::cout<<"---------------------------------------------------------\n";
      std::cout<<"Timestep | "<<i<<" , "<<time[i]<<"\n";
      std::cout<<"Rainfall [mm/h], PET [mm/h] = "<<precipitation[i]<<" , "<<PET[i]<<"\n";
    }
    
    lgar_model.SetValue(var_name_precip, &precipitation[i]);
    lgar_model.SetValue(var_name_pet, &PET[i]);

    lgar_model.Update(); // Update model

    if (!is_IO_supress) {
      int num_wetting_fronts =  lgar_model.get_model()->lgar_bmi_params.num_wetting_fronts;
      
      double *soil_moisture_wetting_front = new double[num_wetting_fronts];
      double *soil_thickness_wetting_front = new double[num_wetting_fronts];
      
      lgar_model.GetValue(var_name_wf,&soil_moisture_wetting_front[0]);
      lgar_model.GetValue(var_name_thickness_wf,&soil_thickness_wetting_front[0]);
      
      // write bmi output variables to file
      fprintf(outdata_fptr,"%s,",time[i].c_str());
      
      for (int j = 0; j < num_output_var; j++) {
	std::string name = output_var_names[j];
	double value = 0.0;
	lgar_model.GetValue(name,&value);
	fprintf(outdata_fptr,"%6.10f",value);
	if (j == num_output_var-1)
	  fprintf(outdata_fptr,"\n");
	else
	  fprintf(outdata_fptr,",");
      }
      
      
      // write layers data to file
      fprintf(outlayer_fptr,"# Timestep = %d, %s \n", i, time[i].c_str());
      write_state(outlayer_fptr);
    }
    
  }


  lgar_model.global_mass_balance();

  if (outdata_fptr) {
    fclose(outdata_fptr);
    fclose(outlayer_fptr);
  }

  end_time = clock();
  
  elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC;

  std::cout<<setprecision(4);
  std::cout<<"Time                       =   "<< elapsed <<" sec \n";
   
  return SUCCESS;
}



void
ReadForcingData(std::string config_file, std::vector<std::string>& time, std::vector<double>& precip, std::vector<double>& pet)
{
  // get the forcing file from the config file

  std::ifstream file;
  file.open(config_file);

  if (!file) {
    std::stringstream errMsg;
    errMsg << config_file << " does not exist";
    throw std::runtime_error(errMsg.str());
  }

  std::string forcing_file;
  bool is_forcing_file_set=false;
  
  while (file) {
    std::string line;
    std::string param_key, param_value;

    std::getline(file, line);

    int loc_eq = line.find("=") + 1;
    param_key = line.substr(0, line.find("="));
    param_value = line.substr(loc_eq,line.length());

    if (param_key == "forcing_file") {
      forcing_file = param_value;
      is_forcing_file_set = true;
      break;
    }
  }

  if (!is_forcing_file_set) {
    std::stringstream errMsg;
    errMsg << config_file << " does not provide forcing_file";
    throw std::runtime_error(errMsg.str());
  }
  
  std::ifstream fp;
  fp.open(forcing_file);
  if (!fp) {
    cout<<"file "<<forcing_file<<" doesn't exist. \n";
    abort();
  }
  
  std::string line, cell;
  
  //read first line of strings which contains forcing variables names.
  std::getline(fp, line);

  while (fp) {
    std::getline(fp, line);
    std::stringstream lineStream(line);
    int count = 0;
    while(std::getline(lineStream,cell, ',')) {
      
      if (count == 0) {
	time.push_back(cell);
	count++;
	continue;
      }
      else if (count == 1) {
	precip.push_back(stod(cell));
	count++;
	continue;
      }
      else if (count == 2) {
	pet.push_back(stod(cell));
	count +=1;
	continue;
      }

    }

  }

 
}


extern void write_state(FILE *out){

  struct wetting_front *current = head;

  fprintf(out, "[");
  while(current != NULL)
  {
    if (current == head)
      fprintf(out,"(%lf,%lf,%d,%d,%lf)",current->depth_cm*10., current->theta, current->layer_num,current->front_num, current->psi_cm*10.);
    else
      fprintf(out,"|(%lf,%lf,%d,%d,%lf)",current->depth_cm*10., current->theta, current->layer_num,current->front_num, current->psi_cm*10.);
  current = current->next;
  }
  fprintf(out, "]\n");

}
