/*
  author : Ahmad Jan
  year   : 2022
  email  : ahmad.jan@noaa.gov
  Description: the code runs Lumped Aric/semi-arid model (LASAM) through a BMI
  Input: Precipitation, Potential ET, soil types, initial conditions
  Output: see the list below (also the code writes the variables and state of the wetting fronts
          to two separate files for analysis and visualization
*/


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "fstream"
#include <iomanip>

#include "../bmi/bmi.hxx"
#include "../include/all.hxx"
#include "../include/bmi_lgar.hxx"

// module finds forcing file name in the config file
std::string GetForcingFile(std::string config_file);

// module read forcings (precipitation and PET)
void ReadForcingData(std::string config_file, std::vector<std::string>& time, std::vector<double>& precip, std::vector<double>& pet);


#define SUCCESS 0

int main(int argc, char *argv[])
{

  BmiLGAR model_state;

  bool is_IO_supress = false; // if true no output files will be written

  if (argc != 2) {
    printf("Usage: ./build/xlgar CONFIGURATION_FILE \n");
    printf("Run the LASAM (Lumped Arid/semi-aric Model through its BMI with a configuration file.\n");
    printf("Outputs are written to files `variables_data.csv and layers_data.csv`.\n");
    return SUCCESS;
  }

  clock_t start_time, end_time;
  double elapsed;
  start_time = clock();

  model_state.Initialize(argv[1]);


  std::string var_name_precip = "precipitation_rate";
  std::string var_name_pet    = "potential_evapotranspiration_rate";
  std::string var_name_wf     = "soil_moisture_wetting_fronts";
  std::string var_name_thickness_wf = "soil_depth_wetting_fronts";

  int num_output_var = 11;
  std::vector<std::string> output_var_names(num_output_var);
  std::vector<double> output_var_data(num_output_var);

  output_var_names[0]  = "precipitation";
  output_var_names[1]  = "potential_evapotranspiration";
  output_var_names[2]  = "actual_evapotranspiration";
  output_var_names[3]  = "surface_runoff"; // direct surface runoff
  output_var_names[4]  = "giuh_runoff";
  output_var_names[5]  = "soil_storage";
  output_var_names[6]  = "total_discharge";
  output_var_names[7]  = "infiltration";
  output_var_names[8]  = "percolation";
  output_var_names[9]  = "conceptual_reservoir_to_stream_discharge";
  output_var_names[10] = "mass_balance";


  // total number of timesteps

  // get time steps
  double endtime = model_state.GetEndTime();
  double current_time = model_state.GetCurrentTime();
  double timestep = model_state.GetTimeStep();
  int nsteps = int((endtime-current_time)/timestep); // time steps remaining in this invocation

  std::vector<std::string> time;
  std::vector<double> precipitation;
  std::vector<double> PET;

  ReadForcingData(argv[1], time, precipitation, PET);

  assert (nsteps <= int(PET.size()) ); // assertion to ensure that nsteps are less or equal than the input data
  
  if (verbosity.compare("high") == 0 && !is_IO_supress) {
    std::cout<<"Variables are written to file           : \'data_variables.csv\' \n";
    std::cout<<"Wetting fronts state is written to file : \'data_layers.csv\' \n";
  }

  FILE *outdata_fptr = NULL;
  FILE *outlayer_fptr = NULL;
  FILE *outnonvadose_fptr = NULL;
  FILE *outgiuh_fptr = NULL;

  if (!is_IO_supress) {
    outdata_fptr = fopen("data_variables.csv", "w");  // write output variables (e.g. infiltration, storage etc.) to this file pointer
    outlayer_fptr = fopen("data_layers.csv", "w");    // write output layers to this file pointer
    outnonvadose_fptr = fopen("data_non_vadose_state.csv", "w");
    outgiuh_fptr = fopen("data_giuh_state.csv", "w");

    // write heading (variable names)
    fprintf(outdata_fptr,"Time,");
    fprintf(outnonvadose_fptr,
        "Time,CR_fast_storage_cm,CR_slow_storage_cm,volon_timestep_cm,"
        "runoff_in_prev_step,precip_previous_timestep_cm,cache_fluxes,"
        "cache_count,previous_AET,previous_PET,previous_recharge,"
        "accumulated_PET,accumulated_free_drainage,time_s,timesteps\n");
    for (int j = 0; j < num_output_var; j++) {
      fprintf(outdata_fptr,"%s",output_var_names[j].c_str());
      if (j == num_output_var-1)
	fprintf(outdata_fptr,"\n");
      else
      fprintf(outdata_fptr,",");
    }

  }

  // model timestep and forcing timestep are read from a config file in lgar.cxx
  //  double dt = 3600;
  for (int i = 0; i < nsteps; i++) {

    if (verbosity.compare("none") != 0) {
      std::cout<<"===============================================================\n";
      std::cout<<"Real time | "<<time[i]<<"\n";
      std::cout<<"Rainfall [mm/h], PET [mm/h] = "<<precipitation[i]<<" , "<<PET[i]<<"\n";
    }

    model_state.SetValue(var_name_precip, &precipitation[i]);
    model_state.SetValue(var_name_pet, &PET[i]);

    //model_state.UpdateUntil(dt); // Update model

    model_state.Update(); // Update model

    if (!is_IO_supress) {
      int num_wetting_fronts =  model_state.get_model()->lgar_bmi_params.num_wetting_fronts;

      double *soil_moisture_wetting_front = new double[num_wetting_fronts];
      double *soil_thickness_wetting_front = new double[num_wetting_fronts];

      model_state.GetValue(var_name_wf,&soil_moisture_wetting_front[0]);
      model_state.GetValue(var_name_thickness_wf,&soil_thickness_wetting_front[0]);

      // write bmi output variables to file
      fprintf(outdata_fptr,"%s,",time[i].c_str());

      for (int j = 0; j < num_output_var; j++) {
	std::string name = output_var_names[j];
	double value = 0.0;
	model_state.GetValue(name,&value);
	fprintf(outdata_fptr,"%6.15f",value);
	if (j == num_output_var-1)
	  fprintf(outdata_fptr,"\n");
	else
	  fprintf(outdata_fptr,",");
      }

      fprintf(outnonvadose_fptr, "%s,", time[i].c_str());
      write_non_vadose_state(outnonvadose_fptr, model_state.get_model());

      fprintf(outgiuh_fptr, "%s,", time[i].c_str());
      write_giuh_runoff_queue_state(
          outgiuh_fptr,
          model_state.get_giuh_runoff_queue(),
          model_state.get_num_giuh_ordinates());

      // write layers data to file
      fprintf(outlayer_fptr,"# Timestep = %d, %s \n", i, time[i].c_str());
      write_state(outlayer_fptr, model_state.get_model()->head);
      delete [] soil_moisture_wetting_front;
      delete [] soil_thickness_wetting_front;
    }

  }

  // do final mass balance ( inside Finalize() ) and finish the simulation
  model_state.Finalize();
  if (outdata_fptr) {
    fclose(outdata_fptr);
    fclose(outlayer_fptr);
    fclose(outnonvadose_fptr);
    fclose(outgiuh_fptr);
  }

  end_time = clock();

  elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC;

  std::cout<<setprecision(4);
  std::cout<<"Time                        =   "<< elapsed <<" sec \n";

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


extern void write_state(FILE *out, struct wetting_front* head){

  struct wetting_front *current = head;

  fprintf(out, "[");
  while(current != NULL)
  {
    if (current == head)
      fprintf(out,"(%.17g,%.17g,%d,%d,%d,%.17g,%.17g)",
              current->depth_cm, current->theta, current->layer_num,
              current->front_num, current->to_bottom ? 1 : 0,
              current->psi_cm, current->dzdt_cm_per_h);
    else
      fprintf(out,"|(%.17g,%.17g,%d,%d,%d,%.17g,%.17g)",
              current->depth_cm, current->theta, current->layer_num,
              current->front_num, current->to_bottom ? 1 : 0,
              current->psi_cm, current->dzdt_cm_per_h);
  current = current->next;
  }
  fprintf(out, "]\n");

}

extern void write_non_vadose_state(FILE *out, struct model_state* state)
{
  fprintf(out,
          "CR_fast_storage_cm=%.17g,"
          "CR_slow_storage_cm=%.17g,"
          "volon_timestep_cm=%.17g,"
          "runoff_in_prev_step=%d,"
          "precip_previous_timestep_cm=%.17g,"
          "cache_fluxes=%d,"
          "cache_count=%d,"
          "previous_AET=%.17g,"
          "previous_PET=%.17g,"
          "previous_recharge=%.17g,"
          "accumulated_PET=%.17g,"
          "accumulated_free_drainage=%.17g,"
          "time_s=%.17g,"
          "timesteps=%d\n",
          state->lgar_mass_balance.CR_fast_storage_cm,
          state->lgar_mass_balance.CR_slow_storage_cm,
          state->lgar_mass_balance.volon_timestep_cm,
          state->lgar_bmi_params.runoff_in_prev_step ? 1 : 0,
          state->lgar_bmi_params.precip_previous_timestep_cm,
          state->lgar_mass_balance.cache_fluxes ? 1 : 0,
          state->lgar_bmi_params.cache_count,
          state->lgar_mass_balance.previous_AET,
          state->lgar_mass_balance.previous_PET,
          state->lgar_mass_balance.previous_recharge,
          state->lgar_mass_balance.accumulated_PET,
          state->lgar_mass_balance.accumulated_free_drainage,
          state->lgar_bmi_params.time_s,
          state->lgar_bmi_params.timesteps);
}

extern void write_giuh_runoff_queue_state(
    FILE *out,
    const double *giuh_runoff_queue,
    int num_giuh_ordinates)
{
  fprintf(out, "num_giuh_ordinates=%d,queue=[", num_giuh_ordinates);

  for (int i = 0; i <= num_giuh_ordinates; i++) {
    if (i > 0) fprintf(out, ",");
    fprintf(out, "%.17g", giuh_runoff_queue[i]);
  }

  fprintf(out, "]\n");
}
