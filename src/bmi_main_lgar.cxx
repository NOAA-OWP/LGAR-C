#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "fstream"

#include "../bmi/bmi.hxx"
#include "../include/all.hxx"
#include "../include/bmi_lgar.hxx"


std::string GetForcingFile(std::string config_file);

void ReadForcingData(std::string config_file, std::vector<std::string>& time, std::vector<double>& precip, std::vector<double>& pet);

struct wetting_front *head = NULL;  //<- this pointer stores the address in memory of the first member of the linked
                                    //   list containing all the wetting fronts.  That address is called "head".  
                                    //   The contents of struct wetting_front are defined in "all.h"

struct wetting_front *state_previous = NULL;

#define SUCCESS 0
int main(int argc, char *argv[])
{
  
  BmiLGAR lgar_model;
  
  if (argc != 2) {
    printf("Usage: run_bmifrozensoilcxx CONFIGURATION_FILE\n\n");
    printf("Run the frozensoilcxx model through its BMI with a configuration file.\n");
    printf("Output is written to the file `bmifrozensoilcxx.out`.\n");
    return SUCCESS;
  }
  
  FILE *fp = fopen("bmi_file.out", "w");
  fprintf(fp, "Configuration file = %s\n", argv[1]);
  fprintf(fp, "Initializing... ");

  lgar_model.Initialize(argv[1]);
  
  fprintf(fp, "done\n");
  
  {
    std::string model_name;
    model_name = lgar_model.GetComponentName();
    fprintf(fp, "%s\n", model_name.c_str());
  }
  
  {
    std::string var_name_precip = "precipitation";
    std::string var_name_pet = "potential_evapotranspiration";
    
    int grid, rank, *shape;
    double *var_s = NULL;
    double *var_sc = NULL;

    fprintf(fp, "variable = %s\n", var_name_precip.c_str());
    fprintf(fp, "variable = %s\n", var_name_pet.c_str());
    
    grid = lgar_model.GetVarGrid(var_name_precip);

    rank = lgar_model.GetGridRank(grid);
    fprintf(fp, "rank = %d\n", rank);
    shape = new int[rank];
    lgar_model.GetGridShape(grid, shape);

    fprintf(fp, "shape = %d x %d x %d\n", shape[0],1,1);
  }


  std::string filename = GetForcingFile(argv[1]); // forcing file name
  std::cout<<"forcings "<<filename<<"\n";

  std::ifstream forcing_p;
  forcing_p.open(filename);

  std::string line, cell;
  getline(forcing_p,line); // read heading and skip
  //std::stringstream lineStream(line);
  
  int nsteps = 1;//57;
  std::vector<std::string> time;
  std::vector<double> precipitation;
  std::vector<double> PET;
  
  ReadForcingData(argv[1], time, precipitation, PET);

  //for (int i =0; i < 7; i++)
  //  std::cout<<"Time, P, PET = "<<time[i]<<" "<<precipitation[i]<<" "<<PET[i]<<"\n";
  
  for (int i = 0; i < nsteps; i++) {
    std::cout<<"----------------------------------- \n";
    std::cout<<"Timestep | "<<i<<" "<<time[i]<<"\n";
    std::cout<<"P, PET = "<<precipitation[i]<<" "<<PET[i]<<"\n";
    lgar_model.SetValue("precipitation", &precipitation[i]);
    lgar_model.SetValue("potential_evapotranspiration", &PET[i]);

    lgar_model.Update(); // Update model
    /*
    

    ftm_bmi_model.Update(); // Update model

    ftm_bmi_model.GetValue("ice_fraction_schaake",&ice_frac_v);

    ice_fraction[i] = ice_frac_v;
    
    if (golden_test)
      outfile << i+1 << "," <<ice_frac_v << "\n";
    */ 
  }


  lgar_model.global_mass_balance();
 
 
 /* 
  fprintf(fp, "Finalizing... ");

  model.Finalize();
  fprintf(fp, "done\n");
  fclose(fp);
    */
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


std::string
GetForcingFile(std::string config_file)
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
  fp.close();
  
  return forcing_file;
 
}
