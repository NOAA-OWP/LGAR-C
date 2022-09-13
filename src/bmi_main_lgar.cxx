#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>

#include "../bmi/bmi.hxx"
#include "../include/all.h"
#include "../include/bmi_lgar.hxx"
//#include "../include/lgar.h"
//#include "../include/all.h"

struct wetting_front *head = NULL;  //<- this pointer stores the address in memory of the first member of the linked
                                    //   list containing all the wetting fronts.  That address is called "head".  
                                    //   The contents of struct wetting_front are defined in "all.h"

struct wetting_front *head_previous = NULL;

#define SUCCESS 0
int main(int argc, char *argv[])
{
  /*
  BmiSoilMoistureProfile model;
  
  if (argc != 2) {
    printf("Usage: run_bmifrozensoilcxx CONFIGURATION_FILE\n\n");
    printf("Run the frozensoilcxx model through its BMI with a configuration file.\n");
    printf("Output is written to the file `bmifrozensoilcxx.out`.\n");
    return SUCCESS;
  }

  FILE *fp = fopen("bmi_file.out", "w");
  fprintf(fp, "Configuration file = %s\n", argv[1]);
  fprintf(fp, "Initializing... ");
  
  model.Initialize(argv[1]);
  
  fprintf(fp, "done\n");

  {
    std::string model_name;
    model_name = model.GetComponentName();
    fprintf(fp, "%s\n", model_name.c_str());
  }

  {
    std::string var_name_s = "soil_storage";
    std::string var_name_sc = "soil_storage_change";
    std::string var_name_wt = "soil_water_table";
    std::string var_name_smc = "soil_moisture_profile";
    std::string var_name_smcl = "soil_moisture_layered";
    std::string var_name_smc_bmi = "soil_storage_model";
    
    int grid, rank, *shape;
    double *var_s = NULL;
    double *var_sc = NULL;
    double *water_table_thickness_bmi = NULL;

    fprintf(fp, "variable = %s\n", var_name_s.c_str());
    fprintf(fp, "variable = %s\n", var_name_sc.c_str());
    fprintf(fp, "variable = %s\n", var_name_wt.c_str());
    fprintf(fp, "variable = %s\n", var_name_smc.c_str());
    fprintf(fp, "variable = %s\n", var_name_smcl.c_str());
    fprintf(fp, "variable = %s\n", var_name_smc_bmi.c_str());
    
    grid = model.GetVarGrid(var_name_smc);

    rank = model.GetGridRank(grid);
    fprintf(fp, "rank = %d\n", rank);
    shape = new int[rank];
    model.GetGridShape(grid, shape);

    fprintf(fp, "shape = %d x %d x %d\n", shape[0],1,1);
  */
    /****************************************************************************/
    // unit test data for conceptual soil reservoir
    //double SMCT[] ={0.32207, 0.333438, 0.367336, 0.439}; // soil_moisture_profile
    double SMCT[] = {0.3375634,0.34234789,0.34752731,0.35330897,0.35974973,0.3669139,0.37517657,0.38464055,0.39596964,0.40977129,0.42703622,0.439,0.439,0.439,0.439,0.439,0.439,0.439,0.439,0.439};
    double water_table_thickness = 0.490438; // in meters
    enum option { Conceptual = 1, Layered = 2};
    /****************************************************************************/
    // Set values
    double storage_m = 0.8; //0.526328;
    double storage_change_m = -0.000472;
    double *storage_m_ptr = &storage_m;
    double *storage_change_m_ptr = &storage_change_m;
    double smc_layers[] = {0.25, 0.15, 0.1, 0.12};
    /*
    int soil_moisture_profile_option;

    model.GetValue(var_name_smc_bmi,&soil_moisture_profile_option);

    model.SetValue(var_name_s,storage_m_ptr);

    model.SetValue(var_name_sc,storage_change_m_ptr);

    model.SetValue(var_name_smcl,&smc_layers[0]);
    
    var_s = (double *)model.GetValuePtr(var_name_s);
    var_sc = (double *)model.GetValuePtr(var_name_sc);

    std::cout<<"soil_storage_model: "<<soil_moisture_profile_option<<"\n";
    std::cout<<"storage: "<<*var_s<<"\n";
    std::cout<<"storage change: "<<*var_sc<<"\n";
    
    model.Update();

    water_table_thickness_bmi = (double *)model.GetValuePtr(var_name_wt);
    
    if (soil_moisture_profile_option == Conceptual) {
      std::cout<<"Check: water table thickness = "<<water_table_thickness<<" | water table thickness bmi = "<<*water_table_thickness_bmi<<"\n";
      double diff = std::fabs(water_table_thickness - *water_table_thickness_bmi);
      assert (diff < 1e-4);
    }
    
    // Get values
    double *var_smc = new double[shape[0]];
    
    model.GetValue(var_name_smc,&var_smc[0]);
      
    if (soil_moisture_profile_option == Conceptual) {
      std::cout<<"Referance value | Simulated value | Difference \n";
      for (int i=0; i < shape[0]; i++) {
	std::cout<< left << setw(18) << var_smc[i]
		 << setw(18) << SMCT[i]
		 << setw(1) << abs(var_smc[i] - SMCT[i])<<"\n";
	assert (abs(var_smc[i] - SMCT[i]) < 1.E-6);     
	fprintf(fp, "%6.4e", var_smc[i]);
	fprintf(fp, "\n");
      }
    }
  }
  
  fprintf(fp, "Finalizing... ");

  model.Finalize();
  fprintf(fp, "done\n");
  fclose(fp);
    */
  return SUCCESS;
}
