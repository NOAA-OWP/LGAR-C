/*
  author: Ahmad Jan (ahmad.jan@noaa.gov)
  date:  November 9, 2022
  - Includes unit test for bmi components and run model for a timestep to compute depth and soil moisture of wetting fronts
  - compares bmi input/output names, units, and memory allocation against benchmark
  - test passed criteria : test most of the bmi functions, i.e., bmi functions get/set expected values/behavior
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <iomanip> // std::setw
#include "../bmi/bmi.hxx"
#include "../include/bmi_lgar.hxx"

#define FAILURE 0
#define VERBOSITY 1

#define GREEN "\033[32m"
#define RED   "\033[31m"
#define BLUE  "\033[34m"
#define RESET "\033[0m"

#define SUCCESS 0

int main(int argc, char *argv[])
{
  BmiLGAR model, model_calib;

   if (argc != 2) {
    printf("Usage: ../build/lasam_unitest configs/unittest.txt \n");
    printf("Run the LASAM (Lumped Arid/Semi-arid Model through its BMI with a configuration file.\n");
    printf("Outputs are written to files `variables_data.csv and layers_data.csv`.\n");
    return SUCCESS;
  }

  std::cout<<"\n**************** BEGIN LASAM BMI UNIT TEST *******************\n";

  model.Initialize(argv[1]);
  model_calib.Initialize(argv[1]);
  
  // The following variables and names are benchmark values and names, any (unintended/inconsistent) change to the bmi or model will lead to test failure.
  std::cout<<"\n**************** TEST VALUES ************************************\n";
  int num_layers         = 3;       // total number of layers
  int num_wetting_fronts = 3;       // total number of wetting fronts
  bool test_status       = true;    // unit test status flag, if test fail the flag turns false
  int num_input_vars     = 3;       // total number of bmi input variables
  int num_output_vars    = 15;      // total number of bmi output variables

  // *************************************************************************************
  // names of the bmi input/output variables and the corresponding sizes, with units of input variables
  std::vector <std::string> var_name_input = {"precipitation_rate", "potential_evapotranspiration_rate",
					      "soil_temperature_profile"};

  std::vector <std::string> var_name_output = {"soil_moisture_wetting_fronts", "soil_depth_layers",
					       "soil_depth_wetting_fronts", "soil_num_wetting_fronts",
					       "precipitation", "potential_evapotranspiration",
					       "actual_evapotranspiration", "surface_runoff",
					       "giuh_runoff", "soil_storage", "total_discharge",
					       "infiltration", "percolation", "groundwater_to_stream_recharge",
					       "mass_balance"};

  int nbytes_input[] = {sizeof(double), sizeof(double), sizeof(double)};
  int nbytes_output[] = {int(num_wetting_fronts * sizeof(double)), int(num_layers * sizeof(double)),
			 int(num_wetting_fronts * sizeof(double)), sizeof(int), sizeof(double),
			 sizeof(double), sizeof(double), sizeof(double), sizeof(double), sizeof(double),
			 sizeof(double), sizeof(double), sizeof(double), sizeof(double), sizeof(double)};

  std::vector<std::string> bmi_units = {"mm h^-1", "mm h^-1", "K"};
  // *************************************************************************************

  // screen outout
  std::cout<<"Number layers:             "<< num_layers <<"\n";
  std::cout<<"Number of wetting fronts:  "<< num_layers <<"\n";
  std::cout<<"Number of input vars:      "<< num_input_vars <<"\n";
  std::cout<<"Number of output vars:     "<< num_output_vars <<"\n";

  std::cout<<"\nPulling information from BMI\n************************************\n";

  // Test get_component_name()
  std::string model_name = model.GetComponentName();
  if (VERBOSITY)
    std::cout<<"Model name: "<< model_name <<"\n";



  // *************************************************************************************
  // input and output variables checks
  int count_in = 0;
  int count_out = 0;
  std::vector<std::string> names_in;
  std::vector<std::string> names_out;

  // Test GetInputItemCount
  count_in = model.GetInputItemCount();

  if (VERBOSITY)
    std::cout<<"Input item count: "<< count_in <<"\n";

  if (count_in == num_input_vars)
    test_status &= true;
  else {
    test_status &= false;
    std::string passed = test_status == true ? "Yes" : "No";
    std::cout<<"Test passed: "<<passed<<"\n";
    std::stringstream errMsg;
    errMsg << "Number of input variables are different. "<< count_in << " != "<< num_input_vars << "\n";
    throw std::runtime_error(errMsg.str());
  }


  names_in = model.GetInputVarNames(); // call to BMI GetInputVarNames
  if (VERBOSITY) {
    std::cout<<"Input variable names \n";
    for (int i=0; i<count_in; i++)
      std::cout<<i<<" "<<names_in[i]<<"\n";
  }

  std::cout<<"**************************************** \n";
  // Test GetOutputItemCount
  count_out = model.GetOutputItemCount();
  if (VERBOSITY)
    std::cout<<"Output item count: "<< count_out<<"\n";

  // Test GetOutputVarNames
  names_out = model.GetOutputVarNames();
  if (VERBOSITY) {
    std::cout<<"Output variable names "<<names_out.size()<<"\n";
    for (int i=0; i<count_out; i++)
      std::cout<<i<<" "<<names_out[i]<<"\n";
  }
  if (count_out == num_output_vars)
    test_status &= true;
  else {
    test_status &= false;
    std::string passed = test_status == true ? "Yes" : "No";
    std::cout<<"Test passed: "<<passed<<"\n";
    std::stringstream errMsg;
    errMsg << "Number of output variables are different. "<< count_out <<" != "<< num_output_vars <<"\n";
    throw std::runtime_error(errMsg.str());
  }

  // *************************************************************************************

  // Test BMI: VARIABLE INFORMATION FUNCTIONS
  std::cout<<"\n**************** TEST BMI VARIABLE INFORMATION FUNCTIONS\n***************************\n";

  int grid, itemsize, nbytes;
  std::string location;
  std::string units;
  std::string vartype;

  // Loop over input variables and test important BMI functions
  for (int i=0; i<count_in; i++) {
    std::string var_name = names_in[i];
    if (VERBOSITY)
      std::cout<<"Input var_name: "<< var_name <<"\n";

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // get variable grid; should be non-negative.
    grid = model.GetVarGrid(var_name);
    if (VERBOSITY)
      std::cout<<"Grid: "<< grid <<"\n";

    if (grid >=0)
      test_status &= true;
    else {
      test_status &= false;
      std::string passed = test_status == true ? "Yes" : "No";
      std::cout<<"Test passed: "<<passed<<"\n";
      std::stringstream errMsg;
      errMsg << "grid < 0 \n";
      throw std::runtime_error(errMsg.str());
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    // get variable item size; should be sizeof(int) or sizeof(double)
    itemsize = model.GetVarItemsize(var_name);
    if (VERBOSITY)
      std::cout<<"Itemsize: "<< itemsize <<"\n";
    if (itemsize >0)
      test_status &= true;
    else {
      test_status &= false;
      std::string passed = test_status == true ? "Yes" : "No";
      std::cout<<"Test passed: "<< passed <<"\n";
      std::stringstream errMsg;
      errMsg << "itemsize < 0 \n";
      throw std::runtime_error(errMsg.str());
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    // Test get_var_location(); not sure if this is ever used
    location = model.GetVarLocation(var_name);
    if ( location == "") return FAILURE;
    if (VERBOSITY)
      std::cout<<" location: "<< location<<"\n";
    if (location == "")
      test_status &= false;
    else
      test_status &= true;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    // get input variable units and compare with the expected units (benchmark)
    units = model.GetVarUnits(var_name);
    if (VERBOSITY)
      std::cout<<" units: ["<< units <<"]\n";

    if (units != bmi_units[i]) {
      test_status &= false;
      std::string passed = test_status == true ? "Yes" : "No";
      std::cout<<"Test passed: "<<passed<<"\n";
      std::stringstream errMsg;
      errMsg << "units don't match "<< bmi_units[i] << " "<< units<<"\n";
      throw std::runtime_error(errMsg.str());
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    // Test get_var_type()
    vartype = model.GetVarType(var_name);
    if (VERBOSITY)
      std::cout<<" type: "<< vartype <<"\n";
    if (location == "")
      test_status &= false;
    else
      test_status &= true;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    // get_var_nbytes()
    nbytes = model.GetVarNbytes(var_name);
    if (nbytes == 0) return FAILURE;
    if (VERBOSITY)
      std::cout<<" nbytes: "<< nbytes <<"\n";

    if (var_name == var_name_input[i]) {
      if (nbytes == nbytes_input[i])
	test_status &= true;
      else {
	test_status &= false;
	std::string passed = test_status == true ? "Yes" : "No";
	std::cout<<"Test passed: "<<passed<<"\n";
	std::stringstream errMsg;
	errMsg << "Number of bytes for input var "<< var_name << " are "<< nbytes <<", but should be "<<nbytes_input[i]<<"\n";
	throw std::runtime_error(errMsg.str());
      }
    }
    else {
      std::stringstream errMsg;
      errMsg << "Input variable name"<< var_name<<" should be: soil_storage or soil_storage_change or soil_moisture_layered \n";
      throw std::runtime_error(errMsg.str());

    }
  }

  if (VERBOSITY)
    std::cout<<"\n*****************************************\n";

  for (int i=0; i<count_out; i++) {
    std::string var_name = names_out[i];
    if (VERBOSITY)
      std::cout<<"Output var_name: "<< var_name <<"\n";

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    // Test get_var_grid()
    grid = model.GetVarGrid(var_name);
    std::cout<<grid<<"\n";
    if (grid == -1) return -1;
    if (VERBOSITY)
      std::cout<<"Grid: "<< grid <<"\n";

    if (grid >=0)
      test_status &= true;
    else
      test_status &= false;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Test get_var_itemsize()
    itemsize = model.GetVarItemsize(var_name);
    if (itemsize == 0) return FAILURE;
    if (VERBOSITY)
      std::cout<<"Itemsize: "<< itemsize <<"\n";

    if (itemsize >0)
      test_status &= true;
    else
      test_status &= false;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    // Test get_var_location()
    location = model.GetVarLocation(var_name);
    if ( location == "") return FAILURE;
    if (VERBOSITY)
      std::cout<<" location:"<< location<<"\n";

    if (location == "")
      test_status &= false;
    else
      test_status &= true;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    // Test get_var_units()
    units = model.GetVarUnits(var_name);
    if (VERBOSITY)
      std::cout<<" units: ["<< units <<"]\n";

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    // Test get_var_type()
    vartype = model.GetVarType(var_name);
    if (vartype == "") return FAILURE;
    if (VERBOSITY)
      std::cout<<" type: "<< vartype <<"\n";

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    // get_var_nbytes()
    nbytes = model.GetVarNbytes(var_name);
    if (nbytes == 0) return FAILURE;

    if (VERBOSITY)
      std::cout<<" nbytes: "<< nbytes<<"\n";

    if (var_name == var_name_output[i]){
      if (nbytes == nbytes_output[i])
	test_status &= true;
      else {
	test_status &= false;
	std::string passed = test_status == true ? "Yes" : "No";
	std::cout<<"Test passed: "<<passed<<"\n";
	std::stringstream errMsg;
	errMsg << "Number of bytes for output var"<<var_name<< " should be "<<nbytes_output[i]<<"\n";
	throw std::runtime_error(errMsg.str());
      }
    }
    else {
      std::stringstream errMsg;
      errMsg << "Variable name "<< var_name<<" is not listed in the bmi names of the unit test. \n";
      throw std::runtime_error(errMsg.str());

    }

  }

  // Test BMI: MODEL GRID FUNCTIONS
  std::cout<<"\n \n**************** TEST BMI GRID FUNCTIONS***********************\n";
  int grid_id[] = {0,1,2,3};
  int grid_size_test[] = {1,1,num_layers,num_wetting_fronts};
  int grid_rank, grid_size;
  std::string grid_type;

  for (int i=0; i< 4; i++) {
    if (VERBOSITY)
      std::cout<<"Grid id "<< grid_id[i] <<"\n";

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_grid_rank()
    grid_rank = model.GetGridRank(grid_id[i]);
    if (grid_rank == FAILURE) return FAILURE;
    if (VERBOSITY)
      std::cout<<" rank: "<<grid_rank<<"\n";

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // Test get_grid_size
    grid_size = model.GetGridSize(grid_id[i]);
    if (grid_size == grid_size_test[i]) {
      test_status &= true;
      if (VERBOSITY)
	std::cout<<" grid size: "<<grid_size<<"\n";
    }
    else {
      test_status &= false;
      std::string passed = test_status == true ? "Yes" : "No";
      std::cout<<"Test passed: "<<passed<<"\n";
      std::stringstream errMsg;
      errMsg << "Grid size of should be "<<num_layers<<" or "<< num_wetting_fronts<<"\n";
      throw std::runtime_error(errMsg.str());
    }

  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/

  std::cout<<GREEN<<"\n";
  std::string passed = test_status > 0 ? "Yes" : "No";
  std::cout<<"\n| *************************************** \n";
  std::cout<<"| All tests passed until this point: "<<passed<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<RESET<<"\n";
  
  assert (test_status == true);

  // Test BMI: GET VALUE FUNCTIONS
  std::cout<<"\n\n************** TEST BMI GETTER SETTER FUNCTIONS********************************\n";

  std::cout<<"********** Input variables ***************** \n";
  // Loop through both input and output variables and call get/set_value_*()
  for (int i=0; i<count_in; i++) {
    std::string var_name = names_in[i];
    std::cout<<"Variable name: "<< var_name <<"\n";

    double *var = new double[1];
    double *dest = new double[1];
    int indices[] = {0};
    int len = 1;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    // Test get_value() at each timestep
    model.GetValue(var_name, &(var[0]));
    std::cout<<" Get value: "<< var[0] <<"\n";

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    // Test get_value_at_indices()
    model.GetValueAtIndices(var_name, dest, indices, len);
    std::cout<<" Get value at indices: " << dest[0]<<"\n";

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    // Test get_value_ptr()
    double *var_ptr = new double[1];
    var_ptr = (double*) model.GetValuePtr(var_name);
    std::cout<<" Get value ptr: "<<*var_ptr<<"\n";

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
    // Test BMI set_value_at_indices()
    double dest_new[] = {0.0};

    if (var_name == "precipitation_rate")
      dest_new[0] = 1.896; // in mm/hr
    else if (var_name == "potential_evapotranspiration_rate")
      dest_new[0] = 0.104; // in mm/hr

    double *dest_new_up = new double[1];

    model.SetValueAtIndices(var_name, &indices[0], len, &dest_new[0]);

    std::cout<<" Set value at indices: "<<dest_new[0]<<"\n";
    // get_value_at_indices to see if changed
    model.GetValueAtIndices(var_name, dest_new_up,  &indices[0], len);
    std::cout<<" Get value at indices: "<<dest_new_up[0]<<"\n";
    if (dest_new[0] == dest_new_up[0])
      test_status &= true;
    else
      test_status &= false;

  }

  passed = test_status > 0 ? "Yes" : "No";
  std::cout<<GREEN<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<"| All tests passed until this point: "<<passed<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<RESET<<"\n";

  std::cout<<"************* Output variables ***************** \n";
  model.Update();

  // Benchmark values of wetting fronts depth and moisture (b is for benchmark)
  //std::vector<double> depth_wf_b = {1.873813, 44.00,175.0, 200.0}; // in cm
  std::vector<double> depth_wf_b = {4.55355239489608365, 44.00,175.0, 200.0}; // in cm
  std::vector<double> theta_wf_b = {0.21371581122514613, 0.17270389607163267, 0.25211383152603861, 0.17959348005962811};

  int m_to_cm = 100;
  int m_to_mm = 1000;
  // note model outputs depths in meters

  int num_wf_base = 4; // number of wetting fronts after the rainfall
  // computed values (c is for computed) ; 4 = number of computed wetting fronts
  double *depth_wf_c = new double[num_wf_base];
  double *theta_wf_c = new double[num_wf_base];

  for (int i=0; i<count_out; i++) {

    std::string var_name = names_out[i];
    if (var_name == "soil_moisture_wetting_fronts") {
      std::cout<<"variable name: "<< var_name <<" "<<test_status<<"\n";

      model.GetValue(var_name, &theta_wf_c[0]);
      std::cout<<" Get value: "<< theta_wf_c[0] <<"\n";

      for (int k=0; k<num_wf_base; k++)
	if (fabs(theta_wf_b[k] - theta_wf_c[k]) < 0.0001)
	  test_status &= true;
	else
	  test_status &= false;

    }
    else if (var_name == "soil_depth_wetting_fronts") {
      std::cout<<"variable name: "<< var_name <<" "<<test_status<<"\n";

      model.GetValue(var_name, &depth_wf_c[0]);

      for (int k=0; k<num_wf_base; k++)
	if (fabs(depth_wf_b[k] - depth_wf_c[k]*m_to_cm) < 0.0001)
	  test_status &= true;
	else
	  test_status &= false;

    }
  }

  std::cout<<"\n"<<RED<<"Comparison: "<<RESET<<"Depths of the wetting fronts. \n";
  std::cout<<"Referance value | Simulated value | Difference \n";
  for (int i=0; i < num_wf_base; i++) {
    std::cout<< left << setw(18) << depth_wf_b[i]
	     << setw(18) << depth_wf_c[i] * m_to_cm
	     << setw(18) << abs(depth_wf_b[i] - depth_wf_c[i]*m_to_cm)<<"\n";
    assert (abs(depth_wf_b[i] - depth_wf_c[i]*m_to_cm) < 1.E-5);
  }

  std::cout<<"\n"<<RED<<"Comparison: "<<RESET<<"Moisture of the wetting fronts. \n";
  std::cout<<"Referance value | Simulated value | Difference \n";
  for (int i=0; i < num_wf_base; i++) {
    std::cout<< left << setw(18) << theta_wf_b[i]
	     << setw(18) << theta_wf_c[i]
	     << setw(1) << abs(theta_wf_b[i] - theta_wf_c[i])<<"\n";
    assert (abs(theta_wf_b[i] - theta_wf_c[i]) < 1.E-5);
  }

  passed = test_status > 0 ? "Yes" : "No";


  // check total infiltration, AET, and PET.
  double infiltration_check_mm = 1.896;  // in mm
  double AET_check_mm          = 0.02980092620558239; // in mm
  double PET_check_mm          = 0.104; // in mm
  double infiltration_computed = 0.0;
  double PET_computed          = 0.0;
  double AET_computed          = 0.0;

  model.GetValue("infiltration", &infiltration_computed);
  model.GetValue("potential_evapotranspiration", &PET_computed);
  model.GetValue("actual_evapotranspiration", &AET_computed);

  std::cout<<GREEN<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<"| All BMI Tests passed? "<< passed <<"\n";
  std::cout<<"| Infiltration [mm] : (benchmark vs computed) | "<< infiltration_check_mm <<" vs "
	   << infiltration_computed * m_to_mm <<"\n";
  std::cout<<"| PET [mm]          : (benchmark vs computed) | "<< PET_check_mm <<" vs "<< PET_computed * m_to_mm <<"\n";
  std::cout<<"| AET [mm]          : (benchmark vs computed) | "<< AET_check_mm <<" vs "<< AET_computed * m_to_mm <<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<RESET<<"\n";

  assert (test_status == true);
  
  // to print global mass balance
  //model.Finalize();

  if (fabs(infiltration_check_mm - infiltration_computed * m_to_mm) > 1.E-5) {
    std::stringstream errMsg;
    errMsg << "Error between benchmark and simulated infiltration is "<<
      fabs(infiltration_check_mm - infiltration_computed * m_to_cm) << " which is unexpected. \n";
    throw std::runtime_error(errMsg.str());
  }

  if (fabs(PET_check_mm - PET_computed * m_to_mm) > 1.E-5) {
    std::stringstream errMsg;
    errMsg << "Error between benchmark and simulated PET is "<< fabs(PET_check_mm - PET_computed * m_to_mm)
	   << " which is unexpected. \n";
    throw std::runtime_error(errMsg.str());
  }

  if (fabs(AET_check_mm - AET_computed * m_to_mm) > 1.E-5) {
    std::stringstream errMsg;
    errMsg << "Error between benchmark and simulated AET is "<< fabs(AET_check_mm - AET_computed * m_to_mm)
	   << " which is unexpected. \n";
     printf("benchmark AET: %.17lf \n", AET_check_mm);
     printf("computed AET: %.17lf \n", AET_computed*m_to_mm);

    throw std::runtime_error(errMsg.str());
  }


  std::cout<<GREEN<<"\n";
  std::cout<<"| *************************************** \n";
  std::cout<<"| LASAM Calibration test \n";
  
  double rain_precip = 1.2; // mm/hr
  double evapotran   = 3.0; // mm/hr
  

  // Testing Calibratable parameters
  double *smcmax   = new double[num_layers];
  double *vg_n     = new double[num_layers];
  double *vg_alpha = new double[num_layers];
  double *Ksat     = new double[num_layers];
  double field_capacity;
  double ponded_depth_max;

  double smcmax_set[]   = {0.3513, 0.3773, 0.3617};
  double vg_n_set[]     = {1.44260592334, 1.14965918354, 1.39051695249};
  double vg_alpha_set[] = {0.0021297, 0.0073272, 0.0027454};
  double Ksat_set[]     = {0.446, 0.0743, 0.415};
  double field_capacity_set = 103.3;
  double ponded_depth_max_set = 1.0;
 
  // Get the initial values set through the config file
  model_calib.GetValue("smcmax", &smcmax[0]);
  model_calib.GetValue("van_genuchten_n", &vg_n[0]);
  model_calib.GetValue("van_genuchten_alpha", &vg_alpha[0]);
  model_calib.GetValue("hydraulic_conductivity", &Ksat[0]);
  model_calib.GetValue("field_capacity", &field_capacity);
  model_calib.GetValue("ponded_depth_max", &ponded_depth_max);
  
  for (int i=0; i < num_layers; i++)
    std::cout<<"| Initial values: layer = "<< i+1 <<", smcmax = "<< smcmax[i]
	     <<", vg_n = "<< vg_n[i] <<", vg_alpha = " << vg_alpha[i]
	     <<", Ksat = "<< Ksat[i] <<"\n";
  printf("field_capacity: %lf \n", field_capacity);
  printf("ponded_depth_max: %lf \n", ponded_depth_max);

  // set the new values
  model_calib.SetValue("smcmax", &smcmax_set[0]);
  model_calib.SetValue("van_genuchten_n", &vg_n_set[0]);
  model_calib.SetValue("van_genuchten_alpha", &vg_alpha_set[0]);
  model_calib.SetValue("hydraulic_conductivity", &Ksat_set[0]);
  model_calib.SetValue("field_capacity", &field_capacity_set);
  model_calib.SetValue("ponded_depth_max", &ponded_depth_max_set);
 
  // get the new/updated values
  model_calib.GetValue("smcmax", &smcmax[0]);
  model_calib.GetValue("van_genuchten_n", &vg_n[0]);
  model_calib.GetValue("van_genuchten_alpha", &vg_alpha[0]);
  model_calib.GetValue("hydraulic_conductivity", &Ksat[0]);
  model_calib.GetValue("field_capacity", &field_capacity);
  model_calib.GetValue("ponded_depth_max", &ponded_depth_max);
 
  if (fabs(ponded_depth_max  - ponded_depth_max_set) > 1.E-5) {
    std::stringstream errMsg;
    errMsg << "Mismatch between ponded_depth_max calibrated values set and get "<< ponded_depth_max_set<<" "<< ponded_depth_max
      << " which is unexpected. \n";
    throw std::runtime_error(errMsg.str());
  }

  if (fabs(field_capacity  - field_capacity_set) > 1.E-5) {
    std::stringstream errMsg;
    errMsg << "Mismatch between field_capacity calibrated values set and get "<< field_capacity_set<<" "<< field_capacity
      << " which is unexpected. \n";
    throw std::runtime_error(errMsg.str());
  }

  for (int i=0; i < num_layers; i++) {
    
    if (fabs(smcmax[i]  - smcmax_set[i]) > 1.E-5) {
      std::stringstream errMsg;
      errMsg << "Mismatch between smcmax calibrated values set and get "<< smcmax_set[i]<<" "<< smcmax[i]
	     << " which is unexpected. \n";
      throw std::runtime_error(errMsg.str());
    }
    
    if (fabs(vg_n[i]  - vg_n_set[i]) > 1.E-5) {
      std::stringstream errMsg;
      errMsg << "Mismatch between vg_n calibrated values set and get "<< vg_n_set[i]<<" "<< vg_n[i]
	     << " which is unexpected. \n";
      throw std::runtime_error(errMsg.str());
    }
    
    if (fabs(vg_alpha[i]  - vg_alpha_set[i]) > 1.E-5) {
      std::stringstream errMsg;
      errMsg << "Mismatch between vg_alpha calibrated values set and get "<< vg_alpha_set[i]<<" "<< vg_alpha[i]
	     << " which is unexpected. \n";
      throw std::runtime_error(errMsg.str());
    }
    
    if (fabs(Ksat[i]  - Ksat_set[i]) > 1.E-5) {
      std::stringstream errMsg;
      errMsg << "Mismatch between hydraulic conductivity calibrated values set and get "<< Ksat_set[i]<<" "<< Ksat[i]
	     << " which is unexpected. \n";
      throw std::runtime_error(errMsg.str());
    }
    
  }

  std::cout<<"|  \n";
  for (int i=0; i < num_layers; i++)
    std::cout<<"| Calib. values: layer = "<< i+1 <<", smcmax = "<< smcmax[i]
	     <<", vg_n = "<< vg_n[i] <<", vg_alpha = " << vg_alpha[i]
	     <<", Ksat = "<< Ksat[i] <<"\n";
  printf("field_capacity = %lf \n", field_capacity);
  printf("ponded_depth_max = %lf \n", ponded_depth_max);
  std::cout<<"| *************************************** \n";
  std::cout<<"| LASAM Calibration test passed? YES \n";
  std::cout<<RESET<<"\n";
  // set forcing data for the timestep
  model_calib.SetValue("precipitation_rate", &rain_precip);
  model_calib.SetValue("potential_evapotranspiration_rate", &evapotran);
  model_calib.Update();
  
  //model_calib.Finalize();
  return FAILURE;
}
