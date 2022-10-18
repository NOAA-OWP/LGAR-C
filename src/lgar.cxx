#ifndef LGAR_CXX_INCLUDED
#define LGAR_CXX_INCLUDED

#include "../include/all.hxx"
//#include "../src/aet.cxx"
#include <iostream>
#include <fstream>
//#include <stdio.h>
//#include <stdlib.h>
#include <string.h>
#include <sstream>

using namespace std;

#define VERBOSE 0
/* verbose = 0 (none)
   verbose = 1 (low)
   verbose = 2 (medium)
   verbose = 3 (high)
*/

/*##################################################*/
/*##################################################*/
/*##################################################*/
//--------------------------------------------------
//
//     LL          GGG
//     LL        GGG GGG       AAA  A    RR   RRRR
//     LL       GG    GGG    AAA AAAA     RR RR  RR
//     LL       GG     GG   AA     AA     RRR     
//     LL      GGG    GGG  AAA     AA     RR
//     LL       GG  GG GG   AA     AAA    RR
//     LL        GGGG  GG    AAA  AA A    RR
//     LL              GG      AAAA   AA  RR
//     LL              GG 
//     LLLLLLLL  GG   GG
//                 GGGG  
//
// The lgar specific functions based on soil physics
//-----------------------------------------------------
//
//            SKETCH SHOWING 3 SOIL LAYERS AND 4 WETTING FRONTS
//
//    theta_r                                  theta1             theta_e
//      --------------------------------------------------------------   --depth_cm = 0    -------> theta
//      r                                        f                   e
//      r                                        f                   e
//      r                                --------1                   e   --depth_cm(f1)
//      r                               f         \                  e
//      r     1st soil layer            f          \                 e
//      r                               f wetting front number       e
//      r                               f /                          e
//      r                               f/                           e
//      --------------------------------2-----------------------------  -- depth_cm(f2)
//         r                                      f              e
//         r     2nd soil layer                   f              e
//         r                                      f              e
//  |      ---------------------------------------3---------------      -- depth_cm(f3)
//  |        r                                          f     e
//  |        r                                          f     e
//  |        r       3rd  soil layer                    f     e
//  |        r                                          f     e
//  V        -------------------------------------------4------         -- depth_cm(f4)
// depth


extern void lgar_initialize(string config_file, struct lgar_model_ *model)
{
  
  InitFromConfigFile(config_file, model);
  model->lgar_bmi_params.shape[0] = model->lgar_bmi_params.num_layers;

  // initial number of wetting fronts are same are number of layers
  model->lgar_bmi_params.num_wetting_fronts = model->lgar_bmi_params.num_layers;
  model->lgar_bmi_params.soil_thickness_wetting_fronts = new double[model->lgar_bmi_params.num_wetting_fronts];
  model->lgar_bmi_params.soil_moisture_wetting_fronts = new double[model->lgar_bmi_params.num_wetting_fronts];

  // initialize thickness/depth and soil moisture of wetting fronts (used for model coupling)
  struct wetting_front *current = head;
  for (int i=0; i<model->lgar_bmi_params.num_wetting_fronts; i++) {
    assert (current != NULL);
    model->lgar_bmi_params.soil_moisture_wetting_fronts[i] = current->theta;
    model->lgar_bmi_params.soil_thickness_wetting_fronts[i] = current->depth_cm;
    current = current->next;
  }

}


extern void InitFromConfigFile(string config_file, struct lgar_model_ *model)
{
  
  ifstream fp; //FILE *fp = fopen(config_file.c_str(),"r");
  fp.open(config_file);

  // loop over the variables in the file to see if verbosity is provided, if not default is "none" (prints nothing)
  while (fp) {
    string line;
    string param_key, param_value, param_unit;
    
    getline(fp, line);
    
    int loc_eq = line.find("=") + 1;
    int loc_u = line.find("[");
    param_key = line.substr(0,line.find("="));

    param_value = line.substr(loc_eq,loc_u - loc_eq);
    
    if (param_key == "verbosity") {
      verbosity = param_value;
      if (verbosity.compare("none") != 0) {
	std::cout<<"Verbosity is set to \' "<<verbosity<<"\' \n";
	std::cout<<"          *****         \n";
      }

      fp.clear();
      break;
    }
  }

  
  if (verbosity.compare("none") != 0) {
    std::cout<<"------------- Initialization from config file ---------------------- \n";
  }

  
  bool is_layer_thickness_set = false;
  bool is_initial_psi_set = false;
  bool is_timestep_set = false;
  bool is_forcing_resolution_set = false;
  bool is_layer_soil_type_set = false;
  bool is_wilting_point_psi_cm_set = false;
  bool is_soil_params_file_set = false;
  bool is_max_soil_types_set = false;
  bool is_giuh_ordinates_set = false;
  bool is_verbosity_set = false;
  
  string soil_params_file;


  
  while (fp) {

    string line;
    string param_key, param_value, param_unit;
    
    getline(fp, line);
   
    int loc_eq = line.find("=") + 1;
    int loc_u = line.find("[");
    param_key = line.substr(0,line.find("="));

    bool is_unit = line.find("[") != string::npos;

    if (is_unit)
      param_unit = line.substr(loc_u,line.find("]")+1);
    else
      param_unit = "";

    param_value = line.substr(loc_eq,loc_u - loc_eq);
    
    if (param_key == "layer_thickness") {
      vector<double> vec = ReadVectorData(param_value);
      
      model->lgar_bmi_params.layer_thickness_cm = new double[vec.size()+1];
      model->lgar_bmi_params.cum_layer_thickness_cm = new double[vec.size()+1];

      model->lgar_bmi_params.layer_thickness_cm[0] = 0.0; // the value at index 0 is never used
      // calculate the cumulative (absolute) depth from land surface to bottom of each soil layer
      model->lgar_bmi_params.cum_layer_thickness_cm[0] = 0.0;
      
      for (unsigned int layer=1; layer <= vec.size(); layer++) {
      	model->lgar_bmi_params.layer_thickness_cm[layer] = vec[layer-1];
	model->lgar_bmi_params.cum_layer_thickness_cm[layer] = model->lgar_bmi_params.cum_layer_thickness_cm[layer-1] + vec[layer-1];
      }
      
      model->lgar_bmi_params.num_layers = vec.size();
 
      model->lgar_bmi_params.soil_depth = vec[model->lgar_bmi_params.num_layers-1];
      
      is_layer_thickness_set = true;

      if (verbosity.compare("high") == 0) {
	std::cout<<"Number of layers : "<<model->lgar_bmi_params.num_layers<<"\n";
	for (int i=1; i<=model->lgar_bmi_params.num_layers; i++)
	  std::cout<<"Thickness, cum. depth : "<<model->lgar_bmi_params.layer_thickness_cm[i]<<" , "<<model->lgar_bmi_params.cum_layer_thickness_cm[i]<<"\n";
	std::cout<<"          *****         \n";
      }
      
      continue;
    }
    else if (param_key == "layer_soil_type") {
      vector<double> vec = ReadVectorData(param_value);
      
      model->lgar_bmi_params.layer_soil_type = new int[vec.size()+1];

      // calculate the cumulative (absolute) depth from land surface to bottom of each soil layer
      model->lgar_bmi_params.cum_layer_thickness_cm[0] = 0;
      
      for (unsigned int layer=1; layer <= vec.size(); layer++) {
      	model->lgar_bmi_params.layer_soil_type[layer] = vec[layer-1];
      }

      is_layer_soil_type_set = true;

      continue;
    }
    else if (param_key == "giuh_ordinates") {
      vector<double> vec = ReadVectorData(param_value);
      
      model->lgar_bmi_params.giuh_ordinates = new double[vec.size()+1];
      
      for (unsigned int layer=1; layer <= vec.size(); layer++) {
      	model->lgar_bmi_params.giuh_ordinates[layer] = vec[layer-1];
      }

      model->lgar_bmi_params.num_giuh_ordinates = vec.size();
      
      is_giuh_ordinates_set = true;
      
      if (verbosity.compare("high") == 0) {
	for (int i=1; i<=model->lgar_bmi_params.num_giuh_ordinates; i++)
	  std::cout<<"GIUH ordinates : "<<model->lgar_bmi_params.giuh_ordinates[i]<<"\n";
	
	std::cout<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "initial_psi") {
      model->lgar_bmi_params.initial_psi_cm = stod(param_value);
      is_initial_psi_set = true;
      
      if (verbosity.compare("high") == 0) {
	std::cout<<"Initial Psi : "<<model->lgar_bmi_params.initial_psi_cm<<"\n";
	std::cout<<"          *****         \n";
      }
      
      continue;
    }
    else if (param_key == "max_soil_types") {
      model->lgar_bmi_params.num_soil_types = stod(param_value);
      is_max_soil_types_set = true;
      continue;
    }
    else if (param_key == "soil_params_file") {
      soil_params_file = param_value;
      is_soil_params_file_set = true;

      if (verbosity.compare("high") == 0) {
	std::cout<<"Soil paramaters file : "<<soil_params_file<<"\n";
	std::cout<<"          *****         \n";
      }
      
      continue;
    }
    else if (param_key == "wilting_point_psi") {
      model->lgar_bmi_params.wilting_point_psi_cm = stod(param_value);
      is_wilting_point_psi_cm_set = true;

      if (verbosity.compare("high") == 0) {
	std::cout<<"Wilting point Psi [cm] : "<<model->lgar_bmi_params.wilting_point_psi_cm<<"\n";
	std::cout<<"          *****         \n";
      }
      
      continue;
    }
    else if (param_key == "timestep") {
      model->lgar_bmi_params.timestep_h = stod(param_value);
      
      if (param_unit == "[s]" || param_unit == "[sec]" || param_unit == "") // defalut time unit is seconds
	model->lgar_bmi_params.timestep_h /= 3600; // convert to hours
      else if (param_unit == "[min]" || param_unit == "[minute]")
	model->lgar_bmi_params.timestep_h /= 60; // convert to hours
      else if (param_unit == "[h]" || param_unit == "[hr]") 
	model->lgar_bmi_params.timestep_h /= 1.0; // convert to hours
      
      assert (model->lgar_bmi_params.timestep_h > 0);
      is_timestep_set = true;

      if (verbosity.compare("high") == 0) {
	std::cout<<"Model timestep [hours,seconds]: "<<model->lgar_bmi_params.timestep_h<<" , "<<model->lgar_bmi_params.timestep_h*3600<<"\n";
	std::cout<<"          *****         \n";
      }
      
      continue;
    }
    else if (param_key == "forcing_resolution") {
      model->lgar_bmi_params.forcing_resolution_h = stod(param_value);

      if (param_unit == "[s]" || param_unit == "[sec]" || param_unit == "") // defalut time unit is seconds
	model->lgar_bmi_params.forcing_resolution_h /= 3600; // convert to hours
      else if (param_unit == "[min]" || param_unit == "[minute]")
	model->lgar_bmi_params.forcing_resolution_h /= 60; // convert to hours
      else if (param_unit == "[h]" || param_unit == "[hr]") 
	model->lgar_bmi_params.forcing_resolution_h /= 1.0; // convert to hours
      
      assert (model->lgar_bmi_params.forcing_resolution_h > 0);
      is_forcing_resolution_set = true;

      if (verbosity.compare("high") == 0) {
	std::cout<<"Forcing resolution [hours]: "<<model->lgar_bmi_params.forcing_resolution_h<<"\n";
	std::cout<<"          *****         \n";
      }
      
      continue;
    }
    
  }
  
  fp.close();

  if(!is_max_soil_types_set)
     model->lgar_bmi_params.num_soil_types = 15;          // maximum number of soil types defaults to 15

  if (verbosity.compare("high") == 0) {
    std::cout<<"Maximum number of soil types: "<<model->lgar_bmi_params.num_soil_types<<"\n";
    std::cout<<"          *****         \n";
  }
   
  if(is_soil_params_file_set) {
    //allocate memory to create an array of structures to hold the soils properties data.
    //model->soil_properties = (struct soil_properties_*) malloc((model->lgar_bmi_params.num_layers+1)*sizeof(struct soil_properties_));

    model->soil_properties = new soil_properties_[model->lgar_bmi_params.num_soil_types+1];
    int num_soil_types = model->lgar_bmi_params.num_soil_types;
    double wilting_point_psi_cm = model->lgar_bmi_params.wilting_point_psi_cm;
    lgar_read_vG_param_file(soil_params_file.c_str(), num_soil_types, wilting_point_psi_cm, model->soil_properties);

    if (verbosity.compare("high") == 0) {
      for (int layer=1; layer<=model->lgar_bmi_params.num_layers; layer++) {
	int soil = model->lgar_bmi_params.layer_soil_type[layer];// layer_soil_type[layer];
	std::cout<<"Soil type/name : "<<model->lgar_bmi_params.layer_soil_type[layer]<<" "<<model->soil_properties[soil].soil_name<<"\n";
      }
	std::cout<<"          *****         \n";
    }
  }  
  
  if (!is_layer_thickness_set) {
    stringstream errMsg;
    errMsg << "layer thickness not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }

  if (!is_initial_psi_set) {
    stringstream errMsg;
    errMsg << "initial psi not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }
  
  if (!is_timestep_set) {
    stringstream errMsg;
    errMsg << "timestep not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }

  if(!is_wilting_point_psi_cm_set) {
    stringstream errMsg;
    errMsg << "wilting point psi not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }
  
  if (!is_forcing_resolution_set) {
    stringstream errMsg;
    errMsg << "forcing resolution not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }
  
  if (!is_layer_soil_type_set) {
    stringstream errMsg;
    errMsg << "layer soil type not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }

  // check if the size of the input data is consistent
  //assert (parameters->ncells > 0);
  

  model->lgar_bmi_params.forcing_interval = int(model->lgar_bmi_params.forcing_resolution_h/model->lgar_bmi_params.timestep_h+1.0e-08); // add 1.0e-08 to prevent truncation error

  InitializeWettingFronts(model->lgar_bmi_params.num_layers, model->lgar_bmi_params.initial_psi_cm, model->lgar_bmi_params.layer_soil_type, model->lgar_bmi_params.cum_layer_thickness_cm, model->soil_properties);

  

  if (verbosity.compare("none") != 0) {
    std::cout<<"--- Initial state/conditions --- \n";
    listPrint();
    std::cout<<"          *****         \n";
  }

  // initial mass in the system
  model->lgar_mass_balance.volstart_cm = lgar_calc_mass_bal(model->lgar_bmi_params.num_layers,model->lgar_bmi_params.cum_layer_thickness_cm);

  model->lgar_bmi_params.ponded_depth_cm = 0.0; // initially we start with a dry surface (no surface ponding)
  model->lgar_bmi_params.nint = 120; // hacked, not needed to be an input option
  model->lgar_bmi_params.num_wetting_fronts = model->lgar_bmi_params.num_layers;

  assert (model->lgar_bmi_params.num_layers == listLength());

  if (verbosity.compare("high") == 0) {
    std::cout<<"Initial ponded depth is set to zero. \n";
    std::cout<<"No. of spatial intervals used in trapezoidal integration to compute G : "<<model->lgar_bmi_params.nint<<"\n";
  }
  
  if (verbosity.compare("none") != 0) {
    std::cout<<"------------- Initialization done! ---------------------- \n";
    std::cout<<"--------------------------------------------------------- \n";
  }
}


//###################################################################
// calculate the initial theta and wilting point moisture content
// ffor the soil in each layer from initial_psi and assumed wilting point psi
// and create initial fronts and include them in each of the soil layers
extern void InitializeWettingFronts(int num_layers, double initial_psi_cm, int *layer_soil_type, double *cum_layer_thickness_cm, struct soil_properties_ *soil_properties)
{
  int soil;
  int front=0;
  double Se, theta_init;
  bool bottom_flag;
  struct wetting_front *current;
  
  for(int layer=1;layer<=num_layers;layer++) {
    front++;
    
    soil = layer_soil_type[layer];
    theta_init = calc_theta_from_h(initial_psi_cm,soil_properties[soil].vg_alpha_per_cm,
				   soil_properties[soil].vg_m,soil_properties[soil].vg_n,
				   soil_properties[soil].theta_e,soil_properties[soil].theta_r);
    
    // the next lines create the initial moisture profile
    bottom_flag=true;  // all intial wetting fronts are in contact with the bottom of the layer they exist in
    // NOTE: The listInsertFront function does lots of stuff.
    current = listInsertFront(cum_layer_thickness_cm[layer],theta_init,front,layer,bottom_flag);
    current->psi_cm = initial_psi_cm;
    Se = calc_Se_from_theta(current->theta,soil_properties[soil].theta_e,soil_properties[soil].theta_r);
    current->K_cm_per_s = calc_K_from_Se(Se, soil_properties[soil].Ksat_cm_per_s, soil_properties[soil].vg_m);  // cm/s
  
  }

}

/*
Reads 1D data from the config file
- used for reading soil discretization (1D)
- used for reading layers depth from the surface if model `layered` is chosen
*/

extern vector<double>
ReadVectorData(string key)
{
  int pos = 0;
  string delimiter = ",";
  vector<double> value(0.0);
  string z1 = key;
  
  if (z1.find(delimiter) == string::npos) {
    double v = stod(z1);
    if (v == 0.0) {
      stringstream errMsg;
      errMsg << "soil_z (depth of soil reservior) should be greater than zero. It it set to "<< v << " in the config file "<< "\n";
      throw runtime_error(errMsg.str());
    }
    
    value.push_back(v);
    
  }
  else {
    while (z1.find(delimiter) != string::npos) {
      pos = z1.find(delimiter);
      string z_v = z1.substr(0, pos);

      value.push_back(stod(z_v.c_str()));
      
      z1.erase(0, pos + delimiter.length());
      if (z1.find(delimiter) == string::npos)
	value.push_back(stod(z1));
    }
  }
  
  return value;
}

/*
extern void lgar_update(struct lgar_model_ *model)
{
  std::cout<<"in lgar update \n";
  listPrint();
  double mm_to_cm = 0.1;

  // local variables for readibility
  int subcycles = model->lgar_bmi_params.forcing_interval;
  
  double precip_timestep_cm;
  double PET_timestep_cm;
  double ponded_depth_cm;
  double AET_timestep_cm;
  double volstart_timestep_cm;
  double volend_timestep_cm = 0.0; // this should not be reset to 0.0 in the for loop
  double volin_timestep_cm;
  double volrunoff_timestep_cm;
  double volrech_timestep_cm;
  double precip_previous_timestep_cm;
  double num_layers;
  double timestep_h = model->lgar_bmi_params.timestep_h;
  int nint = model->lgar_bmi_params.nint;
  double wilting_point_psi_cm = model->lgar_bmi_params.wilting_point_psi_cm;
  double AET_thresh_Theta = 0.85;    // scaled soil moisture (0-1) above which AET=PET (fix later!)
  double AET_expon = 1.0;  // exponent that allows curvature of the rising portion of the Budyko curve (fix later!)

  
  for (int cycle=0; cycle < subcycles; cycle++) {

    state_previous = NULL;
    state_previous = listCopy(head);
    
    precip_timestep_cm = model->lgar_bmi_params.precipitation_cm * mm_to_cm / double(subcycles); // rate; cm/hour
    PET_timestep_cm = model->lgar_bmi_params.PET_cm * mm_to_cm / double(subcycles);
    ponded_depth_cm = precip_timestep_cm * timestep_h;
    AET_timestep_cm = 0.0;
    volstart_timestep_cm = 0.0;
    volin_timestep_cm =0.0;

    volrunoff_timestep_cm = 0.0;
    volrech_timestep_cm = 0.0;
    
    precip_previous_timestep_cm = model->lgar_bmi_params.precip_previous_timestep_cm;
    
    num_layers = model->lgar_bmi_params.num_layers;
    double delta_theta;   // the width of a front, such that its volume=depth*delta_theta
    double dry_depth;
    
    
    if (PET_timestep_cm>0) {
      // Calculate AET from PET and root zone soil moisture.  Note PET was reduced iff raining
      
      AET_timestep_cm = calc_aet(PET_timestep_cm, timestep_h, wilting_point_psi_cm, model->soil_properties, model->lgar_bmi_params.layer_soil_type, AET_thresh_Theta, AET_expon);
    }
    
    
    // put local variables to state timstep variables
    model->lgar_mass_balance.volprecip_timestep_cm = precip_timestep_cm * timestep_h; // volume in cm
    model->lgar_mass_balance.volPET_timestep_cm = PET_timestep_cm;
    
    // put local variables to state global variables
    model->lgar_mass_balance.volprecip_cm += precip_timestep_cm * timestep_h;
    model->lgar_mass_balance.volPET_cm += fmax(PET_timestep_cm,0.0); // ensures non-negative PET
    volstart_timestep_cm = lgar_calc_mass_bal(num_layers,model->lgar_bmi_params.cum_layer_thickness_cm);
    
    
    int soil_num = model->lgar_bmi_params.layer_soil_type[head->layer_num];
    double theta_e = model->soil_properties[soil_num].theta_e;
    bool is_top_wf_saturated = head->theta >= theta_e ? true : false;
    bool create_surficial_front = (precip_previous_timestep_cm == 0.0 && precip_timestep_cm >0.0);
    
    double mass_source_to_soil_timestep = 0.0;
    
    int wf_free_drainage_demand = wetting_front_free_drainage();
    
    //if the follow is true, that would mean there is no wetting front in the top layer to accept the water, must create one.
    if(create_surficial_front && !is_top_wf_saturated)  {
      
      double temp_pd = 0.0; // necessary to assign zero precip due to the creation of new wetting front; AET will still be taken out of the layers
      
      lgar_move_wetting_fronts(&temp_pd,timestep_h, wf_free_drainage_demand, volend_timestep_cm, num_layers, &AET_timestep_cm, model->lgar_bmi_params.cum_layer_thickness_cm, model->lgar_bmi_params.layer_soil_type, model->soil_properties);
      
      dry_depth = lgar_calc_dry_depth(nint, timestep_h, model->lgar_bmi_params.layer_soil_type, model->soil_properties, model->lgar_bmi_params.cum_layer_thickness_cm,&delta_theta);
      
      double theta1 = head->theta;
      lgar_create_surfacial_front(&ponded_depth_cm, &volin_timestep_cm, dry_depth, theta1, model->lgar_bmi_params.layer_soil_type, model->soil_properties, model->lgar_bmi_params.cum_layer_thickness_cm, nint, timestep_h);
      
      state_previous = NULL;
      state_previous = listCopy(head);
      
      model->lgar_mass_balance.volin_cm += volin_timestep_cm;
      
    }

    //listPrint();

    if (ponded_depth_cm > 0 && !create_surficial_front) {
      
      volrunoff_timestep_cm = lgar_insert_water(&ponded_depth_cm, &volin_timestep_cm, precip_timestep_cm, dry_depth, nint, timestep_h, wf_free_drainage_demand, model->lgar_bmi_params.layer_soil_type, model->soil_properties, model->lgar_bmi_params.cum_layer_thickness_cm);
      
      model->lgar_mass_balance.volin_cm += volin_timestep_cm;
      model->lgar_mass_balance.volrunoff_timestep_cm = volrunoff_timestep_cm;
      model->lgar_mass_balance.volrunoff_cm += volrunoff_timestep_cm;
      volrech_timestep_cm = volin_timestep_cm; // this gets updated later, probably not needed here
      model->lgar_mass_balance.volon_cm = ponded_depth_cm;
      //printf("Mass in = %lf %lf %lf \n", volin_timestep_cm, volrech_timestep_cm, volrunoff_timestep_cm);
      if (volrunoff_timestep_cm < 0) abort();  
    }
    else {
      //printf("wetting front created = %lf %d \n", ponded_depth_cm ,!create_surficial_front );
      double hp_cm_max = 0.0; //h_p_max = 0.0;
      
      if (ponded_depth_cm < hp_cm_max) {
	model->lgar_mass_balance.volrunoff_cm += 0.0;
	model->lgar_mass_balance.volon_cm = ponded_depth_cm;
	ponded_depth_cm = 0.0;
	model->lgar_mass_balance.volrunoff_timestep_cm = 0.0;
      }
      else {
	model->lgar_mass_balance.volrunoff_timestep_cm = ponded_depth_cm - hp_cm_max;
	model->lgar_mass_balance.volrunoff_cm += (ponded_depth_cm - hp_cm_max);
	model->lgar_mass_balance.volon_cm = hp_cm_max;
	ponded_depth_cm = hp_cm_max;
	
      }
    }
    
    
    if (!create_surficial_front) {
      lgar_move_wetting_fronts(&volin_timestep_cm,timestep_h, wf_free_drainage_demand, volend_timestep_cm, num_layers, &AET_timestep_cm, model->lgar_bmi_params.cum_layer_thickness_cm, model->lgar_bmi_params.layer_soil_type, model->soil_properties);
      
      // this is the volume of water leaving through the bottom
      volrech_timestep_cm = volin_timestep_cm;
      model->lgar_mass_balance.volrech_timestep_cm = volrech_timestep_cm;
      //printf("Mass in x = %lf %lf \n", volin_timestep_cm, volrech_timestep_cm);
    }
    
    
    int num_dzdt_calculated = lgar_dzdt_calc(nint, model->lgar_bmi_params.layer_soil_type, model->soil_properties, model->lgar_bmi_params.cum_layer_thickness_cm, ponded_depth_cm);
    
    
    model->lgar_mass_balance.volAET_timestep_cm = AET_timestep_cm;
    model->lgar_mass_balance.volAET_cm += AET_timestep_cm;
    model->lgar_mass_balance.volrech_timestep_cm = volrech_timestep_cm;
    model->lgar_mass_balance.volrech_cm += volrech_timestep_cm;
    
    volend_timestep_cm = lgar_calc_mass_bal(num_layers,model->lgar_bmi_params.cum_layer_thickness_cm);
    
    model->lgar_mass_balance.volend_timestep_cm = volend_timestep_cm;
    model->lgar_mass_balance.volend_cm = volend_timestep_cm;
    model->lgar_bmi_params.precip_previous_timestep_cm = precip_timestep_cm;
    
    
    double local_mb = volstart_timestep_cm + model->lgar_mass_balance.volprecip_timestep_cm -  model->lgar_mass_balance.volrunoff_timestep_cm - model->lgar_mass_balance.volAET_timestep_cm - model->lgar_mass_balance.volon_cm - model->lgar_mass_balance.volrech_timestep_cm - volend_timestep_cm;
    
    bool debug_flag = true;
    if(VERBOSE > -1) {
      printf("local mass balance = %0.10e %0.10e %0.10e %0.10e %0.10e %0.10e \n", local_mb, volstart_timestep_cm, model->lgar_mass_balance.volprecip_timestep_cm, volrunoff_timestep_cm,AET_timestep_cm, model->lgar_mass_balance.volend_timestep_cm);
    }
    
    if (fabs(local_mb) >1e-7) {
      printf("Local mass balance (in this timestep) is %.6e, larger than expected, needs some debugging...\n ",local_mb);
      abort();
    }
    
    //listPrint();
    if (head->depth_cm <= 0.0) {
      printf("negative depth = %lf \n", head->depth_cm);
      abort();
    }
  }  
  
}
*/


extern void lgar_global_mass_balance(struct lgar_model_ *model, double *giuh_runoff_queue_cm)
{
  double volstart = model->lgar_mass_balance.volstart_cm;
  double volprecip = model->lgar_mass_balance.volprecip_cm;
  double volrunoff = model->lgar_mass_balance.volrunoff_cm;
  double volAET = model->lgar_mass_balance.volAET_cm;
  double volPET = model->lgar_mass_balance.volPET_cm;
  double volon = model->lgar_mass_balance.volon_cm;
  double volin = model->lgar_mass_balance.volin_cm;
  double volrech = model->lgar_mass_balance.volrech_cm;
  double volend = model->lgar_mass_balance.volend_cm;
  double volend_giuh = model->lgar_mass_balance.volrunoff_giuh_cm;
  double volend_giuh_cm = 0.0;
  double total_Q_cm = model->lgar_mass_balance.volQ_cm;
  //check if the giuh queue have some water left at the end of simulaiton; needs to be included in the global mass balance
  // hold on; this is probably not needed as we have volrunoff in the balance; revist AJK
  for(int i=1; i <= model->lgar_bmi_params.num_giuh_ordinates; i++)
    volend_giuh_cm += giuh_runoff_queue_cm[i];
    
    
  double global_error_cm = volstart + volprecip - volrunoff - volAET - volon -volrech -volend;
  
  
 printf("\n---------------------- Simulation Summary  ------------------------ \n");
 //printf("Time (sec)                 = %6.10f \n", elapsed);
 printf("-------------------------- Mass balance ------------------- \n");
 printf("initial water in soil      = %14.10f cm\n", volstart);
 printf("total precipitation input  = %14.10f cm\n", volprecip);
 printf("total infiltration         = %14.10f cm\n", volin);
 printf("final water in soil        = %14.10f cm\n", volend);
 printf("water remaining on surface = %14.10f cm\n", volon);
 printf("surface runoff             = %14.10f cm\n", volrunoff);
 printf("giuh runoff                = %14.10f cm\n", volrunoff);
 printf("total percolation          = %14.10f cm\n", volrech);
 printf("total AET                  = %14.10f cm\n", volAET);
 printf("total PET                  = %14.10f cm\n", volPET);
 printf("total discharge (Q)        = %14.10f cm\n", total_Q_cm);
 printf("global balance             =   %.6e cm\n", global_error_cm);
}

// calculate precipitation mass that will be added to the layers
extern double potential_infiltration()
{
  // in the main -- this should only be called for the condition (self.precip_data[i]>0)&(self.precip_data[i-1]>0)
  return 0.0;

}

extern int wetting_front_free_drainage() {
  
  int wf_that_supplies_free_drainage_demand = 1;
  struct wetting_front *current;

  //current = head;
  int number_of_wetting_fronts = listLength();
  
  for(current = head; current != NULL; current = current->next) 
  {
    if (current->next != NULL) {
      if (current->layer_num == current->next->layer_num)
	break;
      else
	wf_that_supplies_free_drainage_demand++;

    }
  }

  if (wf_that_supplies_free_drainage_demand > number_of_wetting_fronts)
    wf_that_supplies_free_drainage_demand--;

#if VERBOSE > 1
  printf("wettign_front_free_drainage, layer_num = %d \n", wf_that_supplies_free_drainage_demand);
#endif
  return  wf_that_supplies_free_drainage_demand;
}

// it moves wetting fronts, merge wetting fronts and does the mass balance correction if needed
extern void lgar_move_wetting_fronts(double *ponded_depth_cm, double time_step_s, int wf_free_drainage_demand, double old_mass,
			     int num_layers, double *AET_demand_cm, double *cum_layer_thickness_cm, int *soil_type,
			     struct soil_properties_ *soil_properties)
{

  if (verbosity.compare("high") == 0) {
    printf("State before moving wetting fronts...\n");
    listPrint();
  }
  
  struct wetting_front *current;
  struct wetting_front *next;
  struct wetting_front *previous;

  struct wetting_front *current_old;
  struct wetting_front *next_old;
  struct wetting_front *previous_old;

  struct wetting_front *head_old;
  head_old = listFindFront(1,state_previous);

  double mass_before_move = lgar_calc_mass_bal(0,cum_layer_thickness_cm);
  
  previous = head;
  double theta,Se,theta_e,theta_r;
  double delta_theta;
  double vg_a, vg_m, vg_n, Ks_cm_per_s;
  double psi_cm, K_cm_per_s;
  int layer_num, soil_num;

  // for previous state
  double theta_old, delta_theta_old, K_cm_per_s_old, psi_cm_old;
  int layer_num_old;
  
  int number_of_wetting_fronts = listLength();
  double time_step = time_step_s;
  
  current = head;

  int last_wetting_front_index = number_of_wetting_fronts;
  int layer_number_above, layer_number_below;

  double precip_mass_to_add = (*ponded_depth_cm);

  *ponded_depth_cm = 0.0;

  /* ************************************************************ */
  // main loop advancing all wetting fronts and dooing the mass balance

  //lgar_move_wetting_fronts()
  for (int l = number_of_wetting_fronts; l != 0; l--) {

    if (verbosity.compare("high") == 0) {
      printf("Looping over wetting front = %d %d \n", l);
    }
    
    if (l == 1 && number_of_wetting_fronts >0) {
      current = listFindFront(l,NULL);
      next = listFindFront(l+1,NULL);
      previous = NULL;
      
      current_old = listFindFront(l,state_previous);
      next_old = listFindFront(l+1,state_previous);
      previous_old = NULL;
    }
    else if (l < number_of_wetting_fronts) {
      current = listFindFront(l,NULL);
      next = listFindFront(l+1,NULL);
      previous = listFindFront(l-1,NULL);
      
      current_old = listFindFront(l,state_previous);
      next_old = listFindFront(l+1,state_previous);
      previous_old = listFindFront(l-1,state_previous);
    }
    else if (l == number_of_wetting_fronts) {
      current = listFindFront(l,NULL);
      next = NULL;
      previous = listFindFront(l-1,NULL);
      
      current_old = listFindFront(l,state_previous);
      next_old = NULL;
      previous_old = listFindFront(l-1,state_previous);
    }

    //
    layer_num   = current->layer_num;
    soil_num    = soil_type[layer_num];
    theta_e     = soil_properties[soil_num].theta_e;
    theta_r     = soil_properties[soil_num].theta_r;
    vg_a        = soil_properties[soil_num].vg_alpha_per_cm;
    vg_m        = soil_properties[soil_num].vg_m;
    vg_n        = soil_properties[soil_num].vg_n;
    theta       = current->theta;
    Se          = calc_Se_from_theta(theta, theta_e, theta_r);
    psi_cm      = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
    K_cm_per_s  = calc_K_from_Se(Se, Ks_cm_per_s, vg_m);
    
    layer_number_above = (l == 1) ? layer_num : previous->layer_num;

    layer_number_below = (l == last_wetting_front_index) ? layer_num + 1 : next->layer_num;
    
    if (verbosity.compare("high") == 0) {
       printf ("****************** Cases ***************** \n");
       printf ("Layers (wf, layer, above, below) == %d %d %d %d \n", l ,layer_num, layer_number_above, layer_number_below);
    }
    
    double free_drainage_demand = 0.0;
    double actual_ET_demand = *AET_demand_cm;
    
    // case to check if the wetting front is at the interface, i.e. deepest wetting front within a layer
    // todo. this condition can be replace by current->to_depth = FALSE && l<last_wetting_front_index
    /*************************************************************************************/
    if ( (l<last_wetting_front_index) && (layer_number_below!=layer_num) ) {
      if (verbosity.compare("high") == 0) {
	printf("case (deepest wetting front) : layer_num_below (%d) != layer_num (%d) \n", layer_num, layer_number_below);
	//listPrint();
      }
      
      current->theta = calc_theta_from_h(next->psi_cm, vg_a,vg_m, vg_n, theta_e, theta_r);
      current->psi_cm = next->psi_cm;
    }

    // case to check if the number of wetting fronst are equal to the number of layers, i.e., one wetting front per layer
    /*************************************************************************************/
    
    if (l == number_of_wetting_fronts && layer_number_below != layer_num && number_of_wetting_fronts == num_layers) {

      if (verbosity.compare("high") == 0) {
	printf("case (number_of_wetting_fronts equal to num_layers) : l (%d) == num_layers (%d) == num_wetting_fronts(%d) \n", l, num_layers,number_of_wetting_fronts);
      }
      
      // local variables
      double vg_a_k, vg_m_k, vg_n_k, Ks_cm_per_s_k;
      double theta_e_k, theta_r_k;
      
      current->depth_cm += current->dzdt_cm_per_s * time_step_s; // this is probably not needed, as dz/dt = 0 for the deepest wetting front
      
      double *delta_thetas = (double *) malloc(sizeof(double)*(layer_num+1));
      double *delta_thickness = (double *) malloc(sizeof(double)*(layer_num+1));
      
      double psi_cm_old = current_old->psi_cm;
      double psi_cm_below_old = 0.0;
      
      double psi_cm = current->psi_cm;
      double psi_cm_below = 0.0;
      
      double prior_mass = (current_old->depth_cm - cum_layer_thickness_cm[layer_num-1]) * (current_old->theta-0.0); // 0.0 = next_old->theta
      
      double new_mass = (current->depth_cm - cum_layer_thickness_cm[layer_num-1]) * (current->theta-0.0); // 0.0 = next->theta;
      
      for (int k=1; k<layer_num; k++) {
	int soil_num_k  = soil_type[k];
	theta_e_k = soil_properties[soil_num_k].theta_e;
	theta_r_k = soil_properties[soil_num_k].theta_r;
	vg_a_k    = soil_properties[soil_num_k].vg_alpha_per_cm;
	vg_m_k    = soil_properties[soil_num_k].vg_m;
	vg_n_k    = soil_properties[soil_num_k].vg_n;
	
	double theta_old             = calc_theta_from_h(psi_cm_old, vg_a_k, vg_m_k, vg_n_k, theta_e_k,theta_r_k);
	double theta_below_old       = 0.0;
	double local_delta_theta_old = theta_old - theta_below_old;
	double layer_thickness       = cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k-1];
	
	prior_mass += (layer_thickness * local_delta_theta_old);
	
	double theta       = calc_theta_from_h(psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
	double theta_below = 0.0;
	
	new_mass += layer_thickness * (theta - theta_below);
	delta_thetas[k] = theta_below;
	delta_thickness[k] = layer_thickness;
      }
      
      delta_thetas[layer_num] = 0.0; 
      delta_thickness[layer_num] = current->depth_cm - cum_layer_thickness_cm[layer_num-1];
      
      double free_drainage_demand = 0;
      
      if (wf_free_drainage_demand == l)
	prior_mass += precip_mass_to_add - (free_drainage_demand + actual_ET_demand);
      
      
      double theta_new = lgar_theta_mass_balance(layer_num, soil_num, psi_cm, new_mass, prior_mass, current->depth_cm, delta_thetas, delta_thickness, soil_type, soil_properties);

      current->theta = fmin(theta_new, theta_e);
      
      double Se = calc_Se_from_theta(current->theta,theta_e,theta_r);
      current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
      
    }


    // case to check if the wetting fronst is within the layer
    /*************************************************************************************/
    
    if ( (l < last_wetting_front_index) && (layer_number_below == layer_num) ) {
      
      if (verbosity.compare("high") == 0) {
	printf("case (wetting front within layer) : layer_num (%d) == layer_num_below (%d) \n", layer_num,layer_number_below);
      }
      
      if (layer_num == 1) {
	
	double free_drainage_demand = 0;
	double prior_mass = current_old->depth_cm * (current_old->theta -  next_old->theta);
	
	if (wf_free_drainage_demand == l)
	  prior_mass += precip_mass_to_add - (free_drainage_demand + actual_ET_demand);
	
	current->depth_cm +=  current->dzdt_cm_per_s * time_step_s;
	
	if (current->dzdt_cm_per_s == 0.0 && current->to_bottom == FALSE) // a new front was just created, so don't update it.
	  current->theta = current->theta;
	else
	  current->theta = fmin(theta_e, prior_mass/current->depth_cm + next->theta);
	
      }
      else {
	
	double vg_a_k, vg_m_k, vg_n_k, Ks_cm_per_s_k;
	double theta_e_k, theta_r_k;
	
	current->depth_cm +=  current->dzdt_cm_per_s * time_step_s;
	
	double *delta_thetas = (double *)malloc(sizeof(double)*(layer_num+1));
	double *delta_thickness = (double *)malloc(sizeof(double)*(layer_num+1));
	
	
	double psi_cm_old = current_old->psi_cm;
	double psi_cm_below_old = current_old->next->psi_cm;
	
	double psi_cm = current->psi_cm;
	double psi_cm_below = next->psi_cm;
	
	double prior_mass = (current_old->depth_cm - cum_layer_thickness_cm[layer_num-1]) * (current_old->theta-next_old->theta);
	
	double new_mass = (current->depth_cm - cum_layer_thickness_cm[layer_num-1]) * (current->theta-next->theta);
	
	for (int k=1; k<layer_num; k++) {
	  
	  theta_e_k = soil_properties[soil_num-k].theta_e;
	  theta_r_k = soil_properties[soil_num-k].theta_r;
	  vg_a_k    = soil_properties[soil_num-k].vg_alpha_per_cm;
	  vg_m_k    = soil_properties[soil_num-k].vg_m;
	  vg_n_k    = soil_properties[soil_num-k].vg_n;
	  
	  double theta_old = calc_theta_from_h(psi_cm_old, vg_a_k, vg_m_k, vg_n_k, theta_e_k,theta_r_k);	  
	  double theta_below_old = calc_theta_from_h(psi_cm_below_old, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
	  double local_delta_theta_old = theta_old - theta_below_old;
	  double layer_thickness = (cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k-1]);
	  
	  prior_mass += (layer_thickness * local_delta_theta_old);
	  
	  //-------------------------------------------
	  // do the same for current state
	  double theta = calc_theta_from_h(psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
	  
	  double theta_below = calc_theta_from_h(psi_cm_below, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
	  
	  new_mass += layer_thickness * (theta - theta_below);
	  
	  delta_thetas[k] = theta_below;
	  delta_thickness[k] = layer_thickness;
	}
	
	delta_thetas[layer_num] = next->theta;
	delta_thickness[layer_num] = current->depth_cm - cum_layer_thickness_cm[layer_num-1];
	
	double free_drainage_demand = 0;
	
	if (wf_free_drainage_demand == l)
	  prior_mass += precip_mass_to_add - (free_drainage_demand + actual_ET_demand);
	
	double theta_new = lgar_theta_mass_balance(layer_num, soil_num, psi_cm, new_mass, prior_mass, current->depth_cm, delta_thetas, delta_thickness, soil_type, soil_properties);
	
	current->theta = fmin(theta_new, theta_e);
	
      }
      
      
      double Se = calc_Se_from_theta(current->theta,theta_e,theta_r);
      current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
    }
    
    if (verbosity.compare("high") == 0) {
      printf("*********** Cases for mass balance of wetting fronts are done!! ************** \n");
      listPrint();
    }

    
    // if f_p (predicted infiltration) causes theta > theta_e, mass correction is needed. depth of the wetting front with theta > theta_e is increased to close the mass balance
    // this should be moved out of here to a subroutine; add a call to that subroutine
    if (l == 1) {

      int soil_num_k1  = soil_type[wf_free_drainage_demand];
      double theta_e_k1 = soil_properties[soil_num_k1].theta_e;
	
      struct wetting_front *wf_free_drainage = listFindFront(wf_free_drainage_demand,NULL);
      
      double mass_timestep = (old_mass + precip_mass_to_add) - (actual_ET_demand+free_drainage_demand);

      assert (old_mass > 0.0);
      
      if (wf_free_drainage->theta == theta_e_k1) {

	double current_mass = lgar_calc_mass_bal(0,cum_layer_thickness_cm);
	
	double mass_balance_error = fabs(current_mass - mass_timestep); // mass error
	
	double factor = 1.0;
	bool switched = false;
	double tolerance = 1e-12;
	
	// check if the difference is less than the tolerance
	if (mass_balance_error <= tolerance) {
	  // return current_mass;
	}
	
	double psi_k = wf_free_drainage->depth_cm;
	
	while (mass_balance_error > tolerance) {
	  
	  if (current_mass<mass_timestep) {
	    psi_k = psi_k + 0.01 * factor;
	    switched = false;
	  }
	  else {
	    if (!switched) {
	      switched = true;
	      factor = factor * 0.001;
	    }
	    psi_k = psi_k - 0.01 * factor;
	    
	  }
	  
	  wf_free_drainage->depth_cm = psi_k;
	  
	  current_mass = lgar_calc_mass_bal(0,cum_layer_thickness_cm);
	  mass_balance_error = fabs(current_mass - mass_timestep);
	  
	}	
	
      }
    }
    
    // **************************** MERGING WETTING FRONT ********************************

    if (verbosity.compare("high") == 0) {
      
      if (next != NULL) {
	printf("********** Merging wetting fronts ********** \n");
	//listPrint();
      }
      else
	printf("********** Merging not needed ********** \n");
    }
    
    if (next != NULL) {
      // merge the wetting fronts and returns water leaving through the bottom boundary
      *ponded_depth_cm = lgar_merge_wetting_fronts(num_layers, current, cum_layer_thickness_cm, soil_type, soil_properties);
    }

   
  }
  /*******************************************************************/
  // end of the for loop
  /*******************************************************************/

  // reset current to head to fix any mass balance issues and dry-over-wet wetting fronts conditions
  double mass_change = 0.0;
  
  lgar_fix_wet_over_dry_fronts(&mass_change, cum_layer_thickness_cm, soil_type, soil_properties);

  if (verbosity.compare("high") == 0) {
    printf ("mass change/adjustment (dry_over_wet case) = %lf \n", mass_change);
  }

  
  // adjust AET based on updated mass balance, this ensures to water/mass is lost from the system
  
  if (*AET_demand_cm >0 && mass_change !=0) {
    if (mass_change > 0)
      *AET_demand_cm -= mass_change;
    else {
      *AET_demand_cm += mass_change; // this is not going to happen
      printf("This should not happen!");
      abort();
    }
  }


  // do a general mass check again
  
  double mass_after_move =lgar_calc_mass_bal(0,cum_layer_thickness_cm);

  double mass_correction = mass_before_move+precip_mass_to_add - (mass_after_move+ *AET_demand_cm) ;

  *AET_demand_cm += mass_correction;


  /***********************************************/
  // make sure all psi values are updated
  current = head;
  
  for (int l=1; l != listLength(); l++) {

    int soil_num_k    = soil_type[current->layer_num];
	
    double theta_e_k   = soil_properties[soil_num_k].theta_e;
    double theta_r_k   = soil_properties[soil_num_k].theta_r;
    double vg_a_k      = soil_properties[soil_num_k].vg_alpha_per_cm;
    double vg_m_k      = soil_properties[soil_num_k].vg_m;
    double vg_n_k      = soil_properties[soil_num_k].vg_n;
    double Ksat_cm_per_s_k  = soil_properties[soil_num_k].Ksat_cm_per_s;
    double Se = calc_Se_from_theta(current->theta,theta_e_k,theta_r_k);
    current->psi_cm = calc_h_from_Se(Se, vg_a_k, vg_m_k, vg_n_k);
    
    current->K_cm_per_s = calc_K_from_Se(Se, Ksat_cm_per_s_k, vg_m_k);

    current = current->next;
  }

  if (verbosity.compare("high") == 0) {
    printf("Moving/merging wetting fronts done... \n");
  }
  
}



extern double lgar_merge_wetting_fronts(int num_layers, struct wetting_front *current, double* cum_layer_thickness_cm, int *soil_type, struct soil_properties_ *soil_properties)
{
  double theta,Se,theta_e,theta_r;
  double delta_theta;
  double vg_a, vg_m, vg_n, Ks_cm_per_s;
  double psi_cm, K_cm_per_s;
  int layer_num, soil_num;
  double bottom_flux_cm=0.0;
  
    double column_depth = cum_layer_thickness_cm[num_layers];
    //layer_num = current->layer_num;

    layer_num   = current->layer_num;
    soil_num    = soil_type[layer_num];
    theta_e     = soil_properties[soil_num].theta_e;
    theta_r     = soil_properties[soil_num].theta_r;
    vg_a        = soil_properties[soil_num].vg_alpha_per_cm;
    vg_m        = soil_properties[soil_num].vg_m;
    vg_n        = soil_properties[soil_num].vg_n;
    theta       = current->theta;
    Se          = calc_Se_from_theta(theta, theta_e, theta_r);
    psi_cm      = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
    K_cm_per_s  = calc_K_from_Se(Se, Ks_cm_per_s, vg_m);
    
    struct wetting_front *next = current->next; //listFindFront(l+2,NULL); // next->next;
    struct wetting_front *next_to_next = current->next->next; //listFindFront(l+2,NULL); // next->next;
      
      // case a wetting front passing another wetting front within a layer
      /**********************************************************/
      
      if ( (current->depth_cm > next->depth_cm) && (current->layer_num == next->layer_num) && !next->to_bottom) {

	double current_mass_this_layer = current->depth_cm * (current->theta - next->theta) + next->depth_cm*(next->theta - next_to_next->theta);
	current->depth_cm = current_mass_this_layer / (current->theta - next_to_next->theta);

	double theta_e_k        = soil_properties[soil_num].theta_e;
	double theta_r_k        = soil_properties[soil_num].theta_r;
	double vg_a_k           = soil_properties[soil_num].vg_alpha_per_cm;
	double vg_m_k           = soil_properties[soil_num].vg_m;
	double vg_n_k           = soil_properties[soil_num].vg_n;
	double Ksat_cm_per_s_k  = soil_properties[soil_num].Ksat_cm_per_s;
	double Se_k             = calc_Se_from_theta(current->theta,theta_e_k,theta_r_k);
	current->psi_cm         = calc_h_from_Se(Se_k, vg_a_k, vg_m_k, vg_n_k);
	current->K_cm_per_s     = calc_K_from_Se(Se_k, Ksat_cm_per_s_k, vg_m_k); // AJ - K_temp in python version for 1st layer

#if VERBOSE > 1
	printf ("Deleting wetting front (before)... \n");
	listPrint();
#endif
	listDeleteFront(next->front_num);

#if VERBOSE > 1
	printf ("Deleting wetting front (after) ... \n");
	listPrint();
#endif

      }

      // case wetting front passing the layer depth
      /**********************************************************/
      // or use next_to_next != NULL
      if (current->depth_cm > cum_layer_thickness_cm[layer_num] && current->depth_cm < column_depth) {
	double current_theta = fmin(theta_e, current->theta);
	double overshot_depth = current->depth_cm - next->depth_cm;

	// Given theta return psi_cm, Se or K_cm_per_s etc.
	
	double next_theta_e   = soil_properties[soil_num+1].theta_e;
	double next_theta_r   = soil_properties[soil_num+1].theta_r;
	double next_vg_a      = soil_properties[soil_num+1].vg_alpha_per_cm;
	double next_vg_m      = soil_properties[soil_num+1].vg_m;
	double next_vg_n      = soil_properties[soil_num+1].vg_n;
	double Ksat_cm_per_s  = soil_properties[soil_num+1].Ksat_cm_per_s;

	double Se = calc_Se_from_theta(current->theta,theta_e,theta_r);
	current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);

	current->K_cm_per_s = calc_K_from_Se(Se, Ksat_cm_per_s, vg_m); // AJ check Ksat_cm_per_s here current or next?
	
	double theta_new = calc_theta_from_h(current->psi_cm, next_vg_a, next_vg_m, next_vg_n, next_theta_e, next_theta_r);
	
	double mbal_correction = overshot_depth*(current_theta-next->theta);
	double mbal_Z_correction = mbal_correction/(theta_new - next_to_next->theta);

	double depth_new = cum_layer_thickness_cm[layer_num] + mbal_Z_correction;

	if (depth_new < next_to_next->depth_cm) {
	  next->theta = current->theta;
	  current = listDeleteFront(current->front_num);
	  int num_fronts = listLength();
	  struct wetting_front *current_temp = current;
	  current = listInsertFrontAtDepth(num_fronts, cum_layer_thickness_cm,cum_layer_thickness_cm[layer_num]+mbal_Z_correction, theta_new);
	}
	else {
	  next->theta = current->theta; // assign current theta to next layer theta before deleting the wetting front

	  current = listDeleteFront(current->front_num);

	  current = current->next; // this is points to the next_to_next node (w.r.t the original list)

	  double current_mass_this_layer =  mbal_Z_correction * (theta_new - current->theta) + (current->depth_cm - cum_layer_thickness_cm[layer_num])*(current->theta - current->next->theta);

	  // spread the mass by the space (Mass = Depth * space)
	  // Q : can this depth ever get greater than the next layer depth??
	  current->depth_cm = cum_layer_thickness_cm[layer_num] + current_mass_this_layer / (current->theta - current->next->theta);
	  
	  current->theta = theta_new;

	  int soil_num = soil_type[current->layer_num];
	  
	  double theta_e_k   = soil_properties[soil_num].theta_e;
	  double theta_r_k   = soil_properties[soil_num].theta_r;
	  double vg_a_k      = soil_properties[soil_num].vg_alpha_per_cm;
	  double vg_m_k      = soil_properties[soil_num].vg_m;
	  double vg_n_k      = soil_properties[soil_num].vg_n;
	  double Ksat_cm_per_s_k  = soil_properties[soil_num].Ksat_cm_per_s;
	  
	  double Se_k = calc_Se_from_theta(current->theta,theta_e_k,theta_r_k);
	  current->psi_cm = calc_h_from_Se(Se_k, vg_a_k, vg_m_k, vg_n_k);
	  current->K_cm_per_s = calc_K_from_Se(Se_k, Ksat_cm_per_s_k, vg_m_k);
	}
	
      }
      // case wetting front moving (no merging) within a layer; just update values (check if this is actually needed)
      /**********************************************************/
      else if (current->depth_cm < next->depth_cm && current->to_bottom == false) {
	double theta_e_k   = soil_properties[soil_num].theta_e;
	double theta_r_k   = soil_properties[soil_num].theta_r;
	double vg_a_k      = soil_properties[soil_num].vg_alpha_per_cm;
	double vg_m_k      = soil_properties[soil_num].vg_m;
	double vg_n_k      = soil_properties[soil_num].vg_n;
	double Ksat_cm_per_s_k  = soil_properties[soil_num].Ksat_cm_per_s;

	double Se_k = calc_Se_from_theta(current->theta,theta_e_k,theta_r_k);
	current->psi_cm = calc_h_from_Se(Se_k, vg_a_k, vg_m_k, vg_n_k);
	current->K_cm_per_s = calc_K_from_Se(Se_k, Ksat_cm_per_s_k, vg_m_k);	
      }

      // case wetting front stationary (not moving) and the deepest wetting front in the layer, again just update values
      /**********************************************************/
      else if (current->to_bottom == true) {
	double theta_e_k   = soil_properties[soil_num].theta_e;
	double theta_r_k   = soil_properties[soil_num].theta_r;
	double vg_a_k      = soil_properties[soil_num].vg_alpha_per_cm;
	double vg_m_k      = soil_properties[soil_num].vg_m;
	double vg_n_k      = soil_properties[soil_num].vg_n;
	double Ksat_cm_per_s_k  = soil_properties[soil_num].Ksat_cm_per_s;

	double Se_k = calc_Se_from_theta(current->theta,theta_e_k,theta_r_k);

	current->psi_cm = calc_h_from_Se(Se_k, vg_a_k, vg_m_k, vg_n_k);
	current->K_cm_per_s = calc_K_from_Se(Se_k, Ksat_cm_per_s_k, vg_m_k);
      }
      
      // case wetting front is the deepest one in the last layer (most deepested wetting front in the domain)
      /**********************************************************/
      if (next_to_next == NULL && current->depth_cm > cum_layer_thickness_cm[layer_num]) {

	//  this is the water leaving the system through the bottom of the soil
	//*ponded_depth_cm = (current->theta - next->theta) *  (current->depth_cm - next->depth_cm);
	bottom_flux_cm = (current->theta - next->theta) *  (current->depth_cm - next->depth_cm);
	double theta_e_k   = soil_properties[soil_num].theta_e;
	double theta_r_k   = soil_properties[soil_num].theta_r;
	double vg_a_k      = soil_properties[soil_num].vg_alpha_per_cm;
	double vg_m_k      = soil_properties[soil_num].vg_m;
	double vg_n_k      = soil_properties[soil_num].vg_n;
	double Ksat_cm_per_s_k  = soil_properties[soil_num].Ksat_cm_per_s;
	
	next->theta = current->theta;
	double Se_k = calc_Se_from_theta(current->theta,theta_e_k,theta_r_k);
	next->psi_cm = calc_h_from_Se(Se_k, vg_a_k, vg_m_k, vg_n_k);
	next->K_cm_per_s = calc_K_from_Se(Se_k, Ksat_cm_per_s_k, vg_m_k);
	listDeleteFront(current->front_num);
	
      }

      return bottom_flux_cm;
}


extern void lgar_fix_wet_over_dry_fronts(double *mass_change, double* cum_layer_thickness_cm, int *soil_type, struct soil_properties_ *soil_properties)
{

  struct wetting_front *current;
  struct wetting_front *next;
  current = head;
  next = current->next;
  
  for (int l=1; l != listLength(); l++) {
    
    if (next != NULL) {
      struct wetting_front *next_to_next = listFindFront(l+2,NULL);
      /*
      // this code needs to be tested again for a drier wetting front on top on a wetting front when AET is zero.
      if ( (current->theta < next->theta) && (current->layer_num == next->layer_num)) {
      listPrint();
      printf("A3_theta: merging ............ %d %d %lf %lf %lf \n", l, current->front_num, current->depth_cm, current->theta, next->theta);
      
      // note: in a regular merge the layer passing the one below appears first in the list, however here, in terms of order in the list, the passing layer appears after the the layer that was passed.
      
      //double current_mass_this_layer = current->depth_cm * (current->theta - next->theta) + next->depth_cm*(next->theta - next_to_next->theta);
      //current->depth_cm = current_mass_this_layer / (current->theta - next_to_next->theta);
      
      double current_mass_this_layer = next->depth_cm * (next->theta - current->theta) + current->depth_cm*(current->theta - next_to_next->theta);
      //printf("A1 = %lf %lf %lf \n", next->depth_cm,next->theta, current->theta);
      //printf("A2 = %lf %lf %lf \n", current->depth_cm,current->theta, next_to_next->theta);
      
      
      current->depth_cm = current_mass_this_layer / (next->theta - next_to_next->theta);
      
      current->theta = next->theta;
      
      int soil_numA    = soil_type[current->layer_num];
      
      double theta_e   = soil_properties[soil_numA].theta_e;
      double theta_r   = soil_properties[soil_numA].theta_r;
      double vg_a      = soil_properties[soil_numA].vg_alpha_per_cm;
      double vg_m      = soil_properties[soil_numA].vg_m;
      double vg_n      = soil_properties[soil_numA].vg_n;
      double Ksat_cm_per_s  = soil_properties[soil_numA].Ksat_cm_per_s;
      double Se = calc_Se_from_theta(current->theta,theta_e,theta_r);
      current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
      
      //current->K_cm_per_s = soil_properties[soil_num].Ksat_cm_per_s;
      current->K_cm_per_s = calc_K_from_Se(Se, Ksat_cm_per_s, vg_m); // AJ - K_temp in python version for 1st layer
      printf("A40 \n");
      listDeleteFront(next->front_num);
      //listPrint();
      //abort();
      
      }*/
      
      
      // this part fixes thetas when AET caused upper wetting front drier than the lower wetting front
      /***************************************************/
      
      if ( (current->theta < next->theta) && (current->layer_num == next->layer_num) ) {
	int layer_num_k = current->layer_num;
	double mass_before = lgar_calc_mass_bal(0,cum_layer_thickness_cm);
	
	current = listDeleteFront(current->front_num);
	// this needs to be revised
	if (layer_num_k > 1) {
	  int soil_num_k     = soil_type[current->layer_num];
	  double theta_e_k   = soil_properties[soil_num_k].theta_e;
	  double theta_r_k   = soil_properties[soil_num_k].theta_r;
	  double vg_a_k      = soil_properties[soil_num_k].vg_alpha_per_cm;
	  double vg_m_k      = soil_properties[soil_num_k].vg_m;
	  double vg_n_k      = soil_properties[soil_num_k].vg_n;
	  double Ksat_cm_per_s_k  = soil_properties[soil_num_k].Ksat_cm_per_s;
	  double Se_k = calc_Se_from_theta(current->theta,theta_e_k,theta_r_k);
	  
	  current->psi_cm = calc_h_from_Se(Se_k, vg_a_k, vg_m_k, vg_n_k);
	  
	  struct wetting_front *current_local = head;
	  
	  while (current_local->layer_num < layer_num_k) {
	    int soil_num_k1    = soil_type[current_local->layer_num];
	    theta_e_k   = soil_properties[soil_num_k1].theta_e;
	    theta_r_k   = soil_properties[soil_num_k1].theta_r;
	    vg_a_k      = soil_properties[soil_num_k1].vg_alpha_per_cm;
	    vg_m_k      = soil_properties[soil_num_k1].vg_m;
	    vg_n_k      = soil_properties[soil_num_k1].vg_n;
	    double Se_l = calc_Se_from_theta(current->theta,theta_e_k,theta_r_k);
	    
	    current_local->psi_cm = calc_h_from_Se(Se_l, vg_a_k, vg_m_k, vg_n_k);
	    current_local->theta = calc_theta_from_h(current->psi_cm, vg_a_k, vg_m_k, vg_n_k,theta_e_k,theta_r_k);
	    current_local = current_local->next;
	    }
	}
	
	double mass_after=lgar_calc_mass_bal(0,cum_layer_thickness_cm);
	*mass_change += fabs(mass_after - mass_before);
	
	// note: mass_before is less when we have wetter front over drier front condition, however, lgar_calc_mass_bal returns mass_before > mass_after due to fabs(theta_current - theta_next); for mass_before the functions compuates more than the actual mass; removing fabs in that function might be one option, but for now we are adding fabs to mass_change to make sure we added extra water back to AET after deleting the drier front
	//printf("Mass change = %6.15f %6.15f %6.15f \n", mass_change, mass_after, mass_before);
	
      }
      current = current->next;
      next = current->next;
    }
  }

}

extern double lgar_insert_water(double *ponded_depth_cm, double *volin_this_timestep, double precip_timestep_cm,
				double dry_depth, int nint, double time_step_s, int wf_free_drainage_demand,
				int *soil_type, struct soil_properties_ *soil_properties, double *cum_layer_thickness_cm)
{

  // This subroutine inserts precipitation into the soil
  // note ponded_depth_cm is a pointer.   Access it's value as (*ponded_depth_cm).
  
  int wf_that_supplies_free_drainage_demand = wf_free_drainage_demand;
  
  // local vars
  double theta, theta_e,Se,theta_r;
  double delta_theta;
  double psi_cm, K_cm_per_s;
  double vg_a, vg_m, vg_n,Ksat_cm_per_s;
  double h_min_cm;
  bool debug_flag=TRUE;
  struct wetting_front *current;
  struct wetting_front *current_free_drainage;
  struct wetting_front *current_free_drainage_next;
  int layer_num,soil_num;
  double f_p = 0.0;
  double runoff = 0.0;

  double hp_cm_max=0.0;

  double h_p = fmax(*ponded_depth_cm - precip_timestep_cm * time_step_s,0.0); // water ponded on the surface
  
  current = head;
  current_free_drainage = listFindFront(wf_that_supplies_free_drainage_demand,NULL);
  current_free_drainage_next = listFindFront(wf_that_supplies_free_drainage_demand+1,NULL);

  int number_of_wetting_fronts = listLength();
  double time_step = time_step_s;
  
  int l = number_of_wetting_fronts;
  int last_wetting_front_index = number_of_wetting_fronts;
  int num_layers = 3; // hacked
  int layer_num_fp = current_free_drainage->layer_num;

  
  double Geff;
  
  if (number_of_wetting_fronts == num_layers) {
    Geff = 0.0;
  }
  else {
   
    double theta = current_free_drainage->theta;
    double theta_below = current_free_drainage_next->theta;
    
    soil_num = soil_type[layer_num_fp];
      
    theta_e = soil_properties[soil_num].theta_e;  // rhs of the new front, assumes theta_e as per Peter
    theta_r = soil_properties[soil_num].theta_r;
    h_min_cm = soil_properties[soil_num].h_min_cm;
    vg_a     = soil_properties[soil_num].vg_alpha_per_cm;
    vg_m     = soil_properties[soil_num].vg_m;
    vg_n     = soil_properties[soil_num].vg_n;
    Ksat_cm_per_s = soil_properties[soil_num].Ksat_cm_per_s;

    Se = calc_Se_from_theta(theta,theta_e,theta_r);
    psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
  
    Geff = calc_Geff(theta_below, theta_e, theta_e, theta_r, vg_a, vg_n, vg_m, h_min_cm, Ksat_cm_per_s, nint);

  }
   
  if (layer_num_fp == 1) {
      f_p = Ksat_cm_per_s * (1 + (Geff+ h_p)/current_free_drainage->depth_cm);
  }
  else {
    
    double bottom_sum = (current_free_drainage->depth_cm-cum_layer_thickness_cm[layer_num_fp-1])/Ksat_cm_per_s;
    double xx = bottom_sum;
    
    for (int k = 1; k < layer_num_fp; k++) {
      int soil_num_k = soil_type[layer_num_fp-k];
      double Ksat_cm_per_s_k = soil_properties[soil_num_k].Ksat_cm_per_s;
      
      bottom_sum += (cum_layer_thickness_cm[layer_num_fp - k] - cum_layer_thickness_cm[layer_num_fp - (k+1)])/ Ksat_cm_per_s_k;
    }
    
    f_p = (current_free_drainage->depth_cm + Geff + h_p) / bottom_sum;
    
  }
  
  // checkpoint # AJ
  int soil_num_k = soil_type[head->layer_num];
  double theta_e1 = soil_properties[soil_num_k].theta_e; // saturated theta of top layer
  
  if ((layer_num_fp==num_layers) && (current_free_drainage->theta == theta_e1) && (num_layers == number_of_wetting_fronts)) 
    f_p = 0.0;

  double ponded_depth_temp = *ponded_depth_cm;
    
  double free_drainage_demand = 0;
  
  // checkpoint # AJ
  if ((layer_num_fp==num_layers) && (num_layers == number_of_wetting_fronts))
    ponded_depth_temp = *ponded_depth_cm - f_p*time_step - free_drainage_demand*0;
  else
    ponded_depth_temp = *ponded_depth_cm - f_p*time_step - free_drainage_demand*0;
    
  ponded_depth_temp = fmax(ponded_depth_temp, 0.0);
  
  double fp_cm = f_p * time_step + free_drainage_demand/time_step; // infiltration in cm
  
  double precip_mass_to_add;

  if (hp_cm_max > 0.0 ) {
    
    if (ponded_depth_temp < hp_cm_max) {
      runoff = 0.0;
      *volin_this_timestep = fmin(*ponded_depth_cm, fp_cm);
      *ponded_depth_cm = *ponded_depth_cm - *volin_this_timestep;
      return runoff;
    }
    else if (ponded_depth_temp > hp_cm_max ) {
      runoff = ponded_depth_temp - hp_cm_max;
      *ponded_depth_cm = hp_cm_max;
      *volin_this_timestep=fp_cm;
      
      return runoff;
    }

  }

  else {
    // if it got to this point, no ponding is allowed, either infiltrate or runoff
    // order is important here; assign zero to ponded depth once done using it computing volume in and runoff
    *volin_this_timestep = fmin(*ponded_depth_cm, fp_cm); //
    runoff = *ponded_depth_cm < fp_cm ? 0.0 : (*ponded_depth_cm - *volin_this_timestep);
    *ponded_depth_cm = 0.0;
    
  }

  return runoff;
}

extern void lgar_create_surfacial_front(double *ponded_depth_cm, double *volin, double dry_depth, double theta1,
                                        int *soil_type, struct soil_properties_ *soil_properties, 
                                        double *cum_layer_thickness_cm, int nint,double time_step_s)
{
  // This subroutine is called iff there is no surfacial front, it creates a new front and inserts ponded depth
  // into the soil.  Note ponded_depth_cm is a pointer.   Access it's value as (*ponded_depth_cm).

  // local vars
  double theta_e,Se,theta_r;
  double delta_theta,hp_cm;
  double vg_alpha_per_cm, vg_m, vg_n, Ksat_cm_per_s;
  double h_min_cm;
  double Geff;
  
  bool to_bottom = FALSE;
  struct wetting_front *current;
  int layer_num,soil_num,front_num;
  
  current = head;
  
  layer_num = 1;   // we only create new surfacial fronts in the first layer
  soil_num = soil_type[layer_num];
  front_num = 1;   // we are creating a new surfacial front, which by definition must be front #1
  
  theta_e = soil_properties[soil_num].theta_e;  // rhs of the new front, assumes theta_e as per Peter
  theta_r = soil_properties[soil_num].theta_r;
  delta_theta =  theta_e - theta1;
  h_min_cm = soil_properties[soil_num].h_min_cm;
  
  double theta_new;

  if(dry_depth * delta_theta > (*ponded_depth_cm))  // all the ponded depth enters the soil
    {
      *volin = *ponded_depth_cm;
      theta_new = fmin(theta1 + (*ponded_depth_cm) /dry_depth, theta_e);
      listInsertFirst(dry_depth, theta_new, front_num, layer_num, to_bottom);
      *ponded_depth_cm = 0.0;
      hp_cm =0.0;
    }
  else  // not all ponded depth fits in
    {
      *volin = dry_depth * delta_theta;
      *ponded_depth_cm -= dry_depth * delta_theta;
      theta_new = theta_e; //fmin(theta1 + (*ponded_depth_cm) /dry_depth, theta_e);
      if (dry_depth < cum_layer_thickness_cm[1])
	listInsertFirst(dry_depth, theta_e, front_num, layer_num, to_bottom);
      else
	listInsertFirst(dry_depth, theta_e, front_num, layer_num, 1);
      hp_cm = *ponded_depth_cm;
    } 
  
  current = head;  // must do this again because listInsertFirst() created a new *head
  vg_alpha_per_cm    = soil_properties[soil_num].vg_alpha_per_cm;
  vg_m               = soil_properties[soil_num].vg_m;
  vg_n               = soil_properties[soil_num].vg_n;
  Ksat_cm_per_s      = soil_properties[soil_num].Ksat_cm_per_s;

  Se = calc_Se_from_theta(theta_new,theta_e,theta_r);
  current->psi_cm = calc_h_from_Se(Se, vg_alpha_per_cm , vg_m, vg_n);

  current->K_cm_per_s = calc_K_from_Se(Se, Ksat_cm_per_s, vg_m); // AJ - K_temp in python version for 1st layer

  Geff = calc_Geff(theta_new, theta1, theta_e, theta_r, vg_alpha_per_cm, vg_n, vg_m, h_min_cm, Ksat_cm_per_s, nint);
  h_min_cm = 0.0;

  current->dzdt_cm_per_s = 0.0; //for now assign 0 to dzdt as it will be computed in compute_dzdt function
  
  return;
  
}


extern double lgar_calc_dry_depth(int nint, double time_step_s, int *soil_type, 
                                  struct soil_properties_ *soil_properties, double *cum_layer_thickness_cm,
                                  double *delta_theta)
{
  // This routine calculates the "dry depth" of a newly created wetting front in the top soil layer after
  // a non-rainy period or a big increase in rainrate  on an unsaturated first layer.
  // Note: Calculation of the initial depth of a new wetting front in the first layer uses the concept of "dry depth",
  // described in the 2015 GARTO paper (Lai et al., An efficient and guaranteed stable numerical method ffor
  // continuous modeling of infiltration and redistribution with a shallow dynamic water table).
    
  // local variables
  
  struct wetting_front *current;        
  double theta1,theta2,theta_e,theta_r; 
  double vg_alpha_per_cm,vg_n,vg_m,Ksat_cm_per_s,h_min_cm; 
  double time_step_h;                   
  double tau;                           
  double Geff;                          
  double dry_depth;                     
  int    soil_num;                          
  int    layer_num;                   
  
  current=head;                         
  
  layer_num  = current->layer_num;
  soil_num   = soil_type[layer_num];
  
  // copy values of soil properties into shorter variable names to improve readability
  theta_r         = soil_properties[soil_num].theta_r;
  vg_alpha_per_cm = soil_properties[soil_num].vg_alpha_per_cm; 
  vg_m            = soil_properties[soil_num].vg_m;
  vg_n            = soil_properties[soil_num].vg_n;
  Ksat_cm_per_s   = soil_properties[soil_num].Ksat_cm_per_s;
  h_min_cm        = soil_properties[soil_num].h_min_cm;

  // these are the limits of integration
  theta1   = current->theta;                 // water content of the first (most surficial) existing wetting front
  theta_e  = soil_properties[soil_num].theta_e;
  theta2 = theta_e;
  
  *delta_theta = theta_e - current->theta;  // return the delta_theta value to the calling function
  
  tau       = time_step_s * Ksat_cm_per_s/(theta_e-current->theta); //3600
  
  //if(theta1>theta_e) {fred=1;}
  
  Geff      = calc_Geff(theta1, theta2, theta_e, theta_r, vg_alpha_per_cm, vg_n, vg_m, h_min_cm, Ksat_cm_per_s, nint);
 
  // note that dry depth originally has a factor of 0.5 in front
  dry_depth = 0.5*(tau + sqrt( tau*tau + 4.0*tau*Geff) ); 

  //when dry depth greater than layer 1 thickness, set dry depth to layer 1 thickness
  dry_depth = fmin(cum_layer_thickness_cm[layer_num], dry_depth);
 
 return dry_depth;

}

double lgar_calc_mass_bal(int num_soil_layers, double *cum_layer_thickness)
{
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/* function to calculate the amount of soil moisture  in the profile (cm) */
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

struct wetting_front* current;
struct wetting_front* previous;
struct wetting_front* next;     // Bware, might get a little confusing, because there exists a next->next

double sum=0.0;
double base_depth;
int layer,prev_layer,front;

if(head == NULL) return 0.0;

current=head;
prev_layer=ONE;
do
  {
  layer=current->layer_num;
  base_depth=cum_layer_thickness[layer-1];   // note cum_layer_thickness[0]=0.0;
  
  if(current->next != NULL) {            // this is not the last entry in the list
    next=current->next;
    if(next->layer_num == current->layer_num) {
      sum += (current->depth_cm - base_depth) * (current->theta - next->theta); // note no need for fabs() here otherwise we get more mass for the case dry-over-wet front within a layer
    }
    else {
      sum += (current->depth_cm - base_depth) * current->theta;
    }
  }
  else { // this is the last entry in the list.  This must be the deepest front in the final layer
    layer=current->layer_num;
    sum+=current->theta * (current->depth_cm - base_depth);
  }
  
  current=current->next;
  
  } while(current != NULL );   // putting conditional at end of do looop makes sure it executes once.

 return sum;
}

extern void lgar_read_vG_param_file(char const* vG_param_file_name, int num_soil_types, double wilting_point_psi_cm,
				    struct soil_properties_ *soil_properties)
{
  //###################################################################
  // READ THE SOIL PARAMS ******************************************
  // OPEN FILE TO READ IN THE vG parameters ffor standard soil types
  //###################################################################

  if (verbosity.compare("high") == 0) {
    std::cout<<"Reading Van Genuchton parameters files...\n";
  }
  
  // local vars
  FILE *in_vG_params_fptr = NULL;
  char jstr[256];
  char soil_name[30];
  bool error;
  int length;
  int soil;             // soil counter
  double theta_e,theta_r,vg_n,vg_m,vg_alpha_per_cm,Ksat_cm_per_h;  // shorthand variable names
  double m,p,lambda;
  
  // open the file
  if((in_vG_params_fptr=fopen(vG_param_file_name,"r"))==NULL)
    {printf("Can't open input file named %s. Program stopped.\n",vG_param_file_name); exit(-1);}
  
  fgets(jstr,255,in_vG_params_fptr);   // read the header line and ignore
  
  for(soil=1;soil<=num_soil_types;soil++)  // read the num_soil_types lines with data
    {

      fgets(jstr,255,in_vG_params_fptr);
      sscanf(jstr,"%s %lf %lf %lf %lf %lf %lf",soil_name,&theta_r,&theta_e,&vg_alpha_per_cm,&vg_n,&vg_m,&Ksat_cm_per_h);
      length=strlen(soil_name);
      
      if(length>MAX_SOIL_NAME_CHARS)
	{
	  printf("While reading vG soil parameter file: %s, soil name longer than allowed.  Increase MAX_SOIL_NAME_CHARS\n",
		 vG_param_file_name); 
	  printf("Program stopped.\n");
	  exit;
	}
      strcpy(soil_properties[soil].soil_name,soil_name);
      soil_properties[soil].theta_r         = theta_r;
      soil_properties[soil].theta_e         = theta_e;
      soil_properties[soil].vg_alpha_per_cm = vg_alpha_per_cm; // cm^(-1)
      soil_properties[soil].vg_n            = vg_n;
      soil_properties[soil].vg_m            = vg_m;
      soil_properties[soil].Ksat_cm_per_s   = Ksat_cm_per_h;

      soil_properties[soil].theta_wp = calc_theta_from_h(wilting_point_psi_cm, vg_alpha_per_cm,
							 vg_m, vg_n, theta_e, theta_r);
      
      // Given van Genutchen parameters calculate estimates of Brooks & Corey bc_lambda and bc_psib
      if (1.0 < vg_n)
	{
	  m = 1.0 - 1.0 / vg_n;
	  p = 1.0 + 2.0 / m;
	  soil_properties[soil].bc_lambda  = 2.0 / (p - 3.0);
	  soil_properties[soil].bc_psib_cm = (p + 3.0) * (147.8 + 8.1 * p + 0.092 * p * p) / 
	    (2.0 * vg_alpha_per_cm * p * (p - 1.0) * (55.6 + 7.4 * p + p * p));
	  assert(0.0 < soil_properties[soil].bc_psib_cm);
	}
      else
	{
	  fprintf(stderr, "ERROR: Van Genutchen parameter n must be greater than 1\n");
	  error = TRUE;  // TODO FIXME - what todo in this ccase?
	}
      
      /* this is the effective capillary drive after */
      /* Morel-Seytoux et al. (1996) eqn. 13 or 15 */
      /* psi should not be less than this value.  */
      lambda=soil_properties[soil].bc_lambda;
      soil_properties[soil].h_min_cm = soil_properties[soil].bc_psib_cm*(2.0+3.0/lambda)/(1.0+3.0/lambda);                                                             
      
    }
  fclose(in_vG_params_fptr);      // close the file, we're done with it

  return;
}


extern int lgar_dzdt_calc(int nint, int *soil_type, struct soil_properties_ *soil_properties, 
                          double *cum_layer_thickness_cm, double h_p)   // called derivs() in python
{

/****************************************/
/* code to calculate velocity of fronts */
/****************************************/
  int    nfronts_analyzed = 0;    // this is the return value

  struct wetting_front* current;
  struct wetting_front* next;
  struct wetting_front* previous = NULL;
  
  double Se,theta,vg_alpha_per_cm,vg_n,vg_m,Ksat_cm_per_s,theta_e,theta_r;  // local variables to make things clearer
  double delta_theta;
  double Geff;
  double depth_cm;    // the absolute depth down to a wetting front from the surface
  double h_min_cm;
  double K_cm_per_s;  // unsaturated hydraulic conductivity K(theta) at the RHS of the current wetting front (cm/s)
  double h;           // capillary head psi(theta) at the right-hand side of the current wetting front (cm)
  double theta1, theta2;  // limits of integration on Geff from theta1 to theta2
  double bottom_sum;  // store a running sum of L_n/K(theta_n) n increasing from 1 to N-1, as we go down in layers N
  double Z_cm;        // the depth that a wetting front has advanced in layer n
  double dzdt;
  int    soil_num, layer_num, front_num;
  int    prev_layer_num=0;      // so we can tell when the computation moves into a new layer
  
  
  if(head == NULL) {
    return nfronts_analyzed;  // No wetting fronts!!
  }
  
  // make sure to use previous state values as current state is updated during the timestep (that's how it is done is Peter's python version)
  
  current = head;
  
  do {  // loop through the wetting fronts
    dzdt = 0.0;
    nfronts_analyzed++;
    
    // copy structure elements into shorter variables names to increase readability
    // WETTING FRONT PROPERTIES
    front_num    = current->front_num;    // the front number
    theta        = current->theta;        // water content of this front
    layer_num    = current->layer_num;    // what layer the front is in
    K_cm_per_s   = current->K_cm_per_s;   // K(theta)
    
    if (K_cm_per_s <=0) {
      printf("K is zero: %d %lf \n", layer_num, K_cm_per_s);
      abort();
    }
    
    depth_cm = current->depth_cm;     // absolute Z to this wetting front measured down from land surface
    
    // SOIL PROPERTIES
    soil_num        = soil_type[layer_num];
    vg_alpha_per_cm = soil_properties[soil_num].vg_alpha_per_cm;
    vg_n            = soil_properties[soil_num].vg_n;
    vg_m            = soil_properties[soil_num].vg_m;
    theta_e         = soil_properties[soil_num].theta_e;
    theta_r         = soil_properties[soil_num].theta_r;
    h_min_cm        = soil_properties[soil_num].h_min_cm;
    Ksat_cm_per_s   = soil_properties[soil_num].Ksat_cm_per_s;
    
    if(front_num == 1)  prev_layer_num = layer_num;  // first front must be in first layer
    
    next = current->next;    // the next element in the linked list
    if (next == NULL) break; // we're done calculating dZ/dt's because we're at the end of the list
    
    theta1 = next->theta;
    theta2 = current->theta;
    
    
    bottom_sum = 0.0;  // needed ffor multi-layered dz/dt equation.  Equal to sum from n=1 to N-1 of (L_n/K_n(theta_n))
    
    if(current->to_bottom == TRUE) { // checkpoint # AJ
      if(layer_num > 1)
	current->dzdt_cm_per_s = 0.0;
      else
	current->dzdt_cm_per_s = 0.0;
      
      previous = current;
      current = current->next;  // point to the next link
      continue;                 // go to next front, this one fully penetrates the layer
    }
    else {
      
      if(layer_num > 1) {
	double K_cm_per_s_prev = previous->K_cm_per_s; 
	bottom_sum += (current->depth_cm-cum_layer_thickness_cm[layer_num-1])/K_cm_per_s;
      }
      
    }
    
    if(theta1 > theta2) {
      printf("Calculating dzdt : theta1 > theta2 = (%lf, %lf) ... aborting \n", theta1, theta2);  // this should never happen
      abort();
    }
    
    Geff = calc_Geff(theta1, theta2, theta_e, theta_r, vg_alpha_per_cm, vg_n, vg_m, h_min_cm, Ksat_cm_per_s, nint);
    delta_theta = current->theta - next->theta;
    
    if(current->layer_num == 1) { // this front is in the upper layer	  
	dzdt = 1.0/delta_theta*(Ksat_cm_per_s*(Geff+h_p)/current->depth_cm+current->K_cm_per_s);
    }
    else {  // we are in the second or greater layer
      double denominatorA = bottom_sum;
      
      for (int k = 1; k < layer_num; k++) {
	int soil_numA = soil_type[layer_num-k];
	double theta_prevA = calc_theta_from_h(current->psi_cm, soil_properties[soil_numA].vg_alpha_per_cm, soil_properties[soil_numA].vg_m, soil_properties[soil_numA].vg_n,soil_properties[soil_numA].theta_e,soil_properties[soil_numA].theta_r);
	
	
	double Se_prevA = calc_Se_from_theta(theta_prevA,soil_properties[soil_numA].theta_e,soil_properties[soil_numA].theta_r);
	
	double K_cm_per_s_prevA = calc_K_from_Se(Se_prevA,soil_properties[soil_numA].Ksat_cm_per_s, soil_properties[soil_numA].vg_m);
	
	denominatorA += (cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k-1])/ K_cm_per_s_prevA;
	
      }
      
      double theta_prev = calc_theta_from_h(current->psi_cm, soil_properties[soil_num-1].vg_alpha_per_cm, soil_properties[soil_num-1].vg_m, soil_properties[soil_num-1].vg_n,soil_properties[soil_num-1].theta_e,soil_properties[soil_num-1].theta_r);
      
      double Se_prev = calc_Se_from_theta(theta_prev,soil_properties[soil_num-1].theta_e,soil_properties[soil_num-1].theta_r);
      
      double K_cm_per_s_prev = calc_K_from_Se(Se_prev,soil_properties[soil_num-1].Ksat_cm_per_s, soil_properties[soil_num-1].vg_m);
      
      
      Z_cm = depth_cm - cum_layer_thickness_cm[layer_num -1];
      
      double numerator = depth_cm + (Geff +h_p)* Ksat_cm_per_s / K_cm_per_s;
      double denominator = cum_layer_thickness_cm[layer_num -1] / K_cm_per_s_prev + bottom_sum;
      
      dzdt = 1.0/delta_theta * numerator / denominatorA;
    }
    
    previous = current;
    current->dzdt_cm_per_s = dzdt;
    
    current = current->next;  // point to the next link
    
  } while(current != NULL );   // putting conditional at end of do looop makes sure it executes at least once
  
  
  return nfronts_analyzed;
}


extern double lgar_theta_mass_balance(int layer_num, int soil_num, double psi_cm, double new_mass, double prior_mass, double depth_cm_old, double *delta_theta, double *delta_thickness, int *soil_type, struct soil_properties_ *soil_properties)
{

  double x = psi_cm;
  double delta_mass = fabs(new_mass - prior_mass);
  double tolerance = 1e-12;

  double factor = 1.0;// * delta_mass;
  bool switched = false;

  double theta_of_wf_in_layer1_above, theta_of_wf_in_layer2_above;
  double new_mass_in_layer1_above, new_mass_in_layer2_above, new_mass_in_same_layer;
  double theta = 0; // this will be updated and returned


  // check if the difference is less than the tolerance
  if (delta_mass <= tolerance) {
    theta = calc_theta_from_h(x, soil_properties[soil_num].vg_alpha_per_cm, soil_properties[soil_num].vg_m, soil_properties[soil_num].vg_n,soil_properties[soil_num].theta_e,soil_properties[soil_num].theta_r);
    return theta;
  }

  while (delta_mass > tolerance) {
    
    if (new_mass>prior_mass) {
      x = x + 0.1 * factor;
      switched = false;
    }
    else {
      if (!switched) {
	switched = true;
	factor = factor * 0.1;
      }
      x = x - 0.1 * factor;

    }

    //x = fmax(x,0.0);
    double theta_layer;
    double mass_layers= 0.0;
    
    theta = calc_theta_from_h(x, soil_properties[soil_num].vg_alpha_per_cm, soil_properties[soil_num].vg_m, soil_properties[soil_num].vg_n,soil_properties[soil_num].theta_e,soil_properties[soil_num].theta_r);
    mass_layers += delta_thickness[layer_num] * (theta - delta_theta[layer_num]);
    
    for (int k=1; k<layer_num; k++) {
      int soil_numA =  soil_type[k];

      theta_layer = calc_theta_from_h(x, soil_properties[soil_numA].vg_alpha_per_cm, soil_properties[soil_numA].vg_m, soil_properties[soil_numA].vg_n,soil_properties[soil_numA].theta_e,soil_properties[soil_numA].theta_r);

      mass_layers += delta_thickness[k] * (theta_layer - delta_theta[k]);
    }

    new_mass = mass_layers;
    delta_mass = fabs(new_mass - prior_mass);
    
  }
 
  return theta;
  
}

#endif
