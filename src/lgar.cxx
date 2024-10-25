#ifndef LGAR_CXX_INCLUDED
#define LGAR_CXX_INCLUDED

#include "../include/all.hxx"
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>

using namespace std;


//#####################################################################################
/* authors : Ahmad Jan, Fred Ogden, and Peter La Follette
   year    : 2022
   email   : ahmad.jan@noaa.gov, peter.lafollette@noaa.gov
   - The file constains lgar subroutines */
//#####################################################################################


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



// ############################################################################################
/*
  @param soil_moisture_wetting_fronts  : 1D array of thetas (soil moisture content) per wetting front;
                                         output to other models (e.g. soil freeze-thaw)
  @param soil_depth_wetting_fronts : 1D array of absolute depths of the wetting fronts [meters];
					 output to other models (e.g. soil freeze-thaw)
*/
// ############################################################################################
extern void lgar_initialize(string config_file, struct model_state *state)
{
  int soil;
  
  InitFromConfigFile(config_file, state);
  state->lgar_bmi_params.shape[0] = state->lgar_bmi_params.num_layers;
  state->lgar_bmi_params.shape[1] = state->lgar_bmi_params.num_wetting_fronts;

  // initial number of wetting fronts are same are number of layers
  state->lgar_bmi_params.num_wetting_fronts           = state->lgar_bmi_params.num_layers;
  state->lgar_bmi_params.soil_depth_wetting_fronts    = new double[state->lgar_bmi_params.num_wetting_fronts];
  state->lgar_bmi_params.soil_moisture_wetting_fronts = new double[state->lgar_bmi_params.num_wetting_fronts];

  // initialize array for holding calibratable parameters
  state->lgar_calib_params.theta_e  = new double[state->lgar_bmi_params.num_layers];
  state->lgar_calib_params.theta_r  = new double[state->lgar_bmi_params.num_layers];
  state->lgar_calib_params.vg_n     = new double[state->lgar_bmi_params.num_layers];
  state->lgar_calib_params.vg_alpha = new double[state->lgar_bmi_params.num_layers];
  state->lgar_calib_params.Ksat     = new double[state->lgar_bmi_params.num_layers];
  
  // initialize thickness/depth and soil moisture of wetting fronts (used for model coupling)
  // also initialize calibratable parameters
  state->lgar_calib_params.field_capacity_psi = state->lgar_bmi_params.field_capacity_psi_cm;
  state->lgar_calib_params.ponded_depth_max = state->lgar_bmi_params.ponded_depth_max_cm;

  struct wetting_front *current = state->head;
  for (int i=0; i<state->lgar_bmi_params.num_wetting_fronts; i++) {
    assert (current != NULL);
    
    soil = state->lgar_bmi_params.layer_soil_type[i+1];

    state->lgar_bmi_params.soil_moisture_wetting_fronts[i] = current->theta;
    state->lgar_bmi_params.soil_depth_wetting_fronts[i]    = current->depth_cm * state->units.cm_to_m;

    state->lgar_calib_params.theta_e[i]  = state->soil_properties[soil].theta_e;
    state->lgar_calib_params.theta_r[i]  = state->soil_properties[soil].theta_r;
    state->lgar_calib_params.vg_n[i]     = state->soil_properties[soil].vg_n;
    state->lgar_calib_params.vg_alpha[i] = state->soil_properties[soil].vg_alpha_per_cm;
    state->lgar_calib_params.Ksat[i]     = state->soil_properties[soil].Ksat_cm_per_h;
    
    current = current->next;
  }


  /* initialize bmi input variables to -1.0 (on purpose), this should be assigned (non-negative) and if not,
     the code will throw an error in the Update method */
  state->lgar_bmi_input_params->precipitation_mm_per_h = -1.0;
  state->lgar_bmi_input_params->PET_mm_per_h           = -1.0;

  // initialize all global mass balance variables to zero
  state->lgar_mass_balance.volprecip_cm              = 0.0;
  state->lgar_mass_balance.volin_cm                  = 0.0;
  state->lgar_mass_balance.volend_cm                 = 0.0;
  state->lgar_mass_balance.volAET_cm                 = 0.0;
  state->lgar_mass_balance.volrech_cm                = 0.0;
  state->lgar_mass_balance.volrunoff_cm              = 0.0;
  state->lgar_mass_balance.volrunoff_giuh_cm         = 0.0;
  state->lgar_mass_balance.volQ_cm                   = 0.0;
  state->lgar_mass_balance.volPET_cm                 = 0.0;
  state->lgar_mass_balance.volon_cm                  = 0.0;
  state->lgar_mass_balance.volon_timestep_cm         = 0.0; /* setting volon and precip at the initial time to 0.0
							       as they determine the creation of surficail wetting front */
  state->lgar_bmi_params.precip_previous_timestep_cm = 0.0;
  state->lgar_mass_balance.volQ_gw_timestep_cm       = 0.0; /* setting flux from groundwater_reservoir_to_stream to zero,
							       will be non-zero when groundwater reservoir is added/simulated */
  state->lgar_mass_balance.volchange_calib_cm        = 0.0;
}


// #############################################################################################################################
/*
  Read and initialize values from a configuration file
  @param verbosity              : supress all outputs, file writing if 'none', screen output and write file if 'high'
  @param layer_thickness_cm     : 1D (double) array of layer thicknesses in cm, read from config file
  @param layer_soil_type        : 1D (int) array of layers soil type, read from config file, each integer represent a soil type
  @param num_layers             : number of actual soil layers
  @param num_wetting_fronts     : number of wetting fronts
  @param num_cells_temp         : number of cells of the discretized soil temperature profile
  @param cum_layer_thickness_cm : 1D (double) array of cumulative thickness of layers, allocate memory at run time
  @param soil_depth_cm          : depth of the computational domain (i.e., depth of the last/deepest soil layer from the surface)
  @param initial_psi_cm         : model initial (psi) condition
  @param timestep_h             : model timestep in hours
  @param forcing_resolution_h   : forcing resolution in hours
  @param forcing_interval       : factor equals to forcing_resolution_h/timestep_h (used to determine model subtimestep's forcings)
  @param num_soil_types         : number of soil types; must be less than or equal to MAX_NUM_SOIL_TYPES
  @param AET_cm                 : actual evapotranspiration in cm
  @param soil_temperature       : 1D (double) array of soil temperature [K]; bmi input for coupling lasam to soil freeze thaw model
  @param soil_temperature_z     : 1D (double) array of soil discretization associated with temperature profile [m];
                                  depth from the surface in meters
  @param frozen_factor          : frozen factor causing the hydraulic conductivity to decrease due to frozen soil
                                  (when coupled to soil freeze thaw model)
  @param wilting_point_psi_cm   : wilting point (the amount of water not available for plants or not accessible by plants)
  @param field_capacity_psi_cm  : field capacity, represented with a capillary head (head above which drainage is much faster)
  @param ponded_depth_cm        : amount of water on the surface not available for surface drainage (initialized to zero)
  @param ponded_depth_max cm    : maximum amount of water on the surface not available for surface drainage (default is zero)
  @param nint                   : number of trapezoids used in integrating the Geff function (set to 120)
  @param time_s                 : current time [s] (initially set to zero)
  @param sft_coupled            : model coupling flag. if true, lasam is coupled to soil freeze thaw model; default is uncoupled version
  @param giuh_ordinates         : geomorphological instantaneous unit hydrograph
  @param num_giuh_ordinates     : number of giuh ordinates
  @param TO_enabled             : if set to true, the model uses a multilayer version of the TO recharge model to compute fluxes between a shallow water table and the vadose zone 
*/

// #############################################################################################################################
extern void InitFromConfigFile(string config_file, struct model_state *state)
{

  ifstream fp; //FILE *fp = fopen(config_file.c_str(),"r");
  fp.open(config_file);
  //struct wetting_front* head = state->head;
  
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
	std::cerr<<"Verbosity is set to \' "<<verbosity<<"\' \n";
	std::cerr<<"          *****         \n";
      }

      fp.clear();
      break;
    }
  }

  // seek to beginning of input after searching for 'verbosity'
  fp.clear();
  fp.seekg(0, fp.beg);


  if (verbosity.compare("none") != 0) {
    std::cerr<<"------------- Initialization from config file ---------------------- \n";
  }

  // setting these options to false (defualt) 
  state->lgar_bmi_params.sft_coupled           = false;
  state->lgar_bmi_params.use_closed_form_G     = false;
  state->lgar_bmi_params.adaptive_timestep     = false;
  state->lgar_bmi_params.TO_enabled            = false;
  state->lgar_bmi_params.free_drainage_enabled = false;
  // setting mass balance tolerance to be large by default; this can be specified in the config file
  state->lgar_bmi_params.mbal_tol = 1.E1;
  
  bool is_layer_thickness_set       = false;
  bool is_initial_psi_set           = false;
  bool is_timestep_set              = false;
  bool is_endtime_set               = false;
  bool is_forcing_resolution_set    = false;
  bool is_layer_soil_type_set       = false;
  bool is_wilting_point_psi_cm_set  = false;
  bool is_field_capacity_psi_cm_set = false;
  bool is_soil_params_file_set      = false;
  bool is_max_valid_soil_types_set  = false;
  bool is_giuh_ordinates_set        = false;
  bool is_soil_z_set                = false;
  bool is_ponded_depth_max_cm_set   = false;
  bool is_root_zone_depth_cm_set    = false;

  string soil_params_file;

  // a temporary array to store the original (hourly based) giuh values
  std::vector<double> giuh_ordinates_temp;
 
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

      state->lgar_bmi_params.layer_thickness_cm = new double[vec.size()+1];
      state->lgar_bmi_params.cum_layer_thickness_cm = new double[vec.size()+1];

      state->lgar_bmi_params.layer_thickness_cm[0] = 0.0; // the value at index 0 is never used
      // calculate the cumulative (absolute) depth from land surface to bottom of each soil layer
      state->lgar_bmi_params.cum_layer_thickness_cm[0] = 0.0;

      for (unsigned int layer=1; layer <= vec.size(); layer++) {
      	state->lgar_bmi_params.layer_thickness_cm[layer] = vec[layer-1];
	state->lgar_bmi_params.cum_layer_thickness_cm[layer] = state->lgar_bmi_params.cum_layer_thickness_cm[layer-1] + vec[layer-1];
      }

      state->lgar_bmi_params.num_layers = vec.size();

      state->lgar_bmi_params.soil_depth_cm = state->lgar_bmi_params.cum_layer_thickness_cm[state->lgar_bmi_params.num_layers];
      is_layer_thickness_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Number of layers : "<<state->lgar_bmi_params.num_layers<<"\n";
	for (int i=1; i<=state->lgar_bmi_params.num_layers; i++)
	  std::cerr<<"Thickness, cum. depth : "<<state->lgar_bmi_params.layer_thickness_cm[i]<<" , "
		   <<state->lgar_bmi_params.cum_layer_thickness_cm[i]<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "layer_soil_type") {
      vector<double> vec = ReadVectorData(param_value);

      state->lgar_bmi_params.layer_soil_type = new int[vec.size()+1];

      for (unsigned int layer=1; layer <= vec.size(); layer++)
      	state->lgar_bmi_params.layer_soil_type[layer] = vec[layer-1];

      is_layer_soil_type_set = true;

      continue;
    }
    else if (param_key == "giuh_ordinates") {
      vector<double> vec = ReadVectorData(param_value);

      giuh_ordinates_temp.resize(vec.size()+1);
 
      for (unsigned int i=1; i <= vec.size(); i++)
	giuh_ordinates_temp[i] = vec[i-1];

      is_giuh_ordinates_set = true;

      if (verbosity.compare("high") == 0) {
	for (int i=1; i <= (int)vec.size(); i++)
	  std::cerr<<"GIUH ordinates (hourly) : "<<giuh_ordinates_temp[i]<<"\n";

	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "soil_z") {
      vector<double> vec = ReadVectorData(param_value);

      state->lgar_bmi_params.soil_temperature_z = new double[vec.size()];

      for (unsigned int i=0; i < vec.size(); i++)
      	state->lgar_bmi_params.soil_temperature_z[i] = vec[i];

      state->lgar_bmi_params.num_cells_temp = vec.size();

      is_soil_z_set = true;

      if (verbosity.compare("high") == 0) {
	for (int i=0; i<state->lgar_bmi_params.num_cells_temp; i++)
	  std::cerr<<"Soil z (temperature resolution) : "<<state->lgar_bmi_params.soil_temperature_z[i]<<"\n";

	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "initial_psi") {
      state->lgar_bmi_params.initial_psi_cm = stod(param_value);
      is_initial_psi_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Initial Psi : "<<state->lgar_bmi_params.initial_psi_cm<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "max_valid_soil_types") {
      state->lgar_bmi_params.num_soil_types = std::min(stoi(param_value), MAX_NUM_SOIL_TYPES);
      is_max_valid_soil_types_set = true;
      continue;
    }
    else if (param_key == "soil_params_file") {
      soil_params_file = param_value;
      is_soil_params_file_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Soil paramaters file : "<<soil_params_file<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "wilting_point_psi") {
      state->lgar_bmi_params.wilting_point_psi_cm = stod(param_value);
      is_wilting_point_psi_cm_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Wilting point Psi [cm] : "<<state->lgar_bmi_params.wilting_point_psi_cm<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "field_capacity_psi") {
      state->lgar_bmi_params.field_capacity_psi_cm = stod(param_value);
      is_field_capacity_psi_cm_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Field capacity Psi [cm] : "<<state->lgar_bmi_params.field_capacity_psi_cm<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "root_zone_depth") {
      state->lgar_bmi_params.root_zone_depth_cm = stod(param_value);

      if (state->lgar_bmi_params.root_zone_depth_cm<0.0){
        printf("root_zone_depth is less than 0 \n");
        abort();
      }

      is_root_zone_depth_cm_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"root zone depth [cm] : "<<state->lgar_bmi_params.root_zone_depth_cm<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "use_closed_form_G") { 
      if (param_value == "false") {
        state->lgar_bmi_params.use_closed_form_G = false;
      }
      else if (param_value == "true") {
        state->lgar_bmi_params.use_closed_form_G = true;
      }
      else {
	std::cerr<<"Invalid option: use_closed_form_G must be true or false, or left unspecified (defaulting to false). \n";
        abort();
      }

      continue;
    }
    else if (param_key == "adaptive_timestep") { 
      if ((param_value == "false") || (param_value == "0")) {
        state->lgar_bmi_params.adaptive_timestep = false;
      }
      else if ( (param_value == "true") || (param_value == "1")) {
        state->lgar_bmi_params.adaptive_timestep = true;
      }
      else {
	std::cerr<<"Invalid option: adaptive_timestep must be true or false, or left unspecified (defaulting to false). \n";
        abort();
      }

      continue;
    }
    else if (param_key == "free_drainage_enabled") { 
      if (param_value == "false") {
        state->lgar_bmi_params.free_drainage_enabled = false;
      }
      else if (param_value == "true") {
        state->lgar_bmi_params.free_drainage_enabled = true;
      }
      else {
	std::cerr<<"Invalid option: free_drainage_enabled must be true or false, or left unspecified (defaulting to false). \n";
        abort();
      }

      continue;
    }
    else if (param_key == "mbal_tol") {
      state->lgar_bmi_params.mbal_tol = stod(param_value);

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Mass balance tolerance [cm] : "<<state->lgar_bmi_params.mbal_tol<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "TO_enabled") { 
      if (param_value == "false") {
        state->lgar_bmi_params.TO_enabled = false;
      }
      else if (param_value == "true") {
        state->lgar_bmi_params.TO_enabled = true;
      }
      else {
	std::cerr<<"Invalid option: TO_enabled must be true or false, or left unspecified (defaulting to false). \n";
        abort();
      }

      continue;
    }
    else if (param_key == "timestep") {
      state->lgar_bmi_params.timestep_h = stod(param_value);

      if (param_unit == "[s]" || param_unit == "[sec]" || param_unit == "") // defalut time unit is seconds
	state->lgar_bmi_params.timestep_h /= 3600; // convert to hours
      else if (param_unit == "[min]" || param_unit == "[minute]")
	state->lgar_bmi_params.timestep_h /= 60; // convert to hours
      else if (param_unit == "[h]" || param_unit == "[hr]")
	state->lgar_bmi_params.timestep_h /= 1.0; // convert to hours

  state->lgar_bmi_params.minimum_timestep_h = state->lgar_bmi_params.timestep_h;

      assert (state->lgar_bmi_params.timestep_h > 0);
      is_timestep_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Model timestep [hours,seconds]: "<<state->lgar_bmi_params.timestep_h<<" , "
		 <<state->lgar_bmi_params.timestep_h*3600<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "endtime") {

      if (param_unit == "[s]" || param_unit == "[sec]" || param_unit == "") // defalut time unit is seconds
	state->lgar_bmi_params.endtime_s = stod(param_value);
      else if (param_unit == "[min]" || param_unit == "[minute]")
	state->lgar_bmi_params.endtime_s = stod(param_value) * 60.0;
      else if (param_unit == "[h]" || param_unit == "[hr]")
	state->lgar_bmi_params.endtime_s = stod(param_value) * 3600.0;
      else if (param_unit == "[d]" || param_unit == "[day]")
	state->lgar_bmi_params.endtime_s = stod(param_value) * 86400.0;

      assert (state->lgar_bmi_params.endtime_s > 0);
      is_endtime_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Endtime [days, hours]: "<< state->lgar_bmi_params.endtime_s/86400.0 <<" , "
		 << state->lgar_bmi_params.endtime_s/3600.0<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "forcing_resolution") {
      state->lgar_bmi_params.forcing_resolution_h = stod(param_value);

      if (param_unit == "[s]" || param_unit == "[sec]" || param_unit == "") // defalut time unit is seconds
	state->lgar_bmi_params.forcing_resolution_h /= 3600;                // convert to hours
      else if (param_unit == "[min]" || param_unit == "[minute]")
	state->lgar_bmi_params.forcing_resolution_h /= 60;                 // convert to hours
      else if (param_unit == "[h]" || param_unit == "[hr]")
	state->lgar_bmi_params.forcing_resolution_h /= 1.0;               // convert to hours

      assert (state->lgar_bmi_params.forcing_resolution_h > 0);
      is_forcing_resolution_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Forcing resolution [hours]: "<<state->lgar_bmi_params.forcing_resolution_h<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "sft_coupled") {
      if (param_value == "true") {
	state->lgar_bmi_params.sft_coupled = 1;
      }
      else if (param_value == "false") {
	state->lgar_bmi_params.sft_coupled = 0; // false
      }
      else {
	std::cerr<<"Invalid option: sft_coupled must be true or false. \n";
        abort();
      }
      
      continue;
    }
    else if (param_key == "ponded_depth_max") {
      state->lgar_bmi_params.ponded_depth_max_cm = fmax(stod(param_value), 0.0);
      is_ponded_depth_max_cm_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Maximum ponded depth [cm] : "<<state->lgar_bmi_params.ponded_depth_max_cm<<"\n";
	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "calib_params") {
      if (param_value == "true") {
	state->lgar_bmi_params.calib_params_flag = 1;
      }
      else if (param_value == "false") {
	state->lgar_bmi_params.calib_params_flag = 0; // false
      }
      else {
	std::cerr<<"Invalid option: calib_params must be true or false. \n";
        abort();
      }
      
      continue;
    }
  }

  fp.close();

  if (verbosity.compare("high") == 0) {
    std::string flag = state->lgar_bmi_params.use_closed_form_G == true ? "Yes" : "No";
    std::cerr<<"Using closed_form_G? "<< flag <<"\n";
    std::cerr<<"          *****         \n";
  }

  if (verbosity.compare("high") == 0) {
    std::string flag = state->lgar_bmi_params.sft_coupled == true ? "Yes" : "No";
    std::cerr<<"Coupled to SoilFreezeThaw? "<< flag <<"\n";
    std::cerr<<"          *****         \n";
  }
  
  if(!is_max_valid_soil_types_set)
     state->lgar_bmi_params.num_soil_types = MAX_NUM_SOIL_TYPES;     // maximum number of valid soil types defaults to 15

  if (verbosity.compare("high") == 0) {
    std::cerr<<"Maximum number of soil types: "<<state->lgar_bmi_params.num_soil_types<<"\n";
    std::cerr<<"          *****         \n";
  }

  if (!is_layer_soil_type_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set layer_soil_type. \n";
    throw runtime_error(errMsg.str());
  }
    
  if(is_soil_params_file_set) {
    //allocate memory to create an array of structures to hold the soils properties data.
    //state->soil_properties = (struct soil_properties_*) malloc((state->lgar_bmi_params.num_layers+1)*sizeof(struct soil_properties_));


    state->soil_properties = new soil_properties_[state->lgar_bmi_params.num_soil_types+1];
    int num_soil_types = state->lgar_bmi_params.num_soil_types;
    double wilting_point_psi_cm = state->lgar_bmi_params.wilting_point_psi_cm;
    // int max_num_soil_in_file = lgar_read_vG_param_file(soil_params_file.c_str(), num_soil_types, wilting_point_psi_cm, state->soil_properties);
    lgar_read_vG_param_file(soil_params_file.c_str(), num_soil_types, wilting_point_psi_cm, state->soil_properties);

    // check if soil layers provided are within the range
    state->lgar_bmi_params.is_invalid_soil_type = false; // model not valid for soil types = waterbody, glacier, lava, etc.
    for (int layer=1; layer <= state->lgar_bmi_params.num_layers; layer++) {
      //assert (state->lgar_bmi_params.layer_soil_type[layer] <= state->lgar_bmi_params.num_soil_types);
      //assert (state->lgar_bmi_params.layer_soil_type[layer] <= max_num_soil_in_file);
      if (state->lgar_bmi_params.layer_soil_type[layer] > state->lgar_bmi_params.num_soil_types) {
        state->lgar_bmi_params.is_invalid_soil_type = true;
        if (verbosity.compare("high") == 0) {
          std::cerr << "Invalid soil type: "
              << state->lgar_bmi_params.layer_soil_type[layer]
              <<". Model returns input_precip = ouput_Qout. \n";
        }
        break;
      }
    }

    if (verbosity.compare("high") == 0) {
      for (int layer=1; layer<=state->lgar_bmi_params.num_layers; layer++) {
	int soil = state->lgar_bmi_params.layer_soil_type[layer];
	std::cerr<<"Soil type/name : "<<state->lgar_bmi_params.layer_soil_type[layer]
		 <<" "<<state->soil_properties[soil].soil_name<<"\n";
      }
      std::cerr<<"          *****         \n";
    }
  }
  else {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set soil_params_file. \n";
    throw runtime_error(errMsg.str());
  }
  
  if (!is_layer_thickness_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set layer_thickness. \n";
    throw runtime_error(errMsg.str());
  }

  if (!is_initial_psi_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set initial_psi. \n";
    throw runtime_error(errMsg.str());
  }

  if (!is_timestep_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set timestep. \n";
    throw runtime_error(errMsg.str());
  }

  if (!is_endtime_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set endtime. \n";
    throw runtime_error(errMsg.str());
  }

  if(!is_wilting_point_psi_cm_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set wilting_point_psi. \n Recommended value of 15495.0[cm], corresponding to 15 atm. \n";
    throw runtime_error(errMsg.str());
  }

  if(!is_field_capacity_psi_cm_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set field_capacity_psi. \n Recommended value of 340.9[cm] for most soils, corresponding to 1/3 atm, or 103.3[cm] for sands, corresponding to 1/10 atm. \n";
    throw runtime_error(errMsg.str());
  }

  if(!is_root_zone_depth_cm_set && state->lgar_bmi_params.TO_enabled==true) {
    stringstream errMsg;
    errMsg << "root zone depth not set in the config file while TO mode is enabled "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }

  if (!is_forcing_resolution_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set forcing_resolution. \n";
    throw runtime_error(errMsg.str());
  }

  if (is_giuh_ordinates_set) {
    int factor = int(1.0/state->lgar_bmi_params.forcing_resolution_h);

    state->lgar_bmi_params.num_giuh_ordinates = factor * (giuh_ordinates_temp.size() - 1);
    state->lgar_bmi_params.giuh_ordinates = new double[state->lgar_bmi_params.num_giuh_ordinates+1];
    
    for (int i=0; i<int(giuh_ordinates_temp.size()-1); i++) {
      for (int j=0; j<factor; j++) {
        int index = j + i * factor + 1;
        state->lgar_bmi_params.giuh_ordinates[index] = giuh_ordinates_temp[i+1]/double(factor);
      }
    }
    
    if (verbosity.compare("high") == 0) {
      for (int i=1; i<=state->lgar_bmi_params.num_giuh_ordinates; i++)
	      std::cerr<<"GIUH ordinates (scaled) : "<<state->lgar_bmi_params.giuh_ordinates[i]<<"\n";
      
      std::cerr<<"          *****         \n";
    }
    giuh_ordinates_temp.clear();
  }

  else if (!is_giuh_ordinates_set) {
    stringstream errMsg;
    errMsg << "The configuration file \'" << config_file <<"\' does not set giuh_ordinates. \n";
    throw runtime_error(errMsg.str());
  }

  if (state->lgar_bmi_params.sft_coupled) {
    state->lgar_bmi_params.soil_temperature = new double[state->lgar_bmi_params.num_cells_temp]();
    if (!is_soil_z_set) {
      stringstream errMsg;
      errMsg << "The configuration file \'" << config_file <<"\' does not set soil_z. \n";
      throw runtime_error(errMsg.str());
    }
  }
  else {
    //NJF FIXME these arrays should be allocated based on num_cells_temp...
    state->lgar_bmi_params.soil_temperature   = new double[1]();
    state->lgar_bmi_params.soil_temperature_z = new double[1]();
    state->lgar_bmi_params.num_cells_temp     = 1;
  }

  if (!is_ponded_depth_max_cm_set)
    state->lgar_bmi_params.ponded_depth_max_cm = 0.0; // default maximum ponded depth is set to zero (i.e. no surface ponding)


  state->lgar_bmi_params.forcing_interval = int(state->lgar_bmi_params.forcing_resolution_h/state->lgar_bmi_params.timestep_h+1.0e-08); // add 1.0e-08 to prevent truncation error

  // initialize frozen factor array to 1.
  state->lgar_bmi_params.frozen_factor = new double[state->lgar_bmi_params.num_layers+1];
  for (int i=0; i <= state->lgar_bmi_params.num_layers; i++)
    state->lgar_bmi_params.frozen_factor[i] = 1.0;

  InitializeWettingFronts(state->lgar_bmi_params.TO_enabled, state->lgar_bmi_params.num_layers, state->lgar_bmi_params.initial_psi_cm,
			  state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.cum_layer_thickness_cm, state->lgar_bmi_params.layer_thickness_cm, 
			  state->lgar_bmi_params.frozen_factor, &state->head, state->soil_properties);
  
  if (verbosity.compare("none") != 0) {
    std::cerr<<"--- Initial state/conditions --- \n";
    listPrint(state->head);
    std::cerr<<"          *****         \n";
  }

  // initial mass in the system
  state->lgar_mass_balance.volstart_cm      = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);

  state->lgar_bmi_params.ponded_depth_cm    = 0.0; // initially we start with a dry surface (no surface ponding)
  state->lgar_bmi_params.nint               = 120; // hacked, not needed to be an input option
  state->lgar_bmi_params.num_wetting_fronts = state->lgar_bmi_params.num_layers;

  if (!state->lgar_bmi_params.TO_enabled){
    assert (state->lgar_bmi_params.num_layers == listLength(state->head)); //only necessary for non TO mode; number of initial wetting fronts isn't necessarily equal to the number of layers in LGARTO
  }

  if (verbosity.compare("high") == 0) {
    std::cerr<<"Initial ponded depth is set to zero. \n";
  }

  state->lgar_bmi_input_params     = new lgar_bmi_input_parameters;
  state->lgar_bmi_params.time_s    = 0.0;
  state->lgar_bmi_params.timesteps = 0.0;

  if (verbosity.compare("none") != 0) {
    std::cerr<<"------------- Initialization done! ---------------------- \n";
    std::cerr<<"--------------------------------------------------------- \n";
  }

}

//##############################################################################
/*
  calculates initial theta (soil moisture content) and hydraulic conductivity
  from the prescribed psi value for each of the soil layers
*/
// #############################################################################
extern void InitializeWettingFronts(bool TO_enabled, int num_layers, double initial_psi_cm, int *layer_soil_type, double *cum_layer_thickness_cm,
				    double *layer_thickness_cm, double *frozen_factor, struct wetting_front** head, struct soil_properties_ *soil_properties)
{

  if (TO_enabled==true){ //PTL: LGARTO case, where wetting front psi values are not initialized with inital_psi_cm but instead are based on WF depth (because they are initially hydrostatic)

    int soil;
    int layer=1;
    double Se, theta_init;
    bool bottom_flag;
    double Ksat_cm_per_h;
    struct wetting_front *current;
    int number_of_WFs_per_layer = 4;//4 //used to determine the number of initial WFs per soil layer in LGARTO. For longer simulations, both simulation quality and runtime do not depend heavilty on this, for shorter simulations it can be important. 
    bool switch_to_next_layer_flag = false;
    double prior_psi_cm = cum_layer_thickness_cm[num_layers];
    double new_wf_depth;
    int wf_in_layer = 1;
    //extra_moisture_factor and extra_height_factor, when not equal to 0, initalize the model with wetting fronts that are not at their hydrostatic positions. Using 0 for both is preferred unless you want initial conditions to deviate from hydrostatic. 
    //extra_moisture_factor and extra_height_factor have been useful for testing and development but are note recommended to be changed anymore. 
    double extra_moisture_factor = 0.0;  //20
    double extra_height_factor = 0.0;

    for(int front=1;front<=(num_layers*number_of_WFs_per_layer);front++) {

      soil = layer_soil_type[layer];
      double total_depth = cum_layer_thickness_cm[num_layers];

      if ( (front%number_of_WFs_per_layer)==0 ){ //This is if a WF is at the bottom of its layer, in which case it has to have the same psi value directly below it, unless it is the deepest WF. 
        wf_in_layer = 1; //reset this counter
        initial_psi_cm = (layer_thickness_cm[layer]/number_of_WFs_per_layer) - extra_moisture_factor;
        if (layer<num_layers){
          initial_psi_cm = (total_depth - (cum_layer_thickness_cm[layer-1] + (number_of_WFs_per_layer-1)*layer_thickness_cm[layer+1*0]/number_of_WFs_per_layer)) - extra_moisture_factor;
        }
        new_wf_depth = cum_layer_thickness_cm[layer-1] + layer_thickness_cm[layer];
        prior_psi_cm = initial_psi_cm;
      }
      else{ //This is if a WF is not at the bottom of its layer, in which case its depth will be set based on having an equal spacing of WFs in the layer, and the psi will be set based on what it should be for the hydrostatic position. Note that extra_moisture_factor and extra_height_factor move wetting fronts away from their hydrostatic positions for testing purposes.
        if ( wf_in_layer == 1 ){
          initial_psi_cm = prior_psi_cm;
          }
        else{
          initial_psi_cm = (total_depth - (cum_layer_thickness_cm[layer-1] + (wf_in_layer-1)*layer_thickness_cm[layer]/number_of_WFs_per_layer)) - extra_moisture_factor;
        } 
        new_wf_depth = (cum_layer_thickness_cm[layer-1] + wf_in_layer*layer_thickness_cm[layer]/number_of_WFs_per_layer) - extra_height_factor;
      }

      if (initial_psi_cm<0){
        initial_psi_cm = 0;
      }
      theta_init = calc_theta_from_h(initial_psi_cm,soil_properties[soil].vg_alpha_per_cm,
            soil_properties[soil].vg_m,soil_properties[soil].vg_n,
            soil_properties[soil].theta_e,soil_properties[soil].theta_r);

      if (verbosity.compare("high") == 0) {
        printf("layer, theta, psi, alpha, m, n, theta_e, theta_r = %d, %6.6f, %6.6f, %6.6f, %6.6f, %6.6f, %6.6f, %6.6f \n", layer, theta_init, initial_psi_cm, soil_properties[soil].vg_alpha_per_cm, soil_properties[soil].vg_m,soil_properties[soil].vg_n,soil_properties[soil].theta_e,soil_properties[soil].theta_r);
      }

      if ( (front%number_of_WFs_per_layer)==0 ){
        bottom_flag=true;  
      }
      else{
        bottom_flag=false; //in LGARTO, not all initial WFs extend to the bottom of their layers.
      }

      if (new_wf_depth<cum_layer_thickness_cm[layer-1]){
        new_wf_depth = cum_layer_thickness_cm[layer-1];
      // This code makes it so that if a wetting front would be higher than it is allowed to be for its layer, via accidentally making extra_height_factor too large, then the wetting front is at the top of its layer. To be safe, generally should make sure that extra_height_factor is small enough such that WFs are in their layers without this code.
      } 

      // the next lines create the initial moisture profile
      current = listInsertFront(new_wf_depth,theta_init,front,layer,bottom_flag,head,TO_enabled);

      current->psi_cm = initial_psi_cm;
      Se = calc_Se_from_theta(current->theta,soil_properties[soil].theta_e,soil_properties[soil].theta_r);

      Ksat_cm_per_h = frozen_factor[layer] * soil_properties[soil].Ksat_cm_per_h;
      current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h , soil_properties[soil].vg_m);  // cm/s

      if (switch_to_next_layer_flag){
        layer++;
        switch_to_next_layer_flag = false;
      }
      else{
        wf_in_layer ++;
        if ( ((front+1)%number_of_WFs_per_layer)==0 ){
          switch_to_next_layer_flag = true;
        }
      }
    }
  }

  else{ //case where TO_enabled is false, so all initial wetting fronts will have the same psi value 
    int soil;
    int front = 0;
    double Se, theta_init;
    bool bottom_flag;
    double Ksat_cm_per_h;
    struct wetting_front *current;

    for(int layer=1;layer<=num_layers;layer++) {
      front++;

      soil = layer_soil_type[layer];
      theta_init = calc_theta_from_h(initial_psi_cm,soil_properties[soil].vg_alpha_per_cm,
            soil_properties[soil].vg_m,soil_properties[soil].vg_n,
            soil_properties[soil].theta_e,soil_properties[soil].theta_r);

      if (verbosity.compare("high") == 0) {
        printf("layer, theta, psi, alpha, m, n, theta_e, theta_r = %d, %6.6f, %6.6f, %6.6f, %6.6f, %6.6f, %6.6f, %6.6f \n",
        layer, theta_init, initial_psi_cm, soil_properties[soil].vg_alpha_per_cm, soil_properties[soil].vg_m,
        soil_properties[soil].vg_n,soil_properties[soil].theta_e,soil_properties[soil].theta_r);
      }

      // the next lines create the initial moisture profile
      bottom_flag = true;  // all initial wetting fronts are in contact with the bottom of the layer they exist in
      // NOTE: The listInsertFront function does lots of stuff.

      current = listInsertFront(cum_layer_thickness_cm[layer],theta_init,front,layer,bottom_flag, head, TO_enabled);

      current->psi_cm = initial_psi_cm;
      Se = calc_Se_from_theta(current->theta,soil_properties[soil].theta_e,soil_properties[soil].theta_r);

      Ksat_cm_per_h = frozen_factor[layer] * soil_properties[soil].Ksat_cm_per_h;
      current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h , soil_properties[soil].vg_m);  // cm/s
    }
  }
}

// ##################################################################################
/*
  Reads 1D data from the config file
  - used for reading soil discretization (1D)
  - used for reading layers depth from the surface 
*/
// ##################################################################################
extern vector<double>
ReadVectorData(string key)
{
  int pos = 0;
  string delimiter = ",";
  vector<double> value(0.0);
  string z1 = key;

  if (z1.find(delimiter) == string::npos) {
    double v = stod(z1);
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

// ############################################################################################
/*
  calculates frozen factor based on L. Wang et al. (www.hydrol-earth-syst-sci.net/14/557/2010/)
  uses layered-average soil temperatures and an exponential function to compute frozen fraction
  for each layer
 */
// ############################################################################################
extern void frozen_factor_hydraulic_conductivity(struct lgar_bmi_parameters lgar_bmi_params)
{

  int c = 0, count;
  double layer_temp;
  double factor;

  for (int layer=1; layer<=lgar_bmi_params.num_layers; layer++) {
    layer_temp = 0.0;
    count = 0;

    while (lgar_bmi_params.soil_temperature_z[c] <= lgar_bmi_params.cum_layer_thickness_cm[layer]) {
      layer_temp +=lgar_bmi_params.soil_temperature[c];
      c++;
      count++;

      if (c == lgar_bmi_params.num_cells_temp) // to ensure we don't access out of bound values
	break;
    }

    // assert (layer_temp > 100.0); // just a check to ensure the while loop executes at least once
    assert (count > 0); // just a check to ensure the while loop executes at least once

    layer_temp /= count;  // layer-averaged temperature

    factor = exp(-10 * (273.15 - layer_temp)); /* Eq. 6 (L. Wang et al.,Frozen soil parameterization in a distributed
                                                biosphere hydrological model, www.hydrol-earth-syst-sci.net/14/557/2010/)
						and Eq. 22 (Bao et al., An enthalpy-based frozen model) */

    factor = fmax(fmin(factor,1.0), 0.05); // 0.05 <= factor <= 1.0
    lgar_bmi_params.frozen_factor[layer] = factor;

  }

  if (verbosity.compare("high") == 0) {
    for (int i=1; i <= lgar_bmi_params.num_layers; i++)
      std::cerr<<"frozen factor = "<< lgar_bmi_params.frozen_factor[i]<<"\n";
  }

}

/*
extern void lgar_update(struct model_state *state)
{ if we ever decided to run this version without the bmi then we simply need to copy `update method` from the bmi here.}
*/


// #########################################################################################
/*
  calculates global mass balance at the end of simulation
*/
// #########################################################################################
extern void lgar_global_mass_balance(struct model_state *state, double *giuh_runoff_queue_cm)
{
  double volstart           = state->lgar_mass_balance.volstart_cm;
  double volprecip          = state->lgar_mass_balance.volprecip_cm;
  double volrunoff          = state->lgar_mass_balance.volrunoff_cm;
  double volAET             = state->lgar_mass_balance.volAET_cm;
  double volPET             = state->lgar_mass_balance.volPET_cm;
  double volon              = state->lgar_mass_balance.volon_cm;
  double volin              = state->lgar_mass_balance.volin_cm;
  double volrech            = state->lgar_mass_balance.volrech_cm;
  double volend             = state->lgar_mass_balance.volend_cm;
  double volrunoff_giuh     = state->lgar_mass_balance.volrunoff_giuh_cm;
  double volend_giuh_cm     = 0.0;
  double total_Q_cm         = state->lgar_mass_balance.volQ_cm;
  double volchange_calib_cm = state->lgar_mass_balance.volchange_calib_cm;
  
  //check if the giuh queue have some water left at the end of simulaiton; needs to be included in the global mass balance
  // hold on; this is probably not needed as we have volrunoff in the balance; revist AJK
  for(int i=0; i <= state->lgar_bmi_params.num_giuh_ordinates; i++)
    volend_giuh_cm += giuh_runoff_queue_cm[i];


  double global_error_cm = volstart + volprecip - volrunoff - volAET - volon - volrech - volend + volchange_calib_cm;

  // if (volAET<0.0 || volrech<-10.0){
  //   printf("AET less than 0, or recharge is unrealistically small, likely because free drainage is not handled properly \n");
  //   abort();
  // }

  // else{
    printf("\n********************************************************* \n");
    printf("-------------------- Simulation Summary ----------------- \n");
    //printf("Time (sec)                 = %6.10f \n", elapsed);
    printf("------------------------ Mass balance ------------------- \n");
    printf("Initial water in soil     = %14.10f cm\n", volstart);
    printf("Total precipitation       = %14.10f cm\n", volprecip);
    printf("Total infiltration        = %14.10f cm\n", volin);
    printf("Final water in soil       = %14.10f cm\n", volend);
    printf("Surface ponded water      = %14.10f cm\n", volon);
    printf("Surface runoff            = %14.10f cm\n", volrunoff);
    printf("GIUH runoff               = %14.10f cm\n", volrunoff_giuh);
    printf("GIUH water (in array)     = %14.10f cm\n", volend_giuh_cm);
    printf("Total percolation         = %14.10f cm\n", volrech);
    printf("Total AET                 = %14.10f cm\n", volAET);
    printf("Total PET                 = %14.10f cm\n", volPET);
    printf("Total discharge (Q)       = %14.10f cm\n", total_Q_cm);
    printf("Vol change (calibration)  = %14.10f cm\n", volchange_calib_cm);
    printf("Global balance            =   %.6e cm\n", global_error_cm);
  // }

}

// ############################################################################################
/*
 finds the wetting front that corresponds to psi (head) value closest to zero
 (i.e., saturation in terms of psi). This is the wetting front that experiences infiltration
 and actual ET based on precipitation or ponded head and PET, respectively. For example, the actual ET
 is extracted from this wetting front plus the wetting fronts above it.
 Note: the free_drainage name came from its python version. In early development of LGAR in Python, 
 we initially included free drainage, and free drainage demand was extracted from this wetting front as well.
 */
// ############################################################################################
extern int wetting_front_free_drainage(struct wetting_front* head) {

  int wf_that_supplies_free_drainage_demand = 1;
  struct wetting_front *current;

  //current = head;
  int number_of_wetting_fronts = listLength(head);

  for(current = head; current != NULL; current = current->next)
  {
    if (current->next != NULL) {
      if ((current->layer_num == current->next->layer_num) && (current->is_WF_GW==0))
	break;
      else
	wf_that_supplies_free_drainage_demand++;

    }
  }

  if (wf_that_supplies_free_drainage_demand > number_of_wetting_fronts)
    wf_that_supplies_free_drainage_demand--;

  if (verbosity.compare("high") == 0) {
    printf("wetting_front_free_drainage = %d \n", wf_that_supplies_free_drainage_demand);
  }

  return  wf_that_supplies_free_drainage_demand;
}

// #######################################################################################################
/*
  the function moves wetting fronts, merge wetting fronts and does the mass balance correction when needed
  @param current : wetting front pointing to the current node of the current state
  @param next    : wetting front pointing to the next node of the current state
  @param previous    : wetting front pointing to the previous node of the current state
  @param current_old : wetting front pointing to the current node of the previous state
  @param next_old    : wetting front pointing to the next node of the previous state
  @param head : pointer to the first wetting front in the list of the current state

  Note: '_old' denotes the wetting_front or variables at the previous timestep (or state)
*/
// #######################################################################################################
extern double lgar_move_wetting_fronts(bool TO_enabled, double timestep_h, double *free_drainage_subtimestep_cm, double PET_timestep_cm, double wilting_point_psi_cm, double field_capacity_psi_cm, double root_zone_depth_cm,
             double *volin_cm, int wf_free_drainage_demand, double old_mass, int num_layers, double surf_frac_rz, double *AET_demand_cm, double *cum_layer_thickness_cm,
				     int *soil_type, double *frozen_factor, struct wetting_front** head,
				     struct wetting_front* state_previous, struct soil_properties_ *soil_properties, double *surf_AET_vec)
{

  if (verbosity.compare("high") == 0) {
    printf("State before moving wetting fronts...\n");
    listPrint(*head);
  }

  struct wetting_front *current;
  struct wetting_front *next;
  struct wetting_front *previous;

  struct wetting_front *current_old;
  struct wetting_front *next_old;

  double column_depth = cum_layer_thickness_cm[num_layers];

  previous = *head;
  double theta_e,theta_r;
  double vg_a, vg_m, vg_n;
  int layer_num, soil_num;

  int number_of_wetting_fronts = listLength(*head); //when running in LGAR mode, all WFs are surface WFs. When running in LGARTO mode, not all WFs are surface WFs -- some are groundwater WFs.
  //since there are different rules for advancing GW / surface WFs, and since these result in different blocks of code, we need to know how many WFs are surface WFs, if TO mode is on.
  int number_of_surface_WFs = 0;
  int number_of_TO_WFs_above_surface_WFs = 0;
  if (TO_enabled==true){
    for (int swf=1; swf != listLength(*head); swf++) {
      current = listFindFront(swf, *head, NULL);
      if (!(current->is_WF_GW)){
        number_of_surface_WFs ++;
      }
    }
    for (int TOwf=1; TOwf != listLength(*head); TOwf++) {
      current = listFindFront(TOwf, *head, NULL);
      if (current->is_WF_GW){
        number_of_TO_WFs_above_surface_WFs ++;
      }
      if (current->next->is_WF_GW==0){
        break;
      }
      if ( (TOwf == 1) & (current->is_WF_GW==0) ) {
        break;
      }
    }
  }else
  {
    number_of_surface_WFs = number_of_wetting_fronts;
  }

  if (number_of_surface_WFs == 0){
    number_of_TO_WFs_above_surface_WFs = 0;
  }

  if (TO_enabled==false){
    number_of_TO_WFs_above_surface_WFs = 0;
  }

  current = *head;

  int last_wetting_front_index = number_of_wetting_fronts;
  int layer_num_above, layer_num_below;

  double precip_mass_to_add = (*volin_cm); // water to be added to the soil
  double bottom_boundary_flux_cm_from_surf_WFs = 0.0; //water that was part of surface WFs that leaves the model lower boundary because a surf WF got too deep

  double bottom_boundary_flux_cm = 0.0; // water leaving the system through the bottom boundary
  struct wetting_front *deepest_surf_WF = listFindFront(listLength_surface(*head) + listLength_TO_WFs_above_surface_WFs(*head), *head, NULL);
  double deepest_surf_depth_at_start = 0.0;
  if (deepest_surf_WF!=NULL){
    deepest_surf_depth_at_start = deepest_surf_WF->depth_cm;
  }

  *volin_cm = 0.0; // assuming that all the water can fit in, if not then re-assign the left over water at the end

  /* ************************************************************ */
  // main loop advancing all surface wetting fronts and doing the mass balance
  // loop goes over deepest to top most surface wetting front
  // wf denotes wetting front

  for (int wf = number_of_surface_WFs+number_of_TO_WFs_above_surface_WFs; wf != 0; wf--) {
    current = listFindFront(wf, *head, NULL);

    if (current->is_WF_GW==1){// Just an extra check; this loop should only be for surface WFs
      break;//maybe continue, in case a TO WF gets between surf WFs?
    }

    if (verbosity.compare("high") == 0) {
      printf("Moving |******** Wetting Front = %d *********| \n", wf);
    }

    if (wf == 1 && number_of_wetting_fronts >0) {
      current = listFindFront(wf, *head, NULL);
      next = listFindFront(wf+1, *head, NULL);
      previous = NULL;

      current_old = listFindFront(wf, *head, state_previous);
      next_old = listFindFront(wf+1, *head, state_previous);
    }
    else if (wf < number_of_wetting_fronts) {
      current = listFindFront(wf, *head, NULL);
      next = listFindFront(wf+1, *head, NULL);
      previous = listFindFront(wf-1, *head, NULL);

      current_old = listFindFront(wf, *head, state_previous);
      next_old = listFindFront(wf+1, *head, state_previous);
    }
    else if (wf == number_of_wetting_fronts) {
      current = listFindFront(wf, *head, NULL);
      next = NULL;
      previous = listFindFront(wf-1, *head, NULL);

      current_old = listFindFront(wf, *head, state_previous);
      next_old = NULL;
    }

    layer_num   = current->layer_num;
    soil_num    = soil_type[layer_num];
    theta_e     = soil_properties[soil_num].theta_e;
    theta_r     = soil_properties[soil_num].theta_r;
    vg_a        = soil_properties[soil_num].vg_alpha_per_cm;
    vg_m        = soil_properties[soil_num].vg_m;
    vg_n        = soil_properties[soil_num].vg_n;

    // find indices of above and below layers
    layer_num_above = (wf == 1) ? layer_num : previous->layer_num;
    layer_num_below = (wf == last_wetting_front_index) ? layer_num + 1 : next->layer_num;

    if (verbosity.compare("high") == 0) {
       printf ("Layers (current, above, below) == %d %d %d \n", layer_num, layer_num_above, layer_num_below);
       listPrint(*head);
    }

    // double free_drainage_demand = 0.0;
    double free_drainage_demand = *free_drainage_subtimestep_cm;
    double actual_ET_demand;//big change 12 oct
    if (!TO_enabled){
      actual_ET_demand = *AET_demand_cm;
    }
    else{
      actual_ET_demand = surf_AET_vec[wf];
      *AET_demand_cm += surf_AET_vec[wf];
    }

    // case to check if the wetting front is at the interface, i.e. deepest wetting front within a layer
    // psi of the layer below is already known/updated, so we that psi to compute the theta of the deepest current layer
    // todo. this condition can be replace by current->to_depth = FALSE && l<last_wetting_front_index
    /*             _____________
       layer_above             |
                            ___|
			   |
                   ________|____    <----- wetting fronts at the interface have same psi value
		           |
       layer current    ___|
                       |
                   ____|________
       layer_below     |
                  _____|________
    */
    /*************************************************************************************/
    if ( (wf < last_wetting_front_index) && (layer_num_below != layer_num) ) {
      
      if (verbosity.compare("high") == 0) {
	printf("case (deepest wetting front within layer) : layer_num (%d) != layer_num_below (%d) \n", layer_num, layer_num_below);
      }
      current->theta = calc_theta_from_h(next->psi_cm, vg_a,vg_m, vg_n, theta_e, theta_r);
      current->psi_cm = next->psi_cm;
    }

    // case to check if the number of wetting fronts are equal to the number of layers, i.e., one wetting front per layer
    /*************************************************************************************/
    /* For example, 3 layers and 3 wetting fronts in a state. psi profile is constant, and theta profile is non-uniform due
       to different van Genuchten parameters
                theta profile       psi profile  (constant head)
               _____________       ______________
                         |                   |
               __________|__       __________|___
	               |                     |
               ________|____       __________|___
                   |                         |
               ____|________       __________|___
    */

    if (wf == number_of_wetting_fronts && layer_num_below != layer_num && number_of_wetting_fronts == num_layers) {

      if (verbosity.compare("high") == 0) {
	printf("case (number_of_wetting_fronts equal to num_layers) : l (%d) == num_layers (%d) == num_wetting_fronts(%d) \n", wf, num_layers,number_of_wetting_fronts);
      }

      // local variables
      double vg_a_k, vg_m_k, vg_n_k;
      double theta_e_k, theta_r_k;

      current->depth_cm += current->dzdt_cm_per_h * timestep_h; // this is probably not needed, as dz/dt = 0 for the deepest wetting front

      double *delta_thetas = (double *) malloc(sizeof(double)*(layer_num+1));
      double *delta_thickness = (double *) malloc(sizeof(double)*(layer_num+1));

      double psi_cm_old = current_old->psi_cm;

      double psi_cm = current->psi_cm;

      double prior_mass = (current_old->depth_cm - cum_layer_thickness_cm[layer_num-1]) * (current_old->theta - 0.0); // 0.0 = next_old->theta

      double new_mass = (current->depth_cm - cum_layer_thickness_cm[layer_num-1]) * (current->theta - 0.0); // 0.0 = next->theta;

      for (int k=1; k<layer_num; k++) {
	int soil_num_k  = soil_type[k];
	theta_e_k = soil_properties[soil_num_k].theta_e;
	theta_r_k = soil_properties[soil_num_k].theta_r;
	vg_a_k    = soil_properties[soil_num_k].vg_alpha_per_cm;
	vg_m_k    = soil_properties[soil_num_k].vg_m;
	vg_n_k    = soil_properties[soil_num_k].vg_n;

	// using psi_cm_old for all layers because the psi is constant across layers in this particular case
	double theta_old             = calc_theta_from_h(psi_cm_old, vg_a_k, vg_m_k, vg_n_k, theta_e_k,theta_r_k);
	double theta_below_old       = 0.0;
	double local_delta_theta_old = theta_old - theta_below_old;
	double layer_thickness       = cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k-1];

	prior_mass += (layer_thickness * local_delta_theta_old);

	double theta       = calc_theta_from_h(psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
	double theta_below = 0.0;

	new_mass += layer_thickness * (theta - theta_below);
  //NJF theta_below is always 0, so all delta_thetas are always 0...
  //does this really need a dynamic array in this case???
	delta_thetas[k] = theta_below;
	delta_thickness[k] = layer_thickness;
      }

      delta_thetas[layer_num] = 0.0;
      delta_thickness[layer_num] = current->depth_cm - cum_layer_thickness_cm[layer_num-1];

      // double free_drainage_demand = 0.0;
      double free_drainage_demand = *free_drainage_subtimestep_cm;

      if (TO_enabled){
        if (wf_free_drainage_demand == wf){ 
          prior_mass += precip_mass_to_add - free_drainage_demand;
        }

        if (wf <= (listLength_surface(*head)+listLength_TO_WFs_above_surface_WFs(*head))){
          prior_mass += -(actual_ET_demand);
        }
      }
      else{
        if (wf_free_drainage_demand == wf)
    prior_mass += precip_mass_to_add - (free_drainage_demand + actual_ET_demand);
      }

      // theta mass balance computes new theta that conserves the mass; new theta is assigned to the current wetting front
      double AET_before_AET_demand_cm = *AET_demand_cm;

      double theta_new = lgar_theta_mass_balance(layer_num, soil_num, psi_cm, new_mass, prior_mass, AET_demand_cm,
						 delta_thetas, delta_thickness, soil_type, soil_properties);

      double AET_after_AET_demand_cm = *AET_demand_cm;
      double diff_in_AET = AET_after_AET_demand_cm - AET_before_AET_demand_cm;
      actual_ET_demand += diff_in_AET;
      //done with delta_thetas and delta_thickness, cleanup memory
      free(delta_thetas);
      free(delta_thickness);


      
      current->theta = fmax(theta_r, fmin(theta_new, theta_e));

      double Se = calc_Se_from_theta(current->theta,theta_e,theta_r);
      current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);

      /* note: theta and psi of the current wetting front are updated here based on the wetting front's mass balance,
	 upper wetting fronts will be updated later in the lgar_merge_ module (the place where all state
	 variables are updated before proceeding to the next timestep. */

    }


    // case to check if the 'current' wetting front is within the layer and not at the layer's interface
    // layer_num == layer_num_below means there is another wetting front below the current wetting front
    // and they both belong to the same layer (in simple words, wetting fronts not at the interface)
    // l < last_wetting_front_index means that the current wetting front is not the deepest wetting front in the domain
    /*************************************************************************************/
  
    if ( (wf < last_wetting_front_index) && (layer_num == layer_num_below) ) {


      if (verbosity.compare("high") == 0) {
	printf("case (wetting front within a layer) : layer_num (%d) == layer_num_below (%d) \n", layer_num,layer_num_below);
      }

      // if wetting front is the most surficial wetting front
      if (layer_num == 1) {

	double free_drainage_demand = *free_drainage_subtimestep_cm;
	// prior mass = mass contained in the current old wetting front
	double prior_mass = current_old->depth_cm * (current_old->theta -  next_old->theta);

  if (TO_enabled){
    if (wf_free_drainage_demand == wf){
      prior_mass += precip_mass_to_add - free_drainage_demand;
    }

    if (wf <= (listLength_surface(*head)+listLength_TO_WFs_above_surface_WFs(*head))){
      prior_mass += -( actual_ET_demand );
    }
  }
  else{
    if (wf_free_drainage_demand == wf)
      prior_mass += precip_mass_to_add - (free_drainage_demand + actual_ET_demand);
  }

	current->depth_cm += current->dzdt_cm_per_h * timestep_h;

	/* condition to bound the wetting front depth, if depth of a wf, at this timestep,
	   gets greater than the domain depth, it will be merge anyway as it is passing
	   the layer depth */
	// current->depth_cm = fmin(current->depth_cm, column_depth); //removed because this does not work with LGARTO when a surface WF passes to a lower layer where the next TO WF has a very close theta value to the surface WF, removing does not affect LGAR performance 

	if (current->dzdt_cm_per_h == 0.0 && current->to_bottom == FALSE) // a new front was just created, so don't update it.
	  current->theta = current->theta;
	else {
      if ((prior_mass/current->depth_cm + next->theta)<theta_r){
        //the idea here is that in some cases, the reduction in theta via WF movement or AET will be intense enough such that theta goes below theta_r.
        //it requires a fairly unusual soil, or strong free drainage, which I encountered during random parameter sampling.

        double mass_before_theta_went_below_theta_r = lgar_calc_mass_bal(cum_layer_thickness_cm, *head) - current->depth_cm*(current->theta - (prior_mass/current->depth_cm + next->theta));
        if (!TO_enabled && listLength(*head)>1){//instead of deleting, we rely on how in LGARTO mode this WF will now merge with the WF below.
          current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
          current = next;
          //possible current should be made to next
        }
        double mass_after_theta_went_below_theta_r = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        *AET_demand_cm = *AET_demand_cm - fabs(mass_before_theta_went_below_theta_r - mass_after_theta_went_below_theta_r);
        actual_ET_demand -= fabs(mass_before_theta_went_below_theta_r - mass_after_theta_went_below_theta_r);
      }
      else {//This is the case where theta>theta_r, which will be almost all of the time 
	      current->theta = fmax(theta_r + 0*1.E-9, fmin(theta_e, prior_mass/current->depth_cm + next->theta));
      }
    }

      }
      else {

	/*
	  this note is copied from Python version:
	  "However, calculation of theta via mass balance is a bit trickier. This is because each wetting front
	  in deeper layers can be thought of as extending all the way to the surface, in terms of psi values.
	  For example, a wetting front in layer 2 with a theta value of 0.4 will in reality extend to layer
	  1 with a theta value that is different (usually smaller) due to different soil hydraulic properties.
	  But, the theta value of this extended wetting front is not recorded in current or previous states.
          So, simply from states, the mass balance of a wetting front that, in terms of psi, extends between
	  multiple layers cannot be calculated. Therefore, the theta values that the current wetting front *would*
	  have in above layers is calculated from the psi value of the current wetting front, with the assumption
	  that the hydraulic head of this wetting front is the same all the way up to the surface.

	  - LGAR paper (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022WR033742) has a better description, using diagrams, of the mass balance of wetting fronts
	*/

	double vg_a_k, vg_m_k, vg_n_k;
	double theta_e_k, theta_r_k;

	current->depth_cm += current->dzdt_cm_per_h * timestep_h;

	double *delta_thetas    = (double *)malloc(sizeof(double)*(layer_num+1));
	double *delta_thickness = (double *)malloc(sizeof(double)*(layer_num+1));


	double psi_cm_old = current_old->psi_cm;
	double psi_cm_below_old = current_old->next->psi_cm;

	double psi_cm = current->psi_cm;
	double psi_cm_below = next->psi_cm;

	// mass = delta(depth) * delta(theta)
	//      = difference in current and next wetting front thetas times depth of the current wetting front
	double prior_mass = (current_old->depth_cm - cum_layer_thickness_cm[layer_num-1]) * (current_old->theta - next_old->theta);
	double new_mass = (current->depth_cm - cum_layer_thickness_cm[layer_num-1]) * (current->theta - next->theta);

	// compute mass in the layers above the current wetting front
	// use the psi of the current wetting front and van Genuchten parameters of
	// the respective layers to get the total mass above the current wetting front
	for (int k=1; k<layer_num; k++) {
	  int soil_num_k  = soil_type[k];
	  theta_e_k = soil_properties[soil_num_k].theta_e;
	  theta_r_k = soil_properties[soil_num_k].theta_r;
	  vg_a_k    = soil_properties[soil_num_k].vg_alpha_per_cm;
	  vg_m_k    = soil_properties[soil_num_k].vg_m;
	  vg_n_k    = soil_properties[soil_num_k].vg_n;

	  double theta_old = calc_theta_from_h(psi_cm_old, vg_a_k, vg_m_k, vg_n_k, theta_e_k,theta_r_k);
	  double theta_below_old = calc_theta_from_h(psi_cm_below_old, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
	  double local_delta_theta_old = theta_old - theta_below_old;
	  double layer_thickness = (cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k-1]);

	  prior_mass += (layer_thickness * local_delta_theta_old);

	  //-------------------------------------------
	  // do the same for the current state
	  double theta = calc_theta_from_h(psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);

	  double theta_below = calc_theta_from_h(psi_cm_below, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);

	  new_mass += layer_thickness * (theta - theta_below);

	  delta_thetas[k] = theta_below;
	  delta_thickness[k] = layer_thickness;
	}

	delta_thetas[layer_num] = next->theta;
	delta_thickness[layer_num] = current->depth_cm - cum_layer_thickness_cm[layer_num-1];

	double free_drainage_demand = *free_drainage_subtimestep_cm;
  

  if (TO_enabled){
    if (wf_free_drainage_demand == wf){
      prior_mass += precip_mass_to_add - free_drainage_demand;
    }

    if (wf <= (listLength_surface(*head)+listLength_TO_WFs_above_surface_WFs(*head))){
      prior_mass += -( actual_ET_demand );
    }
  }
  else{
    if (wf_free_drainage_demand == wf)
      prior_mass += precip_mass_to_add - (free_drainage_demand + actual_ET_demand);
  }


  // theta mass balance computes new theta that conserves the mass; new theta is assigned to the current wetting front
  double AET_before_AET_demand_cm = *AET_demand_cm;
	double theta_new = lgar_theta_mass_balance(layer_num, soil_num, psi_cm, new_mass, prior_mass, AET_demand_cm,
						   delta_thetas, delta_thickness, soil_type, soil_properties);
  double AET_after_AET_demand_cm = *AET_demand_cm;
  double diff_in_AET = AET_after_AET_demand_cm - AET_before_AET_demand_cm;
  actual_ET_demand += diff_in_AET;
  //done with delta_thetas and delta_thickness, cleanup memory
  free(delta_thetas);
  free(delta_thickness);

	current->theta = fmax(theta_r, fmin(theta_new, theta_e));

      }

      double Se = calc_Se_from_theta(current->theta,theta_e,theta_r);
      current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);

    }
  

    // if f_p (predicted infiltration) causes theta > theta_e, mass correction is needed.
    // depth of the wetting front is increased to close the mass balance when theta > theta_e.
    // l == 1 is the last iteration (top most wetting front), so do a check on the mass balance)

    int top_most_surface_WF = 1;

    if (TO_enabled==true){
      top_most_surface_WF = listLength_TO_WFs_above_surface_WFs(*head) + 1;
    }

    if ((wf == 1) || (wf==top_most_surface_WF)) { 
      free_drainage_demand = *free_drainage_subtimestep_cm;
      wf_free_drainage_demand = wetting_front_free_drainage(*head);
      struct wetting_front *wf_free_drainage = listFindFront(wf_free_drainage_demand, *head, NULL);

      int soil_num_k1  = soil_type[wf_free_drainage->layer_num];
      double theta_e_k1 = soil_properties[soil_num_k1].theta_e;

      double mass_timestep = (old_mass + precip_mass_to_add) - (*AET_demand_cm + free_drainage_demand); //we want to replace actual_ET_demand with *AET_demand_cm for LGARTO, because in LGAR the AET will be taken from just the top most WF in terms of psi, but in LGARTO AET will be extracted from all surf WFs in the root zone 

      assert (old_mass > 0.0);
      
      if (fabs(wf_free_drainage->theta - theta_e_k1) < 1E-15) {//basically, if the free drainage WF reached saturation 
	
        double current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

        double mass_balance_error = fabs(current_mass - mass_timestep); // mass error

        double factor = 1.0;
        if (wf_free_drainage->next!=NULL){
          if (fabs(wf_free_drainage->theta - wf_free_drainage->next->theta)<0.01){// if two adjacent theta values are quite close, an initial factor of 1.0 might not make the mass balance close within 10000 iterations
            factor = factor / fabs(wf_free_drainage->theta - wf_free_drainage->next->theta);
          }
        }

        bool switched = false;
        double tolerance = 1e-10;

        // check if the difference is less than the tolerance
        if (mass_balance_error <= tolerance) {
          // return current_mass;
        }

        double depth_new = wf_free_drainage->depth_cm;

        // loop to adjust the depth for mass balance
        int iter = 0;
        bool iter_aug_flag = FALSE;
        bool break_flag = FALSE;
        while (fabs(mass_balance_error - tolerance) > 1.E-10) {
          iter++;
          if (iter>1e4) {
            break_flag = TRUE;
            bottom_boundary_flux_cm = bottom_boundary_flux_cm + mass_balance_error;
            // actual_ET_demand += mass_balance_error;

            break;
          }

          if ((iter>1e3) && (!iter_aug_flag)){
            factor = factor * 100;
            iter_aug_flag = TRUE;
          }

          if (current_mass < mass_timestep) {
            depth_new += 0.01 * factor;
            switched = false;
          }
          else {
            if (!switched) {
              switched = true;
              factor = factor * 0.001;
            }
            depth_new -= 0.01 * factor;

          }

          if ( (wf_free_drainage->to_bottom==TRUE) && (wf_free_drainage->layer_num==num_layers) ){
            depth_new = cum_layer_thickness_cm[num_layers];
          }

          wf_free_drainage->depth_cm = depth_new;

          current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
          mass_balance_error = fabs(current_mass - mass_timestep);

        }
        
        //there is a general class of problem where a very small psi value that is greater than 0 will for some but not all soils mathematically yield theta = theta_e, even though theta should be slightly less than theta_e.
        //in layered soils, this can cause a mass balance error. It is fairly rare and only seems to impact cases where the model domain is entirely saturated, which shouldn't happen when LGAR is applied in the correct environment / with sufficient layer thicknesses.
        if (break_flag) {
          current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
          mass_timestep = (old_mass + precip_mass_to_add) - (actual_ET_demand + free_drainage_demand);
          mass_balance_error = mass_timestep - current_mass;
          bottom_boundary_flux_cm += mass_balance_error;
        }

      }
    }
    
  }
  
  /*******************************************************************/
  // end of the for loop
  /*******************************************************************/


  if (verbosity.compare("high") == 0) {
    printf("State after moving but before merging wetting fronts...\n");
    listPrint(*head);
  }

  // ********************** MERGING AND CROSSING WETTING FRONT ****************************
  
  /* In the python version of LGAR, wetting front merging, layer boundary crossing, and lower boundary crossing
     all occur in a loop that happens after wetting fronts have been moved. This prevents the model from crashing,
     as there are rare but possible cases where multiple merging / layer boundary crossing events will happen in
     the same time step. For example, if two wetting fronts cross a layer boundary in the same time step, it will
     be necessary for merging to occur before layer boundary crossing. So, LGAR-C now approaches merging in the
     same way as in LGAR-Py, where wetting fronts are moved, then a new loop does merging for all wetting fronts,
     then a new loop does layer boundary corssing for all wetting fronts, then a new loop does merging again for
     all wetting fronts, and then a new loop does lower boundary crossing for all wetting fronts. this is a thorough
     way to deal with these scenarios. */

  /* Note: we check for dry over wet case before we call merge_wetting_fronts to avoid negative wetting fronts
     due to unknown corner/rare cases */

  double mass_change = 0.0;
  

  int correction_type_surf =  lgarto_correction_type_surf(num_layers, cum_layer_thickness_cm, head);

  while (correction_type_surf!=0){

    if (correction_type_surf==1){
      lgar_merge_wetting_fronts(head, soil_type, frozen_factor, soil_properties);
    }

    if (correction_type_surf==2){
      double mass_corr_loop_start = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
      bool close_psis = correct_close_psis(soil_type, soil_properties, head);
      if (close_psis){
        bottom_boundary_flux_cm += (mass_corr_loop_start - lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
      }
      lgar_wetting_fronts_cross_layer_boundary(num_layers, cum_layer_thickness_cm, soil_type, frozen_factor, head, soil_properties);
    }

    if (correction_type_surf==3){
      bottom_boundary_flux_cm += lgar_wetting_front_cross_domain_boundary(TO_enabled, cum_layer_thickness_cm, soil_type, frozen_factor, head, soil_properties);
          if (isnan(bottom_boundary_flux_cm)){
            bottom_boundary_flux_cm = 0.0;
          }
    }

    if (correction_type_surf==4){
      lgar_fix_dry_over_wet_wetting_fronts(&mass_change, cum_layer_thickness_cm, soil_type, head, soil_properties);
      if (*free_drainage_subtimestep_cm==0.0){
        *AET_demand_cm = *AET_demand_cm - mass_change;
      }
      else {
        *free_drainage_subtimestep_cm = *free_drainage_subtimestep_cm - mass_change;
        // *AET_demand_cm = *AET_demand_cm - mass_change;
        if (*free_drainage_subtimestep_cm<0.0){
          *AET_demand_cm = *AET_demand_cm + *free_drainage_subtimestep_cm;
          *free_drainage_subtimestep_cm = 0.0;
        }
      }
    }

    correction_type_surf =  lgarto_correction_type_surf(num_layers, cum_layer_thickness_cm, head);
    if (verbosity.compare("high") == 0) {
      printf("correction_type_surf at end of iteration in while loop: %d \n", correction_type_surf);
    }
  }

  if (*AET_demand_cm<0.0 && *free_drainage_subtimestep_cm>0.0){
    *free_drainage_subtimestep_cm += *AET_demand_cm;
    *AET_demand_cm = 0.0;
  }

  mass_change = 0.0;

  *volin_cm = bottom_boundary_flux_cm;

  bottom_boundary_flux_cm_from_surf_WFs = bottom_boundary_flux_cm;



  if (verbosity.compare("high") == 0) {
    printf ("mass change/adjustment (dry_over_wet case) = %lf \n", mass_change);
  }


  if (TO_enabled){ //this section of the code deals with moving and then merging TO WFs

    if (verbosity.compare("high") == 0) {
      printf("states before TO WF movement: \n"); 
      listPrint(*head);
    }

    double mass_before_TO_move = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

    /*************************************************************************************/
    //first, have to check rare case where multiple surface WFs advance and both surpass a TO WF, such that a TO WF would be in between 2 surface WFs. In this case it's just easiest to merge the surface WFs.
    if (listLength_surface(*head)>0){
      lgarto_resolve_TO_WF_between_surf_WFs(soil_type, soil_properties, head);
    }
    /*************************************************************************************/


    /*************************************************************************************/
    //next, we calculate the contribution of the TO WFs to AET, and then we extract the AET from TO WFs, making them shorter. 
    double root_zone_depth = root_zone_depth_cm;

    if (verbosity.compare("high") == 0) {
      printf("root_zone_depth: %lf \n", root_zone_depth);
      printf("before calc AET from TO WFs: \n");
      listPrint(*head);
    }

    double cumulative_ET_from_TO_WFs_cm = lgarto_calc_aet_from_TO_WFs(num_layers, deepest_surf_depth_at_start, root_zone_depth, PET_timestep_cm, timestep_h, surf_frac_rz, 
                                                                      wilting_point_psi_cm, field_capacity_psi_cm, soil_type, cum_layer_thickness_cm, soil_properties, head);


    if (verbosity.compare("high") == 0) {
      printf("cumulative_ET_from_TO_WFs_cm: %10.16lf \n", cumulative_ET_from_TO_WFs_cm); 
      printf("after calc AET from TO WFs: \n");
      listPrint(*head);
    }
    
    double AET_in_TO_WFs_cm = 0.0;
    lgarto_ensure_rooting_zone_population(root_zone_depth, PET_timestep_cm, soil_type, soil_properties, head); 

    *AET_demand_cm += cumulative_ET_from_TO_WFs_cm;

    AET_in_TO_WFs_cm = cumulative_ET_from_TO_WFs_cm;
    /*************************************************************************************/


    /*************************************************************************************/
    //Here the dzdt values for GW wetting fronts are used to update TO WF depths. After this, the bottom_boundary_flux_cm is updated based on TO WF movement.
    if (verbosity.compare("high") == 0) {
      printf("before TO WF depth values are updated via dZ/dt: \n");
      listPrint(*head);
      printf("bottom_boundary_flux_cm: %lf \n", bottom_boundary_flux_cm);
    }
    bool TO_WFs_above_surface_WFs_flag = false;
    double bottom_boundary_flux_above_surface_WFs_cm = 0.0;
    int temp_count_surface_WF = 0;
    for (int wf = (listLength(*head)-1); wf != 1; wf--) { 

      if (wf<=0){
        break;
      }
      
      current = listFindFront(wf,*head, NULL);

      if (current==NULL){
        printf("current is null \n");
        printf("wf: %d \n", wf);
      }
      next = listFindFront(wf+1,*head, NULL);
      if (current->is_WF_GW==0){
        while (current->is_WF_GW==0){
          temp_count_surface_WF++;
          TO_WFs_above_surface_WFs_flag = true;

          if (wf-temp_count_surface_WF<=0){
            break;
          }
          current = listFindFront(wf-temp_count_surface_WF, *head, NULL);
        }
        wf = wf - temp_count_surface_WF;
      }

      if (wf-temp_count_surface_WF<=0){
        break;
      }

      double delta_depth;
      delta_depth = current->dzdt_cm_per_h * timestep_h;
      double delta_theta;
      if (current->layer_num == next->layer_num){
        delta_theta = next->theta - current->theta;
      }
      else{
        double equiv_next_theta;
        int soil_num = soil_type[current->layer_num];

        double theta_e   = soil_properties[soil_num].theta_e;
        double theta_r   = soil_properties[soil_num].theta_r;
        double vg_a      = soil_properties[soil_num].vg_alpha_per_cm;
        double vg_m      = soil_properties[soil_num].vg_m;
        double vg_n      = soil_properties[soil_num].vg_n;

        equiv_next_theta = calc_theta_from_h(next->psi_cm, vg_a, vg_m, vg_n, theta_e, theta_r);

        delta_theta = equiv_next_theta - current->theta;
        
      }

      if (TO_WFs_above_surface_WFs_flag==false){
        current->depth_cm += delta_depth;
      }

      if (TO_WFs_above_surface_WFs_flag){
        bottom_boundary_flux_above_surface_WFs_cm = bottom_boundary_flux_above_surface_WFs_cm + delta_depth * delta_theta; //30 July 2024 test takeout
      }
      bottom_boundary_flux_cm = bottom_boundary_flux_cm + delta_depth * delta_theta; //the + is due to the fact that this is a flux out of the model 

    }

    if (verbosity.compare("high") == 0) {
      printf("after TO WF depth values are updated via dZ/dt: \n");
      listPrint(*head);
      printf("associated mass: %lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
      printf("bottom_boundary_flux_cm: %lf \n", bottom_boundary_flux_cm);
    }

    if (isnan(bottom_boundary_flux_cm)){
      printf("bottom_boundary_flux_cm_temp is not a number \n");
      listPrint(*head);
      abort();
    }

    if (listLength_surface(*head)>0){
      bottom_boundary_flux_cm = lgarto_extract_TO_GW_flux_from_surface_WFs(&bottom_boundary_flux_above_surface_WFs_cm, bottom_boundary_flux_cm, AET_demand_cm, cum_layer_thickness_cm, soil_type, soil_properties, head);
      lgar_global_theta_update(bottom_boundary_flux_above_surface_WFs_cm, soil_type, soil_properties, head); 
      lgar_global_psi_update(soil_type, soil_properties, head);
    }
    /*************************************************************************************/


    /*************************************************************************************/
    //makes it so that the WF density is never too low; creates WFs near the bottom of the root zone. 
    //The idea behind this code is that TO WF merging will naturally lead to a large gap in psi between adjacent wetting fronts towards the bottom of the root zone when AET is significant.
    //Without this code, recharge calculations in dry soils are less accurate 
    current = *head;
    next = current->next;
    bool in_order = TRUE;
    for (int wf = 1; wf != (listLength(*head)); wf++){
      if (current->depth_cm>next->depth_cm){
        in_order = FALSE;
      }
      current = next;
      next = current->next;
      if (next==NULL){
        break;
      }
    }

    current = *head;
    next = current->next;
    double mass_added = 0.0;
    
    for (int wf = 1; wf != (listLength(*head)); wf++){
      struct wetting_front *new_corrective_WF;
      if (next==NULL){
        break;
      }
      if (next->depth_cm>cum_layer_thickness_cm[next->layer_num]){//there is a rare case where TO layer boundary crossing is necessary, but this code occurs before that. So, in the event that TO WF layer boundary crossing is necessary, then TO WF density ensuring at the root zone depth is deferred to another subtimestep.
        break;
      }
      if ( (current->is_WF_GW==1) && (current->depth_cm>0.0) && (current->psi_cm>(next->psi_cm + 0.99*cum_layer_thickness_cm[num_layers]) ) && (current->to_bottom==FALSE) && (in_order) && (current->psi_cm > cum_layer_thickness_cm[num_layers]) && (current->psi_cm<5*cum_layer_thickness_cm[num_layers]) ){
        //high accuracy version: change the factor 0.99 above to 0.1
        //the factor 0.99 in the line above has also been 0.5, it's not terribly important. Basically it will determine how frequently a new WF will be inserted, becasue a larger factor here will make it so that a larger gap is psi is required before a new WF is inserted.
        //reducing this factor (values as low as 0.05 have been explored) theoretically increases both computational expense and accuracy, but in practice for year long simulations the impact seems fairly small.
        double first_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

        int soil_num = soil_type[current->layer_num];

        double theta_e   = soil_properties[soil_num].theta_e;
        double theta_r   = soil_properties[soil_num].theta_r;
        double vg_a      = soil_properties[soil_num].vg_alpha_per_cm;
        double vg_m      = soil_properties[soil_num].vg_m;
        double vg_n      = soil_properties[soil_num].vg_n;

        double new_psi = (current->psi_cm + next->psi_cm)*0.5;
        if (new_psi>current->psi_cm){
          new_psi = (current->psi_cm + next->psi_cm)*0.5;
        }
        double new_theta = calc_theta_from_h(new_psi, vg_a, vg_m, vg_n, theta_e, theta_r);

        new_corrective_WF = listInsertFront(current->depth_cm,new_theta,current->front_num+1,current->layer_num,0,head,1);
  
        new_corrective_WF->to_bottom=FALSE;
        new_corrective_WF->psi_cm=new_psi;
        new_corrective_WF->is_WF_GW=TRUE;
        
        current->depth_cm = current->depth_cm - cum_layer_thickness_cm[num_layers]*4.E-4;
        if (current->front_num>1){
          if (current->depth_cm<listFindFront(current->front_num - 1, *head, NULL)->depth_cm){
            current->depth_cm = (current->next->depth_cm + listFindFront(current->front_num - 1, *head, NULL)->depth_cm)*0.5; 
          }
        }
        if (current->depth_cm<0.0){
          current->depth_cm = 0.0;
        }

        mass_added = lgar_calc_mass_bal(cum_layer_thickness_cm, *head) - first_mass;

        break;
      }

        current = next;
        next = current->next;

    }
    bottom_boundary_flux_cm -= mass_added;

    if (isnan(bottom_boundary_flux_cm)){
      printf("bottom_boundary_flux_cm_temp is not a number \n");
      listPrint(*head);
      abort();
    }
    /*************************************************************************************/


    /*************************************************************************************/
    //after TO WFs have moved due to dZ/dt and AET, and the correct component of TO AET has been extracted from surface WFs, then merging of TO and surface WFs can occur, as well as TO-TO merging, layer boundary crossing, etc.
    bool is_dry_over_wet_wf = lgar_check_dry_over_wet_wetting_fronts(*head);
    if (is_dry_over_wet_wf)
    lgar_fix_dry_over_wet_wetting_fronts(&mass_change, cum_layer_thickness_cm, soil_type, head, soil_properties);
    
    bool merged_in_non_top_layer = false; //the name is a bit too short, this is a boolean for surface - TO merging in particular
    int correction_type = lgarto_correction_type(num_layers, cum_layer_thickness_cm, head); //determines which type of correction will be necessary 

    if (verbosity.compare("high") == 0) {
      printf("correction_type before while loop: %d \n", correction_type);
      printf("mass before loop: %lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
    }

    while (correction_type!=0){
      
      if (verbosity.compare("high") == 0) {
        printf("correction_type at start of iteration within loop: %d \n", correction_type);
      }

      if (correction_type==1){
        double mass_before_surf_TO_merge = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        merged_in_non_top_layer = lgar_merge_surface_and_TO_wetting_fronts(merged_in_non_top_layer, num_layers, cum_layer_thickness_cm, head);
        bool did_a_WF_have_negative_depth = lgarto_correct_negative_depths(head);
        lgarto_cleanup_after_surface_TO_merging_in_layer_below_top(merged_in_non_top_layer, soil_type, soil_properties, head);
        bool close_psis = correct_close_psis(soil_type, soil_properties, head);

        if (did_a_WF_have_negative_depth || close_psis){
          bottom_boundary_flux_cm += (mass_before_surf_TO_merge - lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
          if (isnan(bottom_boundary_flux_cm)){
            printf("bottom_boundary_flux_cm_temp is not a number after surface-TO merging \n");
            listPrint(*head);
            abort();
          }
        }

        if (verbosity.compare("high") == 0) {
          printf("after lgar_merge_surface_and_TO_wetting_fronts and associated cleanup: \n");
          listPrint(*head);
          printf("mass after lgar_merge_surface_and_TO_wetting_fronts and associated cleanup: %lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
        }

        correction_type = lgarto_correction_type(num_layers, cum_layer_thickness_cm, head);

      }

      if (correction_type==2){
        double initial_mass = mass_before_TO_move - AET_in_TO_WFs_cm + mass_change - bottom_boundary_flux_cm + bottom_boundary_flux_cm_from_surf_WFs; //hm so maybe better to just look at mass from start of time step even before surfae stuff, and then account for all fluxes. because it's possible there's an error in the mass bal calc for example from having a surface WF that should merge with TO WFs but has not yet
        if (verbosity.compare("high") == 0) {
          printf("listPrint before iterative TO-TO depth merging: \n");
          listPrint(*head);
        }
        bottom_boundary_flux_cm += lgarto_TO_WFs_merge_via_depth(initial_mass, column_depth, cum_layer_thickness_cm, head, soil_type, soil_properties);
        lgar_global_psi_update(soil_type, soil_properties, head);
        correction_type = lgarto_correction_type(num_layers, cum_layer_thickness_cm, head);
        if (verbosity.compare("high") == 0) {
          printf("listPrint after iterative TO-TO depth merging: \n");
          listPrint(*head);
        }
      }

      if (correction_type==3){
        double initial_mass = mass_before_TO_move - AET_in_TO_WFs_cm + mass_change - bottom_boundary_flux_cm + bottom_boundary_flux_cm_from_surf_WFs; //possibly should just be bottom_boundary_flux_cm from TO WFs, technically some could come from surf WFs that got too deep
        bottom_boundary_flux_cm += lgarto_TO_WFs_merge_via_theta(initial_mass, column_depth, cum_layer_thickness_cm, head, soil_type, soil_properties);
            if (isnan(bottom_boundary_flux_cm)){
              printf("bottom_boundary_flux_cm_temp is not a number in correction type 3 (TO WFs merge via theta) \n");
              listPrint(*head);
              abort();
            }

        lgar_global_psi_update(soil_type, soil_properties, head);

        correction_type = lgarto_correction_type(num_layers, cum_layer_thickness_cm, head);

      }

      if (correction_type==4){
        int front_num_with_negative_depth = -1;
        lgar_TO_wetting_fronts_cross_layer_boundary(&front_num_with_negative_depth, num_layers, cum_layer_thickness_cm, soil_type, frozen_factor, soil_properties, head);
        lgar_global_psi_update(soil_type, soil_properties, head);

        if (front_num_with_negative_depth>-1){//in some uncommon cases where two adjacent TO WFs have very close psi values, the correction that sends a TO WF to a higher layer will give it a negative depth. This would lead to a small mass balance error which is corrected here. 
          double init_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
          listDeleteFront(front_num_with_negative_depth, head, soil_type, soil_properties);
          double after_del_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
          bottom_boundary_flux_cm += (init_mass - after_del_mass);
          if (verbosity.compare("high") == 0){
            printf("deleted negative WF as part of moving TO WF from a lower to higher layer ... \n");
            listPrint(*head);
          }
        }
        
        correction_type = lgarto_correction_type(num_layers, cum_layer_thickness_cm, head);
      }

      if (correction_type==5){
        lgar_merge_wetting_fronts(head, soil_type, frozen_factor, soil_properties);
        correction_type = lgarto_correction_type(num_layers, cum_layer_thickness_cm, head);
      }

      if (correction_type==6){
        double mass_corr_loop_start = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        bool close_psis = correct_close_psis(soil_type, soil_properties, head);
        if (close_psis){
          bottom_boundary_flux_cm += (mass_corr_loop_start - lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
        }
        lgar_wetting_fronts_cross_layer_boundary(num_layers, cum_layer_thickness_cm, soil_type, frozen_factor, head, soil_properties);
        correction_type = lgarto_correction_type(num_layers, cum_layer_thickness_cm, head);
      }

      if (correction_type==7){
        bottom_boundary_flux_cm += lgar_wetting_front_cross_domain_boundary(TO_enabled, cum_layer_thickness_cm, soil_type, frozen_factor, head, soil_properties);
          if (isnan(bottom_boundary_flux_cm)){
            printf("bottom_boundary_flux_cm_temp is not a number in correction type 7 (lgar_wetting_front_cross_domain_boundary) \n");
            listPrint(*head);
            abort();
            bottom_boundary_flux_cm = 0.0;
          }
        correction_type = lgarto_correction_type(num_layers, cum_layer_thickness_cm, head);
      }

      if (correction_type==8){//in some rare cases, a TO WF was deleted such that the deepest surface WF (or more) should become TO.
      //PTL: this is probably not necessary anymore, due to the fact that listDeleteFront now handles correcting to_bottom, is_WF_GW, and psi and theta corrections.
      //Still, I'm leaving it in for now. When I have time (currently towards the end of the OWP contract) I will test without this code. 
        struct wetting_front *temp_WF = *head;
        for (int wf = 1; wf != (listLength(*head)); wf++){
          if ((temp_WF==NULL)){
            break;
          }
          if (temp_WF->is_WF_GW==FALSE && temp_WF->next->is_WF_GW==TRUE && temp_WF->theta<temp_WF->next->theta && temp_WF->layer_num==temp_WF->next->layer_num){
            temp_WF->is_WF_GW = 1;
          }
          if (temp_WF->to_bottom==TRUE && temp_WF->is_WF_GW==FALSE && temp_WF->next->is_WF_GW==TRUE){
            temp_WF->is_WF_GW = 1;
          }
          temp_WF = temp_WF->next;
        }
        correction_type = lgarto_correction_type(num_layers, cum_layer_thickness_cm, head);
      }

      if (verbosity.compare("high") == 0) {
        printf("correction_type at end of iteration in while loop: %d \n", correction_type);
      }

    }

    if (verbosity.compare("high") == 0) {
      printf("correction_type at end of loop: %d \n", correction_type);
      printf("mass at end of loop: %lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
    }
    /*************************************************************************************/


    /*************************************************************************************/
    //LGAR and LGARTO create WFs with precipitation events and delete them with merging. When surface WFs get deep or dry enough, they become TO WFs.
    //Therefore, in some cases, there could more or less be a new TO wetting front per precipitation event. These end up being spaced very closely in theta, and this
    //resolution isn't necessary. This fxn keeps the total number of TO WFs relatively small.

    double mass_before_lgarto_consolidate_excessive_fronts = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
    lgarto_consolidate_excessive_fronts(cum_layer_thickness_cm, head, soil_type, soil_properties);
    double mass_after_lgarto_consolidate_excessive_fronts = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
    bottom_boundary_flux_cm += (mass_before_lgarto_consolidate_excessive_fronts - mass_after_lgarto_consolidate_excessive_fronts);

      if (isnan(bottom_boundary_flux_cm)){
        printf("bottom_boundary_flux_cm_temp is not a number after lgarto_consolidate_excessive_fronts \n");
        listPrint(*head);
        abort();
      }
    /*************************************************************************************/

    /*************************************************************************************/
    //this is a check to make sure that TO WFs above surface WFs have a depth of 0 cm.
    //TO WFs are allowed above the block of surface WFs, but all TO WFs above surface WFs must have a depth of 0, meaning these wetting fronts
    //extend from the gorundwater up to the soil surface.
    current = *head;
    int break_count = 0;
    for (int wf=1; wf != listLength(*head); wf++) {
      if ((current->depth_cm>0.0) && (current->is_WF_GW==TRUE) && (listLength_surface(*head)>0) ){
        listPrint(*head);
        printf("A GW WF above a surface WF has a nonzero depth \n");
        printf("current->depth_cm: %lf \n", current->depth_cm);
        printf("current->front_num: %d \n", current->front_num);
        abort();
      }
      if (current->is_WF_GW==FALSE){
        break_count = break_count + 1;
      }
      if (break_count==listLength_surface(*head)){
        break;
      }
      current = current->next;
    }
    /*************************************************************************************/

  } //this is the right bracket for: if TO mode is enabled.
  //After this bracket, some general checks are performed, for LGAR or LGARTO modes.


  /*************************************************************************************/
  //this is code that deletes WFs if they're drier than hydrostatic and at the surface. 
  //this code seems to improve some simulations but not others. Keeping in for now, and I get the idea behind it, but not sure it it's strictly necessary
  current = *head; 
  next = current->next;
  for (int wf=1; wf != listLength(*head); wf++) {
    if ( (current->depth_cm==0.0) && (listLength_surface(*head)>0) && (current->psi_cm>cum_layer_thickness_cm[num_layers]) ){
      current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
    }
    if (next==NULL){
      break;
    }
    current = next;
    next = current->next;
  } 
  /*************************************************************************************/


  *volin_cm = bottom_boundary_flux_cm;


  /*************************************************************************************/
  //fix any mass balance issues and dry-over-wet wetting fronts conditions. Although, is this necessary again?
  //it might not be, given that a correction_type loop now applies to surface WFs. And removing it seems to have no effect on USDA SCAN simulations.
  //still, leaving in for now, and will test removal more thoroughly later
  bool is_dry_over_wet_wf = lgar_check_dry_over_wet_wetting_fronts(*head);
  
  if (is_dry_over_wet_wf){
    lgar_fix_dry_over_wet_wetting_fronts(&mass_change, cum_layer_thickness_cm, soil_type, head, soil_properties);
    lgar_fix_dry_over_wet_wetting_fronts(&mass_change, cum_layer_thickness_cm, soil_type, head, soil_properties); // should really do a while loop that does it until not necessary 
  }

  if (verbosity.compare("high") == 0) {
    printf ("mass change/adjustment (dry_over_wet case) = %lf \n", mass_change);
  }

  /* sometimes merging can cause a slight mass balance error, to close the mass balance, the tiny error is
     swept into AET. This seems to be a reasonable solution rather than considering all possible cases which
     would probably be hard and time consuming. Some of these cases (e.g., rain or AET on totally saturated
     soil) are very rare and the error is small is magnitude. For LGARTO development, might have to revisit,
     but for LGAR this seems to work great.
  */

  if (fabs(mass_change) > 1.0E-7) {
    *AET_demand_cm = *AET_demand_cm - mass_change;
  }
  /*************************************************************************************/

  /***********************************************/
  // make sure all psi values are updated
  // there is a rare error where, after a wetting front crosses a layer boundary to a soil layer with a much more sensitive soil water retention curve, and the wetting front is near but not at saturation,
  // the barely unsatruated wetting front in the new layer will update its psi value as 0. This causes unequal psi values among adjacent wetting fronts in different layers and then a small mass balance error.
  // While that is ok for this soil layer in particular, adjacent wetting fronts above this one with a less sensitive soil water retention curve will yield a non-theta_e value for the psi value that is slightly above 0.
  // A solution is to either not run the following code, or to not run it when the wetting front is very close to saturation with a very sensitive soil water retention curve. Adding code that runs the following only for psi>1, because the problem occurs for very small nonzero psi (for example 0.001) 

  current = *head;

  for (int wf=1; wf != listLength(*head); wf++) {

    if (current->psi_cm>1.0){//with recent fix in soil fxns, value can probably be much closer to 0 but 1 works fine 
      int soil_num_k    = soil_type[current->layer_num];

      double theta_e_k   = soil_properties[soil_num_k].theta_e;
      double theta_r_k   = soil_properties[soil_num_k].theta_r;
      double vg_a_k      = soil_properties[soil_num_k].vg_alpha_per_cm;
      double vg_m_k      = soil_properties[soil_num_k].vg_m;
      double vg_n_k      = soil_properties[soil_num_k].vg_n;

      double Ksat_cm_per_h_k  = frozen_factor[current->layer_num] * soil_properties[soil_num_k].Ksat_cm_per_h;

      double Se = calc_Se_from_theta(current->theta,theta_e_k,theta_r_k);
      current->psi_cm = calc_h_from_Se(Se, vg_a_k, vg_m_k, vg_n_k); 
      current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h_k, vg_m_k);
    }
    current = current->next;
  }

  /***********************************************/


  /***********************************************/
  //Just a check to make sure that, when there is only 1 layer, than the existing wetting front is at the correct depth.
  //This might have been fixed with other debugging related to scenarios with just 1 layer where the wetting front is completely satruated. Not sure this is necessary.
  if (listLength(*head)==1) {
    if (current->depth_cm != cum_layer_thickness_cm[1]) {
      current->depth_cm = cum_layer_thickness_cm[1];
    }
  }
  /***********************************************/

  /***********************************************/
  // code that deletes a WF that is too high
  current = *head;
  next = current->next;
  int starting_list_length = listLength(*head);
  for (int wf=1; wf != starting_list_length; wf++) {
    if (next!=NULL){
      if (current->depth_cm < cum_layer_thickness_cm[current->layer_num - 1]){//There is a rare case when mass conservative wetting front merging will yield a wetting front with a depth above the layer it is in, which can occur when two wetting fronts are extremely close in moisture value. 
                                //In this case the small mass balance error is accounted for in recharge (was previously AET).
        double start_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        if ((next->depth_cm >= 0.0) && (next->to_bottom==FALSE)){
          current->depth_cm = next->depth_cm + 1.E-7;
        }
        else {
          current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
        }
        double end_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        bottom_boundary_flux_cm += (start_mass - end_mass);
      }
      current = next;
      next = current->next;
    }
    if (next==NULL){
      break;
    }
  }
  /***********************************************/

  /***********************************************/
  //apparently deletes GW fronts that go below bottom bdy (surface fronts that did this have already been handled). Was developed in this way because at some point TO WFs were allowed to oscillate around the model lower boundary, no longer allowed ... possibly this should moved into the correction_type code for LGARTO 
  current = listFindFront(listLength(*head) - 1, *head, NULL);
  next = NULL;
  if (listLength(*head)>1){
    next = current->next;
  }
  if (next!=NULL){
    for (int wf = listLength(*head) - 1; wf != -9999; wf--){
      if ( (current->depth_cm>next->depth_cm) && (current->layer_num==next->layer_num) && (current->layer_num==num_layers) ){
        bottom_boundary_flux_cm += (current->depth_cm - next->depth_cm)*(current->theta-next->theta); //wondering if this should be moved to runoff rather than bottom bdy flux //well, it's a lot rarer than you expected, so not too important, and moving to runoff only makes sense if the domain is completely saturated
        next->psi_cm = current->psi_cm;
        next->theta = current->theta;
        current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
        current = next;
        if (verbosity.compare("high") == 0) {
          printf("after deleted front too deep: \n");
          listPrint(*head);
        }
      }

      if (listLength(*head)==1 || (current->front_num == 1)){
        break;
      }

      current = listFindFront(current->front_num - 1, *head, NULL);
      next = current->next;
      if (current->to_bottom==TRUE){
        current->is_WF_GW = next->is_WF_GW;
      }
    }
  }
  *volin_cm = bottom_boundary_flux_cm; 
  /***********************************************/


  /***********************************************/
  if (isnan(bottom_boundary_flux_cm)){
    printf("bottom_boundary_flux_cm_temp is not a number \n");
    listPrint(*head);
    abort();
  }
  /***********************************************/

  /***********************************************/
  //misc code that prevents or causes a crash for various physically impossible scenarios 
  current = *head;
  next = current->next;
  for (int wf = 1; wf != (listLength(*head)); wf++){
    if (next!=NULL){
      if ( (current->to_bottom==1) && (next->layer_num==current->layer_num) ){
        printf("layer nums out of order \n");
        printf("look at front number %d \n",current->front_num);
        listPrint(*head);
        abort(); 
      }
    }
    if (next!=NULL){
      if((current->to_bottom==TRUE) && (next->is_WF_GW==TRUE) && (current->is_WF_GW==TRUE) && (fabs(current->psi_cm - next->psi_cm)>1e-5)){//so this technically works but it's probably not a great solution. Should instead fix unequal psis among layer interfaces when the problem occurs, which is related to when negative WF depths occur and then a WF gets sent to the surface that shouldn't.
        //oh also, this forces mass balance closure, at least before this point, because if there is a mass bal error, its put into *volin
        double current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        next->psi_cm = current->psi_cm;
        int soil_num_k  = soil_type[next->layer_num];
        
        double theta_e_k = soil_properties[soil_num_k].theta_e;
        double theta_r_k = soil_properties[soil_num_k].theta_r;
        double vg_a_k    = soil_properties[soil_num_k].vg_alpha_per_cm;
        double vg_m_k    = soil_properties[soil_num_k].vg_m;
        double vg_n_k    = soil_properties[soil_num_k].vg_n;
        next->theta = calc_theta_from_h(next->psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k); 
        double new_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        double mass_diff = new_mass - current_mass;
        bottom_boundary_flux_cm -= mass_diff;
        *volin_cm = bottom_boundary_flux_cm; 
      }

    }
    if (next!=NULL){
      if ( (current->layer_num!=next->layer_num) && (current->layer_num!=(next->layer_num - 1)) ){
        printf("layer nums out of order \n");
        printf("look at front number %d \n",current->front_num);
        listPrint(*head);
        abort();
      }
    }
    if (isnan(current->depth_cm)){
      printf("WF depth is not a number \n");
      printf("look at front number %d \n",current->front_num);
      listPrint(*head);
      abort();
    }
    if (current->theta<0.0){
      printf("WF moisture is less than 0 \n");
      printf("look at front number %d \n",current->front_num);
      listPrint(*head);
      abort();
    }
    if (isnan(current->psi_cm)){
      printf("WF capillary head is not a number \n");
      printf("look at front number %d \n",current->front_num);
      listPrint(*head);
      abort();
    }
    int soil_num_k  = soil_type[current->layer_num];
    double theta_e_k = soil_properties[soil_num_k].theta_e;
    double theta_r_k = soil_properties[soil_num_k].theta_r;
    if (current->theta<(theta_r_k - 1.E-3) || current->theta>(theta_e_k + 1.E-3)){//machine precision issues can make theta slightly larger than theta_e and that's ok
      printf("a wetting front is either greater than theta_e or less than theta_r for its layer \n");
      printf("look at front number: %d \n", current->front_num);
      listPrint(*head);
      abort();
    }
    current = next;
    next = current->next;
  }
  /***********************************************/

  /***********************************************/
  //very specific bug where the number of WFs is equal to the number of layers and the psi value is not equal among layers
  if (num_layers>1){
    current = *head;
    if ( (listLength(*head)==num_layers) && (current->psi_cm!=current->next->psi_cm) ){
      current->psi_cm = current->next->psi_cm;

      int soil_num_k  = soil_type[current->layer_num];
      
      double theta_e_k = soil_properties[soil_num_k].theta_e;
      double theta_r_k = soil_properties[soil_num_k].theta_r;
      double vg_a_k    = soil_properties[soil_num_k].vg_alpha_per_cm;
      double vg_m_k    = soil_properties[soil_num_k].vg_m;
      double vg_n_k    = soil_properties[soil_num_k].vg_n;
      current->theta = calc_theta_from_h(current->psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
    }
  }
  /***********************************************/


  /***********************************************/
  //code that just deletes a WF if it's too deep. We seem to have this a few times, but this instance is more general
  double mass_before_delete = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  bool delete_flag = false; 
  current = *head;
  next = current->next;
  for (int wf = 1; wf != (listLength(*head)); wf++){
    if (current->depth_cm>cum_layer_thickness_cm[num_layers]){
      current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
      delete_flag = true;
    }
    current = next;
    next = current->next;
    if (next==NULL){
      break;
    }
  }
  if (delete_flag){
    bottom_boundary_flux_cm += (mass_before_delete - lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
    *volin_cm = bottom_boundary_flux_cm; 
  }
  /***********************************************/
   

  if (isnan(bottom_boundary_flux_cm)){
    printf("bottom_boundary_flux_cm_temp is not a number \n");
    listPrint(*head);
    abort();
  }


  /***********************************************/
  //code that prevents depth=0.0 TO WFs from being drier than the shallowest depth>0.0 TO WF. Relevant when there are both surface WFs and depth=0.0 TO WFs.
  current = *head;
  next = current->next;
  struct wetting_front *last_depth_zero_TO_WF = listFindFront(listLength_TO_WFs_above_surface_WFs(*head), *head, NULL);
  struct wetting_front *shallowest_nonzero_depth_TO_WF = listFindFront(listLength_TO_WFs_above_surface_WFs(*head) + listLength_surface(*head) + 1, *head, NULL);
  if (last_depth_zero_TO_WF!=NULL && shallowest_nonzero_depth_TO_WF!=NULL && listLength_surface(*head)>0){
    while (last_depth_zero_TO_WF->psi_cm < shallowest_nonzero_depth_TO_WF->psi_cm){
      last_depth_zero_TO_WF = listDeleteFront(last_depth_zero_TO_WF->front_num, head, soil_type, soil_properties);
      last_depth_zero_TO_WF = listFindFront(listLength_TO_WFs_above_surface_WFs(*head), *head, NULL);
      if (last_depth_zero_TO_WF==NULL){
        break;
      }
    }
  }
  /***********************************************/

  return(bottom_boundary_flux_cm);

}


// ############################################################################################
/*
  the function merges wetting fronts; called from lgar_move_wetting_fronts.
*/
// ############################################################################################

extern void lgar_merge_wetting_fronts(struct wetting_front** head, int *soil_type, double *frozen_factor, 
				      struct soil_properties_ *soil_properties)
{
  struct wetting_front *current;
  struct wetting_front *next;
  struct wetting_front *next_to_next;
  current = *head; 
  next = current->next;

  for (int wf=1; wf != listLength(*head); wf++) {

    if ( (current->is_WF_GW==0) & (next->is_WF_GW==0) ){
    // local variables
    int layer_num, soil_num;

    layer_num   = current->layer_num;
    soil_num    = soil_type[layer_num];

    next = current->next; 
    next_to_next = current->next->next; 

    // case : wetting front passing another wetting front within a layer
    /**********************************************************/
    // 'current->depth_cm > next->depth_cm' ensures that merging is needed
    // 'current->layer_num == next->layer_num' ensures wetting fronts are in the same layer
    // '!next->to_bottom' ensures that the next wetting front is not the deepest wetting front in the layer
    
      if ( (current->depth_cm > next->depth_cm) && (current->theta>next->theta) && (current->layer_num == next->layer_num) && !next->to_bottom) {
        double current_mass_this_layer = current->depth_cm * (current->theta - next->theta) + next->depth_cm*(next->theta - next_to_next->theta);
        current->depth_cm = current_mass_this_layer / (current->theta - next_to_next->theta);

        double theta_e_k        = soil_properties[soil_num].theta_e;
        double theta_r_k        = soil_properties[soil_num].theta_r;
        double vg_a_k           = soil_properties[soil_num].vg_alpha_per_cm;
        double vg_m_k           = soil_properties[soil_num].vg_m;
        double vg_n_k           = soil_properties[soil_num].vg_n;
        double Ksat_cm_per_h_k  = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num];
        double Se_k             = calc_Se_from_theta(current->theta,theta_e_k,theta_r_k);
        current->psi_cm         = calc_h_from_Se(Se_k, vg_a_k, vg_m_k, vg_n_k);
        current->K_cm_per_h     = calc_K_from_Se(Se_k, Ksat_cm_per_h_k, vg_m_k); // AJ - K_temp in python version for 1st layer

        next = listDeleteFront(next->front_num, head, soil_type, soil_properties);

      }
    }
  current = current->next;
  next = current->next;
  if (next==NULL){
    break;
  }
  }
  if (verbosity.compare("high") == 0) {
    printf("after surface WFs merge: \n");
    listPrint(*head);
  }
}


// ############################################################################################
/*
  the function lets wetting fronts of a sufficient depth cross layer boundaries; called from lgar_move_wetting_fronts.
*/
// ############################################################################################

extern void lgar_wetting_fronts_cross_layer_boundary(int num_layers,
						     double* cum_layer_thickness_cm, int *soil_type,
						     double *frozen_factor, struct wetting_front** head,
						     struct soil_properties_ *soil_properties)
{
  struct wetting_front *current;
  struct wetting_front *next;
  struct wetting_front *next_to_next;
  bool cross_necessary = false;
  current = *head; 

  if (verbosity.compare("high") == 0) {
    printf("Inside layer boundary crossing... \n");
    printf("before surface WFs cross layer bdy: \n");
    listPrint(*head);
  }
  for (int wf=1; wf != listLength(*head); wf++) {
    if ( current->is_WF_GW==0 ){

    double theta_e,theta_r;
    double vg_a, vg_m, vg_n;
    int layer_num, soil_num;

    layer_num   = current->layer_num;
    soil_num    = soil_type[layer_num];
    theta_e     = soil_properties[soil_num].theta_e;
    theta_r     = soil_properties[soil_num].theta_r;
    vg_a        = soil_properties[soil_num].vg_alpha_per_cm;
    vg_m        = soil_properties[soil_num].vg_m;
    vg_n        = soil_properties[soil_num].vg_n;
    double Ksat_cm_per_h  = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num]; 

    next = current->next; 
    next_to_next = current->next->next; 

    if (current->depth_cm > cum_layer_thickness_cm[layer_num] && (next->depth_cm == cum_layer_thickness_cm[layer_num]) && (layer_num!=num_layers) && (current->theta>next->theta) ) { //PTL: took out && current->depth_cm <= column_depth, we actually don't want that. can be that a WF passes the one below, and has a huge dZ/dt value such that it ends up deeper than the model domain. replaced with (layer_num!=num_layers)
      cross_necessary = true;
      double current_theta = fmin(theta_e, current->theta);
      double overshot_depth = current->depth_cm - next->depth_cm;
      int soil_num_next = soil_type[layer_num+1];

      double next_theta_e   = soil_properties[soil_num_next].theta_e;
      double next_theta_r   = soil_properties[soil_num_next].theta_r;
      double next_vg_a      = soil_properties[soil_num_next].vg_alpha_per_cm;
      double next_vg_m      = soil_properties[soil_num_next].vg_m;
      double next_vg_n      = soil_properties[soil_num_next].vg_n;

      double Se = calc_Se_from_theta(current->theta,theta_e,theta_r);
      current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);

      current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m); 

      // current psi with van Gunechten properties of the next layer to get new theta
      double theta_new = calc_theta_from_h(current->psi_cm, next_vg_a, next_vg_m, next_vg_n, next_theta_e, next_theta_r);

      double mbal_correction = overshot_depth * (current_theta - next->theta);
      double mbal_Z_correction = mbal_correction / (theta_new - next_to_next->theta); // this is the new wetting front depth in the new layer
      if (next_to_next->theta==next_theta_e && next_to_next->psi_cm>0.0){//the idea here is that next_to_next will have a very small psi value and be in a layer where the psi-theta relationship loses its 1:1 nature, so this prevents an inf WF depth in some rare cases.
        next_to_next->theta -= 1.E-9;
        mbal_Z_correction = mbal_correction / (theta_new - (next_to_next->theta));
      }

      double depth_new = cum_layer_thickness_cm[layer_num] + mbal_Z_correction; // this is the new wetting front absolute depth

        if (isinf(depth_new)){
          theta_new = theta_new - 1.E-12;
          mbal_correction = overshot_depth * (current_theta - next->theta);
          mbal_Z_correction = mbal_correction / (theta_new - next_to_next->theta); // this is the new wetting front depth

          depth_new = cum_layer_thickness_cm[layer_num] + mbal_Z_correction; // this is the new wetting front absolute depth
        }

        if (isnan(depth_new) || isinf(depth_new)){

          depth_new = cum_layer_thickness_cm[num_layers]+1.0;

        }
        current->depth_cm = cum_layer_thickness_cm[layer_num];

        next->theta = theta_new;
        next->psi_cm = current->psi_cm;
        next->depth_cm = depth_new;
        next->layer_num = layer_num + 1;
        next->dzdt_cm_per_h = current->dzdt_cm_per_h;
        current->dzdt_cm_per_h = 0;
        current->to_bottom = TRUE;
        current->is_WF_GW = 0;
        next->to_bottom = FALSE;
        next->is_WF_GW = 0;

        if ( (next->depth_cm < cum_layer_thickness_cm[layer_num]) ){// this is the event where, in a single time step, for example due to AET, a wetting front got deep enough to go to the next layer but also dry enough to have to merge with the next WF.
          current->is_WF_GW = next->next->is_WF_GW;
          next->is_WF_GW = next->next->is_WF_GW;
        }
      }
    }
  current = current->next;
  }

  if (cross_necessary){
    current = *head;
    for (int wf=1; wf != listLength(*head); wf++) {
      if ( (current->to_bottom==TRUE) && (current->next!=NULL) ){
        current->is_WF_GW = current->next->is_WF_GW;
      }
      current = current->next;
    }

    for (int wf = listLength(*head)-1; wf != 0; wf--) {
      struct wetting_front *current_temp = listFindFront(wf, *head, NULL);
      struct wetting_front *next_temp = current_temp->next;
      if ( (current_temp->to_bottom==TRUE) ){
        current_temp->is_WF_GW = next_temp->is_WF_GW;
        current_temp->psi_cm = next_temp->psi_cm;

        int soil_num_k1 = soil_type[current_temp->layer_num]; 
        double theta_e_k   = soil_properties[soil_num_k1].theta_e;
        double theta_r_k   = soil_properties[soil_num_k1].theta_r;
        double vg_a_k      = soil_properties[soil_num_k1].vg_alpha_per_cm;
        double vg_m_k      = soil_properties[soil_num_k1].vg_m;
        double vg_n_k      = soil_properties[soil_num_k1].vg_n;
        current_temp->theta = calc_theta_from_h(current_temp->psi_cm, vg_a_k, vg_m_k, vg_n_k,theta_e_k,theta_r_k);
      }
    }

  }

  if (verbosity.compare("high") == 0) {
    printf("after surface WFs cross layer bdy: \n");
    listPrint(*head);
  }  
}


// ############################################################################################
/*
  the function lets wetting fronts of a sufficient depth interact with the lower boundary; called from lgar_move_wetting_fronts.
*/
// ############################################################################################

extern double lgar_wetting_front_cross_domain_boundary(bool TO_enabled, double* cum_layer_thickness_cm, int *soil_type, double *frozen_factor,
						       struct wetting_front** head, struct soil_properties_ *soil_properties)
{
  struct wetting_front *current;
  struct wetting_front *next;
  struct wetting_front *next_to_next;
  current = *head; 
  double bottom_flux_cm = 0.0;

  for (int wf=1; wf != listLength(*head); wf++) {

    int layer_num, soil_num;
    double bottom_flux_cm_temp=0.0;

    layer_num   = current->layer_num;
    soil_num    = soil_type[layer_num];

    next = current->next; 
    next_to_next = current->next->next; 

    // case : wetting front is the deepest one in the last layer (most deepested wetting front in the domain)
    /**********************************************************/
    bool break_flag = FALSE;
    if (next_to_next == NULL && current->depth_cm > cum_layer_thickness_cm[layer_num] && current->is_WF_GW==0) {
      break_flag = TRUE;
      //  this is the water leaving the system through the bottom of the soil
      bottom_flux_cm_temp = (current->theta - next->theta) *  (current->depth_cm - next->depth_cm);
      double theta_e_k   = soil_properties[soil_num].theta_e;
      double theta_r_k   = soil_properties[soil_num].theta_r;
      double vg_a_k      = soil_properties[soil_num].vg_alpha_per_cm;
      double vg_m_k      = soil_properties[soil_num].vg_m;
      double vg_n_k      = soil_properties[soil_num].vg_n;
      double Ksat_cm_per_h_k  = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num];

      next->theta = current->theta;
      double Se_k = calc_Se_from_theta(current->theta,theta_e_k,theta_r_k);
      next->psi_cm = calc_h_from_Se(Se_k, vg_a_k, vg_m_k, vg_n_k);
      next->K_cm_per_h = calc_K_from_Se(Se_k, Ksat_cm_per_h_k, vg_m_k);
      if (TO_enabled){
        current->next->is_WF_GW = 1;
      }
      current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
    }

    bottom_flux_cm += bottom_flux_cm_temp;
    if (break_flag){
     break;
    }
    current = current->next;

    }
  if (verbosity.compare("high") == 0) {
    printf("after surface WFs cross model lower boundary: \n");
    listPrint(*head);
  }
  return bottom_flux_cm;

}


// ############################################################################################
/* The function allows surface wetting fronts to merge with groundwater WFs when the surface 
   wetting fronts get deep enough. */
// ############################################################################################

extern bool lgar_merge_surface_and_TO_wetting_fronts(bool merged_in_non_top_layer, int num_layers, double* cum_layer_thickness_cm, struct wetting_front** head)
{
  if (verbosity.compare("high") == 0) {
    printf("surface-TO merging start\n");
    printf("current mass: %.17lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
  }
  for (int wf = 1; wf != (listLength(*head)); wf++) { //technically this loop can stop once we're at the first TO WF below the surface WFs
    struct wetting_front *current = listFindFront(wf, *head, NULL);
    struct wetting_front *next = current->next;

    if ( ((current->depth_cm>next->depth_cm) && (current->is_WF_GW==0) && (next->is_WF_GW==1) && (next->to_bottom==0)) ){
      if (current->layer_num>1)
        {merged_in_non_top_layer = true;}
      if (num_layers==1){
        merged_in_non_top_layer = false;
      }

      if (current->theta>next->next->theta){
        if (verbosity.compare("high") == 0) {
          printf("currently surface-TO merging, case where surface WF is wetter. current->front_num and then listprint: \n");
          printf("current->front_num: %d \n",current->front_num);
          printf("lgar_calc_mass_bal: %lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
          listPrint(*head);
        }

        double theta_next_temp = next->theta;
        int to_bottom_next_temp = next->to_bottom;
        double dzdt_cm_per_h_next_temp = next->dzdt_cm_per_h;
        double K_cm_per_h_next_temp = next->K_cm_per_h;
        double psi_cm_next_temp = next->psi_cm;

        next->depth_cm = current->depth_cm + (current->depth_cm - next->depth_cm)*(next->next->theta-next->theta)/fabs(current->theta-next->next->theta);

        next->theta = current->theta;
        next->to_bottom = current->to_bottom;
        next->dzdt_cm_per_h = current->dzdt_cm_per_h;
        next->K_cm_per_h = current->K_cm_per_h;
        next->is_WF_GW = 0;
        next->psi_cm = current->psi_cm;

        current->depth_cm = 0;
        current->theta = theta_next_temp;
        current->to_bottom = to_bottom_next_temp ;
        current->dzdt_cm_per_h = dzdt_cm_per_h_next_temp;
        current->K_cm_per_h = K_cm_per_h_next_temp;
        current->is_WF_GW = 1; 
        current->psi_cm = psi_cm_next_temp;

      }
      else{
        if (verbosity.compare("high") == 0) {
          printf("currently surface-TO merging, case where surface WF is drier. current->front_num then next and then listprint: \n");
          printf("current->front_num: %d \n",current->front_num);
          printf("next->front_num: %d \n",next->front_num);
          listPrint(*head);
        }
        
        current->is_WF_GW = 1;
        double overshot_depth = current->depth_cm - next->depth_cm;
        current->depth_cm = next->depth_cm;
        next->depth_cm = current->depth_cm - (overshot_depth)*(current->theta-next->theta)/fabs(next->next->theta-current->theta);

        double theta_next_temp = next->theta;
        double dzdt_cm_per_h_next_temp = next->dzdt_cm_per_h;
        double K_cm_per_h_next_temp = next->K_cm_per_h;
        double psi_cm_next_temp = next->psi_cm;

        next->theta = current->theta;
        next->dzdt_cm_per_h = current->dzdt_cm_per_h;
        next->K_cm_per_h = current->K_cm_per_h;
        next->psi_cm = current->psi_cm;

        // double depth_temp = current->depth_cm;
        current->depth_cm = 0;
        current->theta = theta_next_temp;
        current->dzdt_cm_per_h = dzdt_cm_per_h_next_temp;
        current->K_cm_per_h = K_cm_per_h_next_temp;
        current->psi_cm = psi_cm_next_temp;

      }
    }
      next = current->next;
  }


  for (int wf = listLength(*head)-1; wf != 0; wf--) {
    struct wetting_front *current = listFindFront(wf, *head, NULL);
    struct wetting_front *next = current->next;
    if ( (current->to_bottom==TRUE) && (next->is_WF_GW!=current->is_WF_GW) ){
      current->is_WF_GW = next->is_WF_GW;
    }
  }

  listSendToTop(*head);

  if (verbosity.compare("high") == 0) {
    printf("surface-TO merging end\n");
    printf("current mass: %.17lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
  }

  return merged_in_non_top_layer;

}


// ############################################################################################
/* The function handles situation of dry over wet wetting fronts
  mainly happen when AET extracts more water from the upper wetting front
  and the front gets drier than the lower wetting front */
// ############################################################################################
extern void lgar_fix_dry_over_wet_wetting_fronts(double *mass_change, double* cum_layer_thickness_cm, int *soil_type,
					 struct wetting_front** head, struct soil_properties_ *soil_properties)
{
  if (verbosity.compare("high") == 0) {
    printf("start of Fix Dry over Wet Wetting Front... \n");
    listPrint(*head);
    printf("mass: %lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
  }

  struct wetting_front *current;
  struct wetting_front *next;
  current = *head;
  next = current->next;
  bool deleted_front = false;

  for (int l=1; l <= listLength(*head); l++) {
    if ((next != NULL) && (current->is_WF_GW==0)) {
      if (next->is_WF_GW==0){
        // this part fixes case of upper theta less than lower theta due to AET extraction
        // also handles the case when the current and next wetting fronts have the same theta
        // and are within the same layer
        /***************************************************/

        if ( (current->theta <= next->theta) && (current->layer_num == next->layer_num) ) {
    double mass_before = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);

    if (current->layer_num == 1 && next->to_bottom==FALSE){
      double depth_above  = 0.0;
      if (current->front_num>1){
        depth_above = listFindFront(current->front_num - 1, *head, NULL)->depth_cm;
      }
      double thickness_current_front = current->depth_cm - depth_above;
      double thickness_next_front = next->depth_cm - current->depth_cm;

      double target_parital_mass = thickness_current_front*current->theta + thickness_next_front*next->theta;

      next->theta = target_parital_mass / (thickness_next_front + thickness_current_front);

      int soil_num_k  = soil_type[current->layer_num];
      double theta_e_k = soil_properties[soil_num_k].theta_e;
      double theta_r_k = soil_properties[soil_num_k].theta_r;
      double vg_a_k    = soil_properties[soil_num_k].vg_alpha_per_cm;
      double vg_m_k    = soil_properties[soil_num_k].vg_m;
      double vg_n_k    = soil_properties[soil_num_k].vg_n;
      double Se = calc_Se_from_theta(current->theta,theta_e_k,theta_r_k);

      next->psi_cm = calc_h_from_Se(Se, vg_a_k, vg_m_k, vg_n_k); 

      if (next->next!=NULL){
        if (listLength_surface(*head)==2 && next->is_WF_GW==FALSE && next->layer_num==next->next->layer_num && next->next->theta>next->theta){
          next->is_WF_GW = TRUE;
        }
      }
    }

    current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
    deleted_front = true;
    double mass_after = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
    *mass_change += fabs(mass_after - mass_before);

    /* note: mass_before is less when we have wetter front over drier front condition,
      however, lgar_calc_mass_bal returns mass_before > mass_after due to fabs(theta_current - theta_next);
      for mass_before the functions compuates more than the actual mass; removing fabs in that function
      might be one option, but for now we are adding fabs to mass_change to make sure we added extra water
      back to AET after deleting the drier front */

        }

      }
    }
    current = current->next;

        if (current == NULL)
    next = NULL;
        else
    next = current->next;
  }

  if (deleted_front){
    current = *head;
    next = current->next;
    for (int l=1; l <= listLength(*head); l++) {
      if (isnan(current->psi_cm)){
        current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
        current = next;
        next = current->next;
      }
      if (next==NULL){
        break;
      }
    }
  }

  if (verbosity.compare("high") == 0) {
    printf("end of Fix Dry over Wet Wetting Front... \n");
    listPrint(*head);
    printf("mass: %lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
  }

}

// ############################################################################################
/* The function handles situation of dry over wet wetting fronts
  mainly happen when AET extracts more water from the upper wetting front
  and the front gets drier than the lower wetting front */
// ############################################################################################
extern bool lgar_check_dry_over_wet_wetting_fronts(struct wetting_front* head)
{
  struct wetting_front *current = head;
  struct wetting_front *next    = current->next;
  int length = listLength(head);
  
  for (int l=1; l <= length; l++) {
    if (next != NULL) {
      
      if ( (current->theta <= next->theta) && (current->layer_num == next->layer_num) && (current->is_WF_GW==FALSE) && (next->is_WF_GW==FALSE) )
	return true;
      
      current = current->next;
      
      if (current == NULL)
	next = NULL;
      else
	next = current->next;
      
    }
    
  }
  
  return false;
}
      
// ############################################################################################
/* The module computes the potential infiltration capacity, fp (in the lgar manuscript),
   potential infiltration capacity = the maximum amount of water that can be inserted into
   the soil depending on the availability of water.
   this module is called when a new superficial wetting front is not created
   in the current timestep, that is precipitation in the current and previous
   timesteps was greater than zero */
// ############################################################################################
extern double lgar_insert_water(bool use_closed_form_G, int nint, double timestep_h, double *free_drainage_subtimestep_cm, double *AET_demand_cm, double *ponded_depth_cm,
				double *volin_this_timestep, double precip_timestep_cm, int wf_free_drainage_demand,
			        int num_layers, double ponded_depth_max_cm, int *soil_type,
				double *cum_layer_thickness_cm, double *frozen_factor,
				struct wetting_front* head, struct soil_properties_ *soil_properties)
{
  // note ponded_depth_cm is a pointer.   Access its value as (*ponded_depth_cm).

  int wf_that_supplies_free_drainage_demand = wf_free_drainage_demand;

  // local vars
  double theta_e, theta_r;
  double vg_a, vg_m, vg_n,Ksat_cm_per_h;
  double h_min_cm;
  struct wetting_front *current;
  struct wetting_front *current_free_drainage;
  struct wetting_front *current_free_drainage_next;
  int soil_num;
  double f_p = 0.0;
  double runoff = 0.0;

  double h_p = fmax(*ponded_depth_cm - precip_timestep_cm * timestep_h, 0.0); // water ponded on the surface

  current = head;
  current_free_drainage      = listFindFront(wf_that_supplies_free_drainage_demand, head, NULL);
  current_free_drainage_next = listFindFront(wf_that_supplies_free_drainage_demand+1, head, NULL);

  int number_of_wetting_fronts = listLength(head);

  //int last_wetting_front_index = number_of_wetting_fronts;
  int layer_num_fp = current_free_drainage->layer_num;


  double Geff;

  if (number_of_wetting_fronts == num_layers) {
    Geff = 0.0; // i.e., case of no capillary suction, dz/dt is also zero for all wetting fronts
    soil_num = soil_type[layer_num_fp];
    Ksat_cm_per_h = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num]; //23 feb 2024
  }
  else {

    double theta_below = current_free_drainage_next->theta;

    soil_num = soil_type[layer_num_fp];

    theta_e = soil_properties[soil_num].theta_e;  // rhs of the new front, assumes theta_e as per Peter
    theta_r = soil_properties[soil_num].theta_r;
    h_min_cm = soil_properties[soil_num].h_min_cm;
    vg_a     = soil_properties[soil_num].vg_alpha_per_cm;
    vg_m     = soil_properties[soil_num].vg_m;
    vg_n     = soil_properties[soil_num].vg_n;
    double lambda = soil_properties[soil_num].bc_lambda;
    double bc_psib_cm = soil_properties[soil_num].bc_psib_cm;
    Ksat_cm_per_h = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num];

    Geff = calc_Geff(use_closed_form_G, theta_below, theta_e, theta_e, theta_r, vg_a, vg_n, vg_m, h_min_cm, Ksat_cm_per_h, nint, lambda, bc_psib_cm); 

  }

  // if the free_drainage wetting front is the top most, then the potential infiltration capacity has the following simple form
  if (layer_num_fp == 1) {
      f_p = Ksat_cm_per_h * (1 + (Geff + h_p)/current_free_drainage->depth_cm);
  }
  else {
    // see the paper "Layered Green and Ampt Infiltration With Redistribution" by La Follette et al. (https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1029/2022WR033742), equations 16 or 19
    double bottom_sum = (current_free_drainage->depth_cm - cum_layer_thickness_cm[layer_num_fp-1])/Ksat_cm_per_h;

    for (int k = 1; k < layer_num_fp; k++) {
      int soil_num_k = soil_type[layer_num_fp-k];
      double Ksat_cm_per_h_k = soil_properties[soil_num_k].Ksat_cm_per_h * frozen_factor[layer_num_fp - k];

      bottom_sum += (cum_layer_thickness_cm[layer_num_fp - k] - cum_layer_thickness_cm[layer_num_fp - (k+1)])/ Ksat_cm_per_h_k;
    }

    f_p = (current_free_drainage->depth_cm / bottom_sum) + ((Geff + h_p)*Ksat_cm_per_h/(current_free_drainage->depth_cm)); //Geff + h_p

  }

  // checkpoint # AJ

  // if free drainge has to be included, which currently we don't, then the following will be set to hydraulic conductivity
  // of the deeepest layer
  if ((layer_num_fp == num_layers) && (current_free_drainage->psi_cm < 1.E-1) && (num_layers == number_of_wetting_fronts))
    f_p = fmin(f_p, *AET_demand_cm/timestep_h); //the idea here is that, if the soil is nearly completely saturated, a little bit of water can still enter if AET is sufficiently large 

  //this code checks if there is enough storage available for infiltrating water. That is, f_p can only be as big as there is room for water. 
  double max_storage = 0.0;
  for (int k = 1; k < num_layers+1; k++) {
    int layer_num = k;
    soil_num = soil_type[layer_num];
    max_storage += soil_properties[soil_num].theta_e * (cum_layer_thickness_cm[k]-cum_layer_thickness_cm[k-1]);
  }// would be more efficient to make max_storage a vairable that is not computed every time lgar_insert_water is called

  // if (f_p*timestep_h + current_mass > max_storage){
  //   f_p = f_p*timestep_h + current_mass - max_storage;
  // }//19 August 2024 test
  /////the interesting thing here is that there are arguments both for and against using the above 3 lines.
  /////for: the vadose zone can't hold more water than it currently has available for storage, but
  /////against: the recharge fluxes have not yet been extracted, so in reality the above code underestimates f_p.
  /////ultimately either choice is defensible and the impact becomes smaller as the time step gets smaller. 

  // if ( (current_mass)/max_storage > 0.99 and also if TO mode is off ){
  //   printf("warning: vadose zone is 99 percent full or greater. If you are using the model in an environment with more precipitation than PET, LGAR is not an appropriate model because its lower boundary condition is no flow. \n ");
  // } //turning off for now; should really print like once per model run or so, not every time step 

  double ponded_depth_temp = *ponded_depth_cm;

  double free_drainage_demand = *free_drainage_subtimestep_cm;
  

  // 'if' condition is not needed ... AJ
  if ((layer_num_fp==num_layers) && (num_layers == number_of_wetting_fronts))
    ponded_depth_temp = *ponded_depth_cm - f_p * timestep_h - 0*free_drainage_demand;
  else
    ponded_depth_temp = *ponded_depth_cm - f_p * timestep_h - 0*free_drainage_demand;

  ponded_depth_temp   = fmax(ponded_depth_temp, 0.0);

  double fp_cm = f_p * timestep_h + 0*free_drainage_demand; // infiltration in cm

  if (num_layers == listLength(head) && (head->psi_cm<1.E-2)){
    ponded_depth_temp = ponded_depth_temp - free_drainage_demand;
    fp_cm = fp_cm + free_drainage_demand;
  }

  if (ponded_depth_max_cm > 0.0 ) {

    if (ponded_depth_temp < ponded_depth_max_cm) {
      runoff = 0.0;
      *volin_this_timestep = fmin(*ponded_depth_cm, fp_cm); 
      *ponded_depth_cm     = *ponded_depth_cm - *volin_this_timestep;

        if (verbosity.compare("high") == 0){
        printf("fp_cm: %lf \n", fp_cm);
        printf("*volin_this_timestep: %lf \n", *volin_this_timestep);
        printf("runoff: %lf \n", runoff);
        printf("*ponded_depth_cm: %lf \n", *ponded_depth_cm);

      }
      return runoff;
    }
    else if (ponded_depth_temp > ponded_depth_max_cm ) {
      runoff = ponded_depth_temp - ponded_depth_max_cm;
      *ponded_depth_cm     = ponded_depth_max_cm;
      *volin_this_timestep = fp_cm;

        if (verbosity.compare("high") == 0){
          printf("fp_cm: %lf \n", fp_cm);
          printf("*volin_this_timestep: %lf \n", *volin_this_timestep);
          printf("runoff: %lf \n", runoff);
          printf("*ponded_depth_cm: %lf \n", *ponded_depth_cm);

        }

      return runoff;
    }

  }

  else {
    // if it got to this point, no ponding is allowed, either infiltrate or runoff
    // order is important here; assign zero to ponded depth once we compute volume in and runoff
    *volin_this_timestep = fmin(*ponded_depth_cm, fp_cm); //
    runoff = *ponded_depth_cm < fp_cm ? 0.0 : (*ponded_depth_cm - *volin_this_timestep);
    *ponded_depth_cm = 0.0;

  }

  if (verbosity.compare("high") == 0){
    printf("fp_cm: %lf \n", fp_cm);
    printf("*volin_this_timestep: %lf \n", *volin_this_timestep);
    printf("runoff: %lf \n", runoff);
    printf("*ponded_depth_cm: %lf \n", *ponded_depth_cm);

  }

  return runoff;
}

// ######################################################################################
/* This subroutine is called iff there is no surfacial front, it creates a new front and
   inserts ponded depth, and will return some amount if can't fit all water
   into the soil.  Note ponded_depth_cm is a pointer.   Access its value as (*ponded_depth_cm). */
// ######################################################################################
extern double lgar_create_surficial_front(bool TO_enabled, int num_layers, double *ponded_depth_cm, double *volin, double dry_depth,
					double theta1, int *soil_type, double *cum_layer_thickness_cm,
					double *frozen_factor, struct wetting_front** head, struct soil_properties_ *soil_properties)
{

  double bottom_boundary_flux_cm_temp = 0.0;
  if (verbosity.compare("high") == 0) {
    printf("before WF creation: \n");
    listPrint(*head);
  }

  // local vars
  double theta_e,Se,theta_r;
  double delta_theta;
  double vg_alpha_per_cm, vg_m, vg_n, Ksat_cm_per_h;
  double prior_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  bool new_wf_deleted_flag = false;

  bool to_bottom = FALSE;
  struct wetting_front *current;
  struct wetting_front *next;
  int layer_num,soil_num,front_num;

  current = *head;

  if (listLength_surface(*head)>0){//when TO mode is on, there are TO WFs that can have a depth of 0 and are not relevant for new wetting front depth calc
    while (current->is_WF_GW==1){
      current = current->next;
    }
  }

  if (current->depth_cm==0){//same idea as above block -- wetting fronts with 0 depth are possible but not relevant for new WF
    while (current->depth_cm==0){
      current = current->next;
    }
  }

  next = current->next;

  layer_num = 1;   // we only create new surfacial fronts in the first layer
  soil_num = soil_type[layer_num];
  front_num = 1;   // we are creating a new surfacial front, which by definition must be front #1

  theta_e = soil_properties[soil_num].theta_e;
  theta_r = soil_properties[soil_num].theta_r;
  delta_theta =  theta_e - theta1;

  double theta_new = 0.0;

  if(dry_depth * delta_theta > (*ponded_depth_cm))  // all the ponded depth enters the soil
    {
      *volin = *ponded_depth_cm;
      theta_new = fmin(theta1 + (*ponded_depth_cm) /dry_depth, theta_e);
      // theta_new = theta_e; 
      listInsertFirst(dry_depth, theta_new, front_num, layer_num, to_bottom, head, false);
      *ponded_depth_cm = 0.0;
      //hp_cm =0.0;
    }
  else  // not all ponded depth fits in
    {
      *volin = dry_depth * delta_theta;
      *ponded_depth_cm -= dry_depth * delta_theta;
      theta_new = theta_e; //fmin(theta1 + (*ponded_depth_cm) /dry_depth, theta_e);
      if (dry_depth < cum_layer_thickness_cm[1])
	listInsertFirst(dry_depth, theta_e, front_num, layer_num, to_bottom, head, false);
      else
  listInsertFirst(dry_depth, theta_e, front_num, layer_num, to_bottom, head, false); //the idea here is that a new WF should never have to_bottom as 1 -- if it needs to merge with the one below it, it will
      //hp_cm = *ponded_depth_cm;
    }

  current = *head;
  while (current->theta!=theta_new){
    current = current->next;
  }
  int new_wf_num = current->front_num;

  current = *head;  // must do this again because listInsertFirst() created a new *head
  vg_alpha_per_cm    = soil_properties[soil_num].vg_alpha_per_cm;
  vg_m               = soil_properties[soil_num].vg_m;
  vg_n               = soil_properties[soil_num].vg_n;
  Ksat_cm_per_h      = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[layer_num];

  Se = calc_Se_from_theta(theta_new,theta_e,theta_r);
  current->psi_cm = calc_h_from_Se(Se, vg_alpha_per_cm , vg_m, vg_n);

  current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m) * frozen_factor[layer_num]; // AJ - K_temp in python version for 1st layer

  current->dzdt_cm_per_h = 0.0; //for now assign 0 to dzdt as it will be computed/updated in lgar_dzdt_calc function

  int current_front_num_temp = current->front_num;
  
  if (current->next!=NULL){// sometimes a new WF immediately has to merge with another WF or cross a layer bdy
    if (current->depth_cm>current->next->depth_cm){
      //do a merge cross merge
      // I'm thinking it might not be necessary -- this should happen later. Leaving in for now however
      lgar_merge_wetting_fronts(head, soil_type, frozen_factor, soil_properties);
      double mass_corr_loop_start = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
      bool close_psis = correct_close_psis(soil_type, soil_properties, head);
      if (close_psis){
        bottom_boundary_flux_cm_temp += (mass_corr_loop_start - lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
      }
      lgar_wetting_fronts_cross_layer_boundary(num_layers, cum_layer_thickness_cm, soil_type, frozen_factor,
                head, soil_properties);
    }
  }


  /*************************************************************************************/
  //the rest of the fxn lgar_create_surficial_front is code dealing with when TO mode is on. 
  //we check if the new WF is deeper than any TO WFs.
  //Then the mass is updated accordingly. 
  /*************************************************************************************/


  /*************************************************************************************/
  //first, the new WF is moved into its correct position based on depth if necessary. It often will be in TO mode.
  current = listFindFront(current_front_num_temp,*head,NULL);
  next = current->next;
  bool had_to_sort = false;
  if (TO_enabled==1 && next!=NULL){
    if (current->depth_cm>next->depth_cm){
      had_to_sort = true;
      listSortFrontsByDepth(*head); 

      lgarto_resolve_TO_WF_between_surf_WFs(soil_type, soil_properties, head); //do the thing here  

    }
    if (verbosity.compare("high") == 0) {
      printf("list after listSortFrontsByDepth: \n");
      listPrint(*head);
    }
  }

  double mass_change = 0.0;
  bool is_dry_over_wet_wf = lgar_check_dry_over_wet_wetting_fronts(*head);
  
  if (is_dry_over_wet_wf)
  lgar_fix_dry_over_wet_wetting_fronts(&mass_change, cum_layer_thickness_cm, soil_type, head, soil_properties);

  /*************************************************************************************/


  /*************************************************************************************/
  //if sorting was necessary, then new WF is allowed to pass all WFs that it should based on its depth, and then its moisture is updated 
  //such that the psis are in the correct order with respect to depth, and the mass balance is conserved. 
  if (had_to_sort==true){
    double target_mass = prior_mass + *volin;
    double max_storage_possible = 0.0;

    for (int k = 1; k < num_layers+1; k++) {
      layer_num = k;
      soil_num = soil_type[layer_num];
      max_storage_possible += soil_properties[soil_num].theta_e * (cum_layer_thickness_cm[k]-cum_layer_thickness_cm[k-1]);
    }

    if (target_mass > max_storage_possible){
      *volin -= target_mass - max_storage_possible;
      *ponded_depth_cm += target_mass - max_storage_possible;
      target_mass = max_storage_possible;
    }

    //if a new WF was deep enough, higher TO WFs should be sent to surface 
    if (verbosity.compare("high") == 0) {
      printf("before sending WFs up, during new WF creation process: \n");
      listPrint(*head);
    }
    current = *head;
    next = current->next;
    for (int wf = 1; wf != (listLength(*head)); wf++) {
      if ( (listLength_surface(*head)>0) && (current->is_WF_GW==1) ){
        current->depth_cm = 0.0; //note that this should never set a to_bottom == TRUE WF to a depth of 0, because the max depth of a new WF is the thickness of the first layer
      }
      else{
        break;
      }
      current = current->next;
      next = current->next;
    }
    if (verbosity.compare("high") == 0) {
      printf("after sending WFs up, during new WF creation process: \n");
      listPrint(*head);
    }


    bool merged_in_non_top_layer = false; //because WF was inserted in the top layer
    lgar_merge_surface_and_TO_wetting_fronts(merged_in_non_top_layer, num_layers, cum_layer_thickness_cm, head);
    lgarto_correct_negative_depths(head);
    if (verbosity.compare("high") == 0) {
      printf("after a surface-TO merge: \n");
      listPrint(*head);
    }

    if (verbosity.compare("high") == 0) {
      printf("before adjust_new_theta: \n");
      listPrint(*head);
    }

    if (listLength_surface(*head)>0){//there is a really rare case where a dry over wet happens before this and the TO WF below is wetter, so there are no surf WFs left
      bottom_boundary_flux_cm_temp += adjust_new_theta(new_wf_num, target_mass, cum_layer_thickness_cm, soil_type, soil_properties, head);
    }
    else {
      new_wf_deleted_flag = true;
    }

    if (verbosity.compare("high") == 0) {
      printf("after adjust_new_theta: \n");
      listPrint(*head);
    }
  }
  /*************************************************************************************/


  /*************************************************************************************/
  //now that mass balance has been corrected, some cleanup is necessary. 
  // in the event that the new WF went deep enough and is dry enough to become a TO WF, it should. It often will.
  // probably no longer necessary due to the fact that listDeleteFront handles such corrections now 
  current = *head;
  next = current->next;
  for (int wf = 1; wf != (listLength(*head)); wf++) {
    if ( (current->is_WF_GW==0) & (next->is_WF_GW==1) & (current->layer_num==next->layer_num) & (current->theta<next->theta) ){
      current->is_WF_GW = 1;
    }
    current = current->next;
    next = current->next;
  }


  // in the event that the new wetting front is saturated and very close in value to the WF below it, it should be deleted
  current = *head;
  next = current->next;
  for (int wf = 1; wf != (listLength(*head)); wf++) {
    if ( (next!=NULL) && (current!=NULL) ){
      if ( (current->is_WF_GW==1) & (next->is_WF_GW==1) & (current->layer_num==next->layer_num) & (next->psi_cm==0) & ( fabs(current->theta - next->theta)<1e-7 ) ){
        current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
      }
    }
    current = next;
    next = current->next;
    if (wf>=listLength(*head)){
      break;
    }
  } 

  lgar_merge_wetting_fronts(head, soil_type, frozen_factor, soil_properties); //in the event that the new WF is deeper than the shallowest surface WF, these have to merge. Should eventually make sure that all merging and layer bdy crossing code run not just once, but run until they are no longer needed 


  if (verbosity.compare("high") == 0) {
    printf("before code that deletes negative WFs: \n");
    listPrint(*head);
  }

  current = *head;
  next = current->next;
  for (int wf=1; wf != listLength(*head); wf++) {
    if (next!=NULL){
      if (current->depth_cm < cum_layer_thickness_cm[current->layer_num - 1]){//There is a rare case when mass conservative wetting front merging will yield a wetting front with a depth above the layer it is in, which can occur when two wetting fronts are extremely close in moisture value. 
        if (verbosity.compare("high") == 0) {
          printf("wf to be deleted: %d \n",wf);
          printf("wf to be deleted based on current->front_num: %d \n",current->front_num);
          printf("depth of that wf: %lf \n",current->depth_cm);
        }
        current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
      }
      current = next;
      next = current->next;
    }
    if (wf >= listLength(*head)){
      break;
    }
  }

  if (verbosity.compare("high") == 0) {
    printf("after WF creation and code that deletes negative WFs: \n");
    listPrint(*head);
  }
  /*************************************************************************************/

  if (new_wf_deleted_flag){
    bottom_boundary_flux_cm_temp += prior_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  }

  return bottom_boundary_flux_cm_temp;
}


// ############################################################################################
/* This routine calculates the "dry depth" of a newly created wetting front in the top soil layer after
   a non-rainy period or a big increase in rainrate on an unsaturated first layer.
   Note: Calculation of the initial depth of a new wetting front in the first layer uses the concept of "dry depth",
   described in the 2015 GARTO paper (Lai et al., An efficient and guaranteed stable numerical method ffor
   continuous modeling of infiltration and redistribution with a shallow dynamic water table). */
   /*In LGARTO, it's pretty often the case that there is precipitation on a model 
   domain with no surface WFs and several GW wetting fronts that are very close to the surface.
   This can cause a new WF that is either too dry, too deep, or both, and a partial solution was 
   to make the dry_depth a fraction of its original value, which also makes sense given that the dry depth as 
   formulated in the 2015 paper mentioned above only applies to the first few moments of a new WF, whereas
   we are using a coarser time step in this model*/
// ############################################################################################
extern double lgar_calc_dry_depth(bool use_closed_form_G, bool TO_enabled, int nint, double timestep_h, double *delta_theta, int *soil_type,
				  double *cum_layer_thickness_cm, double *frozen_factor,
				  struct wetting_front* head, struct soil_properties_ *soil_properties)
{

  // local variables
  struct wetting_front *current;
  double theta1,theta2,theta_e,theta_r;
  double vg_alpha_per_cm,vg_n,vg_m,Ksat_cm_per_h,h_min_cm;
  double tau;
  double Geff;
  double dry_depth;
  int    soil_num;
  int    layer_num;

  current=head;

  if (listLength_surface(head)>0){
    while (current->is_WF_GW==1){
      current = current->next;
    }
  }

  if (current->depth_cm==0.0){
    while (current->depth_cm==0.0){
      current = current->next;
    }
  }

  layer_num  = current->layer_num;
  soil_num   = soil_type[layer_num];

  // copy values of soil properties into shorter variable names to improve readability
  theta_r         = soil_properties[soil_num].theta_r;
  vg_alpha_per_cm = soil_properties[soil_num].vg_alpha_per_cm;
  vg_m            = soil_properties[soil_num].vg_m;
  vg_n            = soil_properties[soil_num].vg_n;
  Ksat_cm_per_h   = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[layer_num];
  h_min_cm        = soil_properties[soil_num].h_min_cm;
  double lambda = soil_properties[soil_num].bc_lambda;
  double bc_psib_cm = soil_properties[soil_num].bc_psib_cm;

  // these are the limits of integration
  theta1   = current->theta;                 // water content of the first (most surficial) existing wetting front
  theta_e  = soil_properties[soil_num].theta_e;
  theta2 = theta_e;

  *delta_theta = theta_e - current->theta;  // return the delta_theta value to the calling function

  tau  = timestep_h * Ksat_cm_per_h/(theta_e-current->theta); //3600

  Geff = calc_Geff(use_closed_form_G, theta1, theta2, theta_e, theta_r, vg_alpha_per_cm, vg_n, vg_m, h_min_cm, Ksat_cm_per_h, nint, lambda, bc_psib_cm); 

  // note that dry depth originally has a factor of 0.5 in front
  dry_depth = 0.5 * (tau + sqrt( tau*tau + 4.0*tau*Geff) );

  //when dry depth greater than layer 1 thickness, set dry depth to layer 1 thickness
  dry_depth = fmin(cum_layer_thickness_cm[layer_num], dry_depth);

  if (TO_enabled){
    dry_depth = dry_depth * 0.2; //in LGARTO, solution quality is sensitive to dry depth, and the dry depth concept as explained in both the 1997 GAR and 2023 LGAR papers probably is only applicable for small time steps. 
  //It is easy to demonstrate that an unaltered dry depth initializes wetting fronts in an unrealistic way in some cases, effectively assuming an enormous dzdt value during the time step for which the WF was created.
  //LGAR (not LGARTO) simulations seem to be mostly insensitive to reducing the dry depth by a factor of 0.1-0.2, while this range of factors massively improves LGARTO simulations and makes wetting fronts have more reasonable initial depths. 
  //Some text on this topic should be in the Discussion section of the LGARTO paper.
  }
  
  return dry_depth;

}

// ###########################################################################
/* function to calculate the amount of soil moisture (total mass of water)
   in the profile (cm) */
// ###########################################################################
double lgar_calc_mass_bal(double *cum_layer_thickness, struct wetting_front* head)
{

  struct wetting_front* current;
  struct wetting_front* next;     // Beware, might get a little confusing, because there exists a next->next

  double sum=0.0;
  double base_depth;
  int layer;

  if(head == NULL) return 0.0;

  current=head;

  do
    {
      layer=current->layer_num;
      base_depth=cum_layer_thickness[layer-1];   // note cum_layer_thickness[0]=0.0;

      if(current->next != NULL) {            // this is not the last entry in the list
	next=current->next;
	if(next->layer_num == current->layer_num){
	  sum += (current->depth_cm - base_depth) * (current->theta - next->theta); // note no need for fabs() here otherwise we get more mass for the case dry-over-wet front within a layer
  }
	else{
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

// ############################################################################################
/* The module reads the soil parameters.
   Open file to read in the van Genuchten parameters for standard soil types*/
// ############################################################################################
extern int lgar_read_vG_param_file(char const* vG_param_file_name, int num_soil_types, double wilting_point_psi_cm,
				    struct soil_properties_ *soil_properties)
{

  if (verbosity.compare("high") == 0) {
    std::cerr<<"Reading van Genuchten parameters files...\n";
  }

  // local vars
  FILE *in_vG_params_fptr = NULL;
  char jstr[256];
  char soil_name[30];
  // bool error;
  int length;
  int num_soils_in_file = 0;             // soil counter
  int soil = 1;
  double theta_e,theta_r,vg_n,vg_m,vg_alpha_per_cm,Ksat_cm_per_h;  // shorthand variable names
  double m,p,lambda;

  // open the file
  if((in_vG_params_fptr=fopen(vG_param_file_name,"r"))==NULL) {
    printf("Can't open input file named %s. Program stopped.\n",vG_param_file_name); exit(-1);
  }

  fgets(jstr,255,in_vG_params_fptr);   // read the header line and ignore

  while (fgets(jstr,255,in_vG_params_fptr) != NULL) {

    sscanf(jstr,"%s %lf %lf %lf %lf %lf",soil_name,&theta_r,&theta_e,&vg_alpha_per_cm,&vg_n,&Ksat_cm_per_h);
    length=strlen(soil_name);

    if(length>MAX_SOIL_NAME_CHARS) {
      printf("While reading vG soil parameter file: %s, soil name longer than allowed.  Increase MAX_SOIL_NAME_CHARS\n",
	     vG_param_file_name);
      printf("Program stopped.\n");
      exit(0);
    }

    strcpy(soil_properties[soil].soil_name,soil_name);
    soil_properties[soil].theta_r         = theta_r;
    soil_properties[soil].theta_e         = theta_e;
    soil_properties[soil].vg_alpha_per_cm = vg_alpha_per_cm; // cm^(-1)
    soil_properties[soil].vg_n            = vg_n;
    soil_properties[soil].vg_m            = 1-1/vg_n;
    vg_m = soil_properties[soil].vg_m;
    soil_properties[soil].Ksat_cm_per_h   = Ksat_cm_per_h;
    soil_properties[soil].theta_wp = calc_theta_from_h(wilting_point_psi_cm, vg_alpha_per_cm,
						       vg_m, vg_n, theta_e, theta_r);

    // Given van Genuchten parameters calculate estimates of Brooks & Corey bc_lambda and bc_psib
    if (1.0 < vg_n) {
      m = 1.0 - 1.0 / vg_n;
      p = 1.0 + 2.0 / m;
      soil_properties[soil].bc_lambda  = 2.0 / (p - 3.0);
      soil_properties[soil].bc_psib_cm = (p + 3.0) * (147.8 + 8.1 * p + 0.092 * p * p) /
	(2.0 * vg_alpha_per_cm * p * (p - 1.0) * (55.6 + 7.4 * p + p * p));
      assert(0.0 < soil_properties[soil].bc_psib_cm);
    }
    else {
      fprintf(stderr, "ERROR: van Genuchten parameter n must be greater than 1\n");
      //error = TRUE;  // TODO FIXME - what todo in this ccase?
    }

    /* this is the effective capillary drive after */
    /* Morel-Seytoux et al. (1996) eqn. 13 or 15 */
    /* psi should not be less than this value.  */
    lambda=soil_properties[soil].bc_lambda;
    soil_properties[soil].h_min_cm = soil_properties[soil].bc_psib_cm*(2.0+3.0/lambda)/(1.0+3.0/lambda);
    num_soils_in_file++;
    soil++;

    // the soil_properties array has allocation of num_soil_types, so break if condition is satisfied
    if (num_soil_types == num_soils_in_file)
      break;
  }

  fclose(in_vG_params_fptr);      // close the file, we're done with it

  return num_soils_in_file;
}

// ############################################################################################
/* code to calculate velocity of fronts
   equations with full description are provided in the lgar paper (currently under review) */
// ############################################################################################
extern void lgar_dzdt_calc(bool use_closed_form_G, int nint, double timestep_h, double h_p, int *soil_type, double *cum_layer_thickness_cm,
			   double *frozen_factor, struct wetting_front* head, struct soil_properties_ *soil_properties, int num_layers)
{
  if (verbosity.compare("high") == 0) {
    std::cerr<<"Calculating dz/dt .... \n";
    printf("list at start of lgar_dzdt_calc: \n");
    listPrint(head);
  }

  struct wetting_front* current;
  struct wetting_front* next;
  struct wetting_front* previous;

  double vg_alpha_per_cm,vg_n,vg_m,Ksat_cm_per_h,theta_e,theta_r;  // local variables to make things clearer
  double delta_theta;
  double Geff;
  double depth_cm;    // the absolute depth down to a wetting front from the surface
  double h_min_cm;
  double K_cm_per_h;  // unsaturated hydraulic conductivity K(theta) at the RHS of the current wetting front (cm/h)
  double theta1, theta2;  // limits of integration on Geff from theta1 to theta2
  double bottom_sum;  // store a running sum of L_n/K(theta_n) n increasing from 1 to N-1, as we go down in layers N
  double dzdt;
  int    soil_num, layer_num;


  if(head == NULL) {
    stringstream errMsg;
    errMsg << "lgar derivative function called for empty list (no wetting front exists) \n";
    throw runtime_error(errMsg.str());
  }

  current = head;

  do {  // loop through the wetting fronts
    if (verbosity.compare("high") == 0) {
      printf("calculating dzdt for front: %d \n",current->front_num);
    }
    dzdt = 0.0;

    if (current->is_WF_GW){ //here we calculate dzdt for wetting fronts that are in contact with the groundwater 
      if(current->to_bottom == TRUE){
        dzdt = 0.0; //could probably be continue but this makes sure that dzdt = 0.0
      }
      else{
        next = current->next; 
        previous = listFindFront(current->front_num-1,head,NULL);
        layer_num    = current->layer_num; 
        if (num_layers == current->layer_num){
          //Because the WF is in the lowest layer, the single layer form of dZ/dt for the TO WF can be used. 
          //note that TO WFs extend from the bottom of the model domain up, whereas surface WFs extend from the top of the model domain down. Therefore, a TO WF spans just 1 layer if it is only in the bottom layer, but a TO WF that is in more than the bottom layer spans more than 1 layer. 
          double D = cum_layer_thickness_cm[num_layers]; //This is the depth to GW. When we update LGARTO to have variable groundwater depth, this can vary; for now it doesn't 

          double avoid_div_by_zero_factor = 1.E-9;
          if ((next->theta - current->theta) == 1.E-9){
            avoid_div_by_zero_factor = 1.E-10;
          }
          dzdt = (next->K_cm_per_h - current->K_cm_per_h) / (next->theta - current->theta + avoid_div_by_zero_factor) * (1 - (next->psi_cm+1e-6) / (D - current->depth_cm + avoid_div_by_zero_factor)); //the small numbers 1e-6 and 1e-9 are designed to prevent psi = 0 and division by 0
          
          if (verbosity.compare("high") == 0) {
              printf("initial dzdt calc for GW WF: %lf \n ",dzdt);
              printf("current->front_num: %d \n ",current->front_num);
              printf(" \n ");
            }

            if ( (dzdt*timestep_h+current->depth_cm) > D ){
              dzdt = -1*(current->depth_cm + D-1e-6)*(1/timestep_h); //you don't technically need this, but if you don't have it, then a TO WF with a psi value of less than 1e-6 (really, for the psi = 0 WF) will oscillate around its hydrostatic position somewhat intensely 
              if (verbosity.compare("high") == 0) {
                printf("dzdt for GW WF, oscillation prevention case: %lf \n ",dzdt);
                printf("current->front_num: %d \n ",current->front_num);
                printf(" \n ");
              }
            }

          if ( (current->is_WF_GW==TRUE) & (next->is_WF_GW==FALSE) ){ //In the event that there are TO WFs above surface WFs, the deepest such TO WF can't use next because it has to skip over all the surface WFs to get to the relevant TO WF for dZ/dt calculation
            struct wetting_front* highest_TO_WF_below_surface_WFs; 
            highest_TO_WF_below_surface_WFs = next;
            while (highest_TO_WF_below_surface_WFs->is_WF_GW!=1){
              highest_TO_WF_below_surface_WFs = highest_TO_WF_below_surface_WFs->next;
            }
            while ( (highest_TO_WF_below_surface_WFs->to_bottom==TRUE) && (highest_TO_WF_below_surface_WFs->next!=NULL) ){// idea is that highest_TO_WF_below_surface_WFs can not be a boundary (to_bottom) WF
              highest_TO_WF_below_surface_WFs = highest_TO_WF_below_surface_WFs->next;
            }
              dzdt = (highest_TO_WF_below_surface_WFs->K_cm_per_h - current->K_cm_per_h) / (highest_TO_WF_below_surface_WFs->theta - current->theta + 1e-9) * (1 - (highest_TO_WF_below_surface_WFs->psi_cm+1e-6) / (D - current->depth_cm + 1e-9));
            if (verbosity.compare("high") == 0) {
              printf("dzdt for GW WF, deepest TO WF above surface WFs case: %lf \n ",dzdt);
              printf("current->front_num: %d \n ",current->front_num);
            }
          }

          if (isnan(dzdt)){//this should be extremely rare and apparently is no longer triggered after new code that prevents division by 0, but keeping in to be safe
            dzdt = 0.0;
            if (verbosity.compare("high") == 0) {
              printf("dzdt case where isnan(dzdt) is true so dzdt is set to 0: %lf \n ",dzdt);
              printf("current->front_num: %d \n ",current->front_num);
              printf(" \n ");
            }
          }
        }
        
        else{
          //Because the WF is above the lowest layer, the multilayer form of dZ/dt for the TO WF must be used. 
          double K_composite = 0; 
          double K_composite_left = 0;
          double denominator = 0;
          double denominator_left = 0;
          double D = cum_layer_thickness_cm[num_layers]; //This is the depth to GW, or the depth of the deepest TO WF. When we update LGARTO to have variable groundwater depth, this can vary; for now it doesn't

          if ( (current->is_WF_GW==1) & (next->is_WF_GW==0) ){ //In the event that there are TO WFs above surface WFs, the deepest such TO WF can't use next because it has to skip over all the surface WFs to get to the relevant TO WF for dZ/dt calculation
            while (next->is_WF_GW!=1){
              next = next->next;
            }
            while ( (next->to_bottom==TRUE) && (next->next!=NULL) ){// the next unique TO WF in terms of psi is necessary for next. If next is a boundary wetting front, this won't be the case, at least if the WF directly after current is a surface WF.
              next = next->next;
            }
          }

          for (int k = num_layers; k > (current->layer_num-1); k--) {
            int soil_num_loc = soil_type[k]; // _loc denotes the soil_num is local to this loop
            double theta_prev_loc = calc_theta_from_h(next->psi_cm, soil_properties[soil_num_loc].vg_alpha_per_cm,
                        soil_properties[soil_num_loc].vg_m,
                        soil_properties[soil_num_loc].vg_n,soil_properties[soil_num_loc].theta_e,
                        soil_properties[soil_num_loc].theta_r);

            double Se_prev_loc = calc_Se_from_theta(theta_prev_loc,soil_properties[soil_num_loc].theta_e,soil_properties[soil_num_loc].theta_r);

            double K_cm_per_h_prev_loc = calc_K_from_Se(Se_prev_loc,soil_properties[soil_num_loc].Ksat_cm_per_h * frozen_factor[k],
                          soil_properties[soil_num_loc].vg_m);

            if (k==current->layer_num){
              denominator += (cum_layer_thickness_cm[k] - current->depth_cm)/ K_cm_per_h_prev_loc;
            }
            else{
              denominator += (cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k-1])/ K_cm_per_h_prev_loc;
            }
          }

          double temp_psi = current->psi_cm;
          double temp_theta = current->theta;
          bool identical_thetas = false;
          if (current->psi_cm == next->psi_cm){
            temp_psi = temp_psi + 0.1;
            identical_thetas = true;
            int soil_num_loc = soil_type[current->layer_num]; // _loc denotes the soil_num is local to this loop
            temp_theta = calc_theta_from_h(temp_psi, soil_properties[soil_num_loc].vg_alpha_per_cm,
                        soil_properties[soil_num_loc].vg_m,
                        soil_properties[soil_num_loc].vg_n,soil_properties[soil_num_loc].theta_e,
                        soil_properties[soil_num_loc].theta_r);
          }
          for (int k = num_layers; k > (current->layer_num-1); k--) {
            int soil_num_loc = soil_type[k]; // _loc denotes the soil_num is local to this loop
            double theta_prev_loc = calc_theta_from_h(temp_psi, soil_properties[soil_num_loc].vg_alpha_per_cm,
                        soil_properties[soil_num_loc].vg_m,
                        soil_properties[soil_num_loc].vg_n,soil_properties[soil_num_loc].theta_e,
                        soil_properties[soil_num_loc].theta_r);

            double Se_prev_loc = calc_Se_from_theta(theta_prev_loc,soil_properties[soil_num_loc].theta_e,soil_properties[soil_num_loc].theta_r);

            double K_cm_per_h_prev_loc = calc_K_from_Se(Se_prev_loc,soil_properties[soil_num_loc].Ksat_cm_per_h * frozen_factor[k],
                          soil_properties[soil_num_loc].vg_m);
            if (k==layer_num){
              if (previous!=NULL){
                denominator_left += (cum_layer_thickness_cm[k] - current->depth_cm)/ K_cm_per_h_prev_loc;
              }
              else{
                denominator_left += (cum_layer_thickness_cm[k])/ K_cm_per_h_prev_loc;
              }
            }
            else{
              denominator_left += (cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k-1])/ K_cm_per_h_prev_loc;
            }
            
          }

          K_composite = (D - current->depth_cm) / (denominator);
          if (previous==NULL){
            K_composite_left = D / (denominator_left);
          }
          else{
            K_composite_left = (D - current->depth_cm) / (denominator_left);
          }

          if (current->layer_num == next->layer_num){
            if (!identical_thetas){
              dzdt = (K_composite - K_composite_left) / (next->theta - current->theta + 1e-9) * (1 - next->psi_cm / (D - current->depth_cm));
            }
            else{
              //At some point in LGARTO's development, there were cases where adjacent WFs would have identical theta values ... I think this has been resolved, but just in case it hasn't this is good for redundancy
              dzdt = (K_composite - K_composite_left) / (next->theta - temp_theta + 1e-9) * (1 - next->psi_cm / (D - current->depth_cm)); 
            }
          }
          else{//this is for the case where there is a block of surface WFs in between TO WFs that extend to the surface and TO WFs that are below surface WFs, in this case you might have next->layer_num!=current->layer_num
              double theta_temp;
              int soil_num_loc = soil_type[current->layer_num]; // _loc denotes the soil_num is local to this loop
              theta_temp = calc_theta_from_h(next->psi_cm, soil_properties[soil_num_loc].vg_alpha_per_cm,
                          soil_properties[soil_num_loc].vg_m,
                          soil_properties[soil_num_loc].vg_n,soil_properties[soil_num_loc].theta_e,
                          soil_properties[soil_num_loc].theta_r);
              if (!identical_thetas){
                dzdt = (K_composite - K_composite_left) / (theta_temp - current->theta + 1e-9) * (1 - next->psi_cm / (D - current->depth_cm));
              }
              else {
                dzdt = (K_composite - K_composite_left) / (next->theta - theta_temp + 1e-9) * (1 - next->psi_cm / (D - current->depth_cm));
              }
          }

        }
      }

      if ( (current->depth_cm + dzdt*timestep_h)<0.0 ){//prevent WFs from going above the soil surface via dzdt         
        dzdt = current->depth_cm/timestep_h; //so this is just one of many ways to solve the problem ... should be somewhat uncommon
      } 

      if (current->next!=NULL){
        if (current->next->theta==current->theta){
          dzdt = 0.0;
        }
      }

    }

    else{//the above code was dzdt for GW WFs. This code is dzdt for surface WFs.

      // copy structure elements into shorter variables names to increase readability
      // WETTING FRONT PROPERTIES
      layer_num    = current->layer_num;    // what layer the front is in
      K_cm_per_h   = current->K_cm_per_h;   // K(theta)

      if (K_cm_per_h < 0) {
        printf("K is negative (layer_num, wf_num, K): %d %d %lf \n", layer_num, current->front_num, K_cm_per_h);
        listPrint(head);
        printf("Is your n value very close to 1? Very small n values can cause K to become 0. \n");
        //The parameter n must physically attain a value greater than 1. However, when n is small, and apparently less than 1.02, sometimes n can make K evaluate to 0, for larger values of psi.
        //So, checking for K_cm_per_h <= 0 has been replaced by checking if K_cm_per_h is negative. K_cm_per_h should never be negative (although perhaps machine precision could make this occur, although we haven't seen it yet), but mathematically can be 0 in some rare cases. 
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
      Ksat_cm_per_h   = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num];
      double lambda = soil_properties[soil_num].bc_lambda;
      double bc_psib_cm = soil_properties[soil_num].bc_psib_cm;

      next = current->next;    // the next element in the linked list
      if (next == NULL) break; // we're done calculating dZ/dt's because we're at the end of the list

      theta1 = next->theta;
      theta2 = current->theta;


      bottom_sum = 0.0;  // needed ffor multi-layered dz/dt equation.  Equal to sum from n=1 to N-1 of (L_n/K_n(theta_n))

      if(current->to_bottom == TRUE) {
        if(layer_num > 1)
    current->dzdt_cm_per_h = 0.0;
        else
    current->dzdt_cm_per_h = 0.0;

        current = current->next;  // point to the next link
        continue;                 // go to next front, this one fully penetrates the layer
      }
      else if(layer_num > 1) {
        bottom_sum += (current->depth_cm-cum_layer_thickness_cm[layer_num-1])/K_cm_per_h;
      }

      if(theta1 > theta2) {
        printf("Calculating dzdt : theta1 > theta2 = (%lf, %lf) ... aborting \n", theta1, theta2);  // this should never happen
        listPrint(head);
        exit(0);
      }

      Geff = calc_Geff(use_closed_form_G, theta1, theta2, theta_e, theta_r, vg_alpha_per_cm, vg_n, vg_m, h_min_cm, Ksat_cm_per_h, nint, lambda, bc_psib_cm); 
      delta_theta = current->theta - next->theta;

      if(current->layer_num == 1) { // this front is in the upper layer
        if (delta_theta > 0){
    dzdt = 1.0/delta_theta*(Ksat_cm_per_h*(Geff+h_p)/current->depth_cm+current->K_cm_per_h);}
        else{
    dzdt = 0.0;}
      }
      else {  // we are in the second or greater layer
        double denominator = bottom_sum;

        for (int k = 1; k < layer_num; k++) {
    int soil_num_loc = soil_type[layer_num-k]; // _loc denotes the soil_num is local to this loop
    double theta_prev_loc = calc_theta_from_h(current->psi_cm, soil_properties[soil_num_loc].vg_alpha_per_cm,
                soil_properties[soil_num_loc].vg_m,
                soil_properties[soil_num_loc].vg_n,soil_properties[soil_num_loc].theta_e,
                soil_properties[soil_num_loc].theta_r);


    double Se_prev_loc = calc_Se_from_theta(theta_prev_loc,soil_properties[soil_num_loc].theta_e,soil_properties[soil_num_loc].theta_r);

    double K_cm_per_h_prev_loc = calc_K_from_Se(Se_prev_loc,soil_properties[soil_num_loc].Ksat_cm_per_h * frozen_factor[layer_num-k],
                  soil_properties[soil_num_loc].vg_m);

    denominator += (cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k-1])/ K_cm_per_h_prev_loc;

        }

        double numerator = depth_cm;// + (Geff +h_p)* Ksat_cm_per_h / K_cm_per_h;

        if (delta_theta > 0)
    dzdt = (1.0/delta_theta) * ( (numerator / denominator) + Ksat_cm_per_h*(Geff+h_p)/depth_cm );
        else
    dzdt = 0.0;
      }
    }

    if ((dzdt == 0.0) && (current->is_WF_GW==FALSE) && (current->to_bottom==FALSE)){
      //in lgar_move, we have: "if (current->dzdt_cm_per_h == 0.0 && current->to_bottom == FALSE) { // a new front was just created, so don't update it."
      //the issue here is that when theta approaches theta_r, then dzdt can in some cases numerically evaluate to 0, even if the wetting front has to_bottom==FALSE.
      //so, there are cases where a WF should be moving very slowly, but not being completely still. 
      dzdt = 1e-9;
    }

    if (dzdt>1E4){//insanity check. Keep in mind that this can be quite relevant, as adjacent WFs that are very close in theta can yield enormous, unrealistic dzdt values
      dzdt = 1E4;
    }

    if (dzdt<-1E4){//insanity check
      dzdt = -1E4;
    }

    double largest_K_s = 0.0;
    for (int ii=1; ii<=num_layers; ii++) {
      int soil_num = soil_type[ii];
      double temp_K_s = soil_properties[soil_num].Ksat_cm_per_h;
      largest_K_s = fmax(temp_K_s, largest_K_s);
    }

    if (dzdt>100*largest_K_s){//insanity check; was 1E4 but now addtionally defining based on K_s
      dzdt = 100*largest_K_s;
    }

    if (dzdt<-100*largest_K_s){//insanity check; was -1E4 
      dzdt = -100*largest_K_s;
    }
  
    current->dzdt_cm_per_h = dzdt;

    if (verbosity.compare("high") == 0){
      printf("dzdt: %.17lf \n", dzdt);
      printf("just finished dzdt for wf: %d \n", current->front_num);
    }

    current = current->next;  // point to the next link

  } while(current != NULL );   // putting conditional at end of do looop makes sure it executes at least once

  if (verbosity.compare("high") == 0) {
    printf("end of dzdt calculations: \n");
    listPrint(head);
  }

}

// ############################################################################################
/* The function does mass balance for a wetting front to get an updated theta.
   The head (psi) value is iteratively altered until the error between prior mass and new mass
   is within a tolerance. */
// ############################################################################################
extern double lgar_theta_mass_balance(int layer_num, int soil_num, double psi_cm, double new_mass,
				      double prior_mass, double *AET_demand_cm, double *delta_theta, double *delta_thickness,
				      int *soil_type, struct soil_properties_ *soil_properties)
{

  double psi_cm_loc = psi_cm; // location psi
  double delta_mass = fabs(new_mass - prior_mass); // mass different between the new and prior
  double tolerance = 1e-12;

  double factor = fmax(1,psi_cm/100);//was 1.0 previously. This code is far faster and seems to avoid loops with >10000 iterations. Low van Genuchten n values can cause this
  bool switched = false; // flag that determines capillary head to be incremented or decremented

  double theta             = 0; // this will be updated and returned
  double psi_cm_loc_prev   = psi_cm_loc;
  double delta_mass_prev   = delta_mass;
  int count_no_mass_change = 0;
  int break_no_mass_change = 5;
  bool wanted_to_saturate_flag = FALSE;
  
  // check if the difference is less than the tolerance
  if (delta_mass <= tolerance) {
    theta = calc_theta_from_h(psi_cm_loc, soil_properties[soil_num].vg_alpha_per_cm,
			      soil_properties[soil_num].vg_m, soil_properties[soil_num].vg_n,
			      soil_properties[soil_num].theta_e,soil_properties[soil_num].theta_r);
    return theta;
  }

  // the loop increments/decrements the capillary head until mass difference between
  // the new and prior is within the tolerance
  int iter = 0;
  bool iter_aug_flag = FALSE;

  while (delta_mass > tolerance) {
    iter++;

    if (iter>1000 && iter_aug_flag==FALSE){
      factor = factor*100;
      iter_aug_flag = TRUE;
    }

    if (new_mass > prior_mass) {
      psi_cm_loc += 0.1 * factor;
      switched = false;
    }
    else {
      if (!switched) {
	switched = true;
	factor = factor * 0.1;
      }
      
      psi_cm_loc_prev = psi_cm_loc;
      psi_cm_loc -= 0.1 * factor;

      if (psi_cm_loc<0.0){
        wanted_to_saturate_flag = TRUE;
      }
      
      if (psi_cm_loc < 0 && psi_cm_loc_prev != 0) {
	/* this is for the extremely rare case when iterative psi_cm_loc calculation temporarily
	   yields a negative value and the actual answer for psi_cm_loc is nonzero. For example
	   when a completely saturated wetting front with a tiny amount of ET should yield a resulting
	   theta that is slightly below saturation. */
        psi_cm_loc = psi_cm_loc_prev * 0.1;
      }
      
    }

    double theta_layer;
    double mass_layers= 0.0;

    theta = calc_theta_from_h(psi_cm_loc, soil_properties[soil_num].vg_alpha_per_cm, soil_properties[soil_num].vg_m,
			      soil_properties[soil_num].vg_n,soil_properties[soil_num].theta_e,
			      soil_properties[soil_num].theta_r);

    mass_layers += delta_thickness[layer_num] * (theta - delta_theta[layer_num]);

    for (int k=1; k<layer_num; k++) {
      int soil_num_loc =  soil_type[k]; // _loc denotes the variable is local to the loop

      theta_layer = calc_theta_from_h(psi_cm_loc, soil_properties[soil_num_loc].vg_alpha_per_cm,
				      soil_properties[soil_num_loc].vg_m, soil_properties[soil_num_loc].vg_n,
				      soil_properties[soil_num_loc].theta_e, soil_properties[soil_num_loc].theta_r);

      mass_layers += delta_thickness[k] * (theta_layer - delta_theta[k]);
    }

    new_mass = mass_layers;
    delta_mass = fabs(new_mass - prior_mass);

    // stop the loop if the error between the current and previous psi is less than 10^-15
    // 1. enough accuracy, 2. the algorithm can't improve the error further,
    // 3. avoid infinite loop, 4. handles case where theta is very close to theta_r and convergence might be possible but would be extremely slow
    // 5. handles a corner case when prior mass is tiny (e.g., <1.E-5)
    
    if (fabs(psi_cm_loc - psi_cm_loc_prev) < 1E-15 && factor < 1E-13) {
      break;
    }

    // another condition to avoid infinite loop when the error does not improve
    if (fabs(delta_mass - delta_mass_prev) < 1E-15)
      count_no_mass_change++;
    else
      count_no_mass_change = 0;

    // break the loop if the mass does not change in the five consecutive iterations.
    if (count_no_mass_change == break_no_mass_change){
      break;
    }
    
    if (psi_cm_loc > 1e7){//there are rare cases where theta is very close to theta_r, and delta_mass - delta_mass_prev will change extremely slowly. Convergence might be possible but the model will take hours to converge rather than seconds. 
    //an alternative solution was to change the threshold in if (fabs(delta_mass - delta_mass_prev) < 1e-15) to 1e-11, but that solution is somewhat slow. 
      break;
    }

    // -ve pressure will return NAN, so terminate the loop if previous psi is way small and current psi is zero
    // the wetting front is almost saturated
    if (psi_cm_loc <= 0 && psi_cm_loc_prev < 0) {
      break;
    }

    delta_mass_prev = delta_mass;

  }

  //There is a rare case where mass balance closure would require that theta<theta_r. 
  //However, the above loop can never increase psi to the point where theta<theta_r, because theta must always be between theta_r and theta_r, because of the van Genuchten model (calc_theta_from_h).
  //If we get to the case where theta<theta_r would be necessary for mass balance closure, then the above loop will break before delta_mass <= tolerance.
  //In this rare case, the remaining mass balance error is put into AET. 
  if ((delta_mass > tolerance) && (!wanted_to_saturate_flag)){//the second condition is necessary because count_no_mass_change == break_no_mass_change in the loop above will trigger when the model approaches saturation; in this event the extra water should go into runoff (handled eslewhere), because the soil saturates, rather than AET
    *AET_demand_cm = *AET_demand_cm - fabs(delta_mass - tolerance);
  }

  /////freeing dynamically allocated memory is now done outside of the function call
  // free(delta_thickness);
  // free(delta_theta);

  return theta;

}


extern void lgarto_resolve_TO_WF_between_surf_WFs(int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head)
{
  if (verbosity.compare("high") == 0) {
    printf("states before lgarto_resolve_TO_WF_between_surf_WFs: \n");
    listPrint(*head);
  }
  struct wetting_front *current;
  current = *head;
  double depth_of_top_most_mobile_TO_WF = 0.0;
  for (int wf = 1; wf != (listLength(*head)); wf++){
    if ( (current->is_WF_GW==1) && (current->depth_cm>0.0) && (current->to_bottom==FALSE) ){
      depth_of_top_most_mobile_TO_WF = current->depth_cm;
      break;
    }
    current = current->next;
  }

  bool surface_WFs_deeper_than_top_mobile_TO_WF = GW_fronts_among_surf_WFs(*head);

  while (surface_WFs_deeper_than_top_mobile_TO_WF && depth_of_top_most_mobile_TO_WF>0.0){//do the thing here 2
    struct wetting_front *current;
    struct wetting_front *next;
    current = *head; 
    next = current->next;

    if (verbosity.compare("high") == 0) {
      printf("Inside merging wetting fronts, case where a TO WF is between surface WFs... \n");
    }
    for (int wf=1; wf != listLength(*head); wf++) {

      if ( (current->is_WF_GW==FALSE) && (next->is_WF_GW==TRUE) ){

      next = current->next; 

      // case : wetting front passing another wetting front within a layer
      /**********************************************************/
      // 'current->depth_cm > next->depth_cm' ensures that merging is needed
      // 'current->layer_num == next->layer_num' ensures wetting fronts are in the same layer
      // '!next->to_bottom' ensures that the next wetting front is not the deepest wetting front in the layer
      
        if ( (current->layer_num == next->layer_num) && next->to_bottom==FALSE) {

          next = listDeleteFront(next->front_num, head, soil_type, soil_properties);
          surface_WFs_deeper_than_top_mobile_TO_WF = GW_fronts_among_surf_WFs(*head);
          break;

        }
      }
    current = current->next;
    next = current->next;
    }
  }
  if (verbosity.compare("high") == 0) {
    printf("states after lgarto_resolve_TO_WF_between_surf_WFs: \n");
    listPrint(*head);
  }
}


extern void lgarto_ensure_rooting_zone_population(double rzd, double PET_timestep_cm, int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head){
  if (verbosity.compare("high") == 0) {
    printf("states before lgarto_ensure_rooting_zone_population: \n");
    listPrint(*head);
  }
  int num_TO_WFs_in_rz = listLength_TO_WFs_in_rz(rzd, *head); 
  struct wetting_front *current;
  current = *head; 

  while (( (num_TO_WFs_in_rz - listLength_surface(*head))<4) && (PET_timestep_cm>0.0) && (listLength_surface(*head)==0) && current->psi_cm<1.E6 ){//the 1.E-6 ensures that we don't get absurdly dry WFs towards surface 
  // while (( (num_TO_WFs_in_rz - listLength_surface(*head))<15) && current->psi_cm<1.E6 ){//the 1.E-6 ensures that we don't get absurdly dry WFs towards surface 
    //high accuracy version: use the version where new wetting fronts are created regardless of PET
    current = *head; 
    // double new_psi = current->psi_cm + 30.0; //this could be a parameter ... or perhas calculated as some fraction of the thickness of the model domain. 
    double new_psi = current->psi_cm + 0.25*rzd; //0.25 or 0.05

   /*The idea is that merging of TO WFs with other TO WFs or with surface WFs at various points can make it so that there are no WFs left in the rooting zone except for the one that extends to the soil surface.
    More commonly however, it can simply be the case that periods of relatively low precipitation and high AET can move all TO WFs in the rooting zone out of the rooting zone.
    In these somewhat common cases, no AET would be extracted from the rooting zone. So a new WF must be created in order to allow for AET from TO WFs in the rooting zone.
    This also has the advantage of allowing for soil that is drier than its hydrostatic moisture would be, which we can expect during periods of intense AET and no precip.
    The actual value of new_psi can be tuned -- perhaps the amount that should be added could be something like 10% of the domain depth. Still, the general question of 
    what number of TO WFs is optimal is a borader one (for example depending on soil type, you might want to make new_psi either much greater or smaller than it is here). 
    */

    int soil_num = soil_type[current->layer_num];

    double theta_e   = soil_properties[soil_num].theta_e;
    double theta_r   = soil_properties[soil_num].theta_r;
    double vg_a      = soil_properties[soil_num].vg_alpha_per_cm;
    double vg_m      = soil_properties[soil_num].vg_m;
    double vg_n      = soil_properties[soil_num].vg_n;

    // current psi with van Genuchten properties of the next layer to get new theta
    double theta_new = calc_theta_from_h(new_psi, vg_a, vg_m, vg_n, theta_e, theta_r);

    listInsertFirst(0.0, theta_new, 1, 1, FALSE, head, true);
    // abort();
    current = *head; 
    current->psi_cm = new_psi;

    num_TO_WFs_in_rz = listLength_TO_WFs_in_rz(rzd, *head);

  }
  
  if (verbosity.compare("high") == 0) {
    printf("states after lgarto_ensure_rooting_zone_population: \n");
    listPrint(*head);
  }
}


extern double lgarto_extract_TO_GW_flux_from_surface_WFs(double *bottom_boundary_flux_above_surface_WFs_cm, double bottom_boundary_flux_cm, double *AET_demand_cm,
                                                         double* cum_layer_thickness_cm, int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head){
  //Previously, bottom_boundary_flux_above_surface_WFs_cm was calculated. This mass must be subtracted from the rightmost surface WF.

  if (*bottom_boundary_flux_above_surface_WFs_cm!=0.0){
    double mass_at_start_of_extraction = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
    double bottom_bdy_flux_at_start_of_extraction = bottom_boundary_flux_cm;
    int wf_free_drainage_demand = wetting_front_free_drainage(*head);

    if (verbosity.compare("high") == 0) {
      printf("before lgarto_extract_TO_GW_flux_from_surface_WFs \n");
      listPrint(*head);
      printf("associated mass: %lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
      printf("bottom_boundary_flux_above_surface_WFs_cm: %lf \n", *bottom_boundary_flux_above_surface_WFs_cm);
    }

    struct wetting_front *current = listFindFront(wf_free_drainage_demand, *head, NULL);
    int wf_from_which_to_extract = current->front_num;
    double original_theta = current->theta;
    double original_psi_cm = current->psi_cm;
    struct wetting_front *next = current->next;

    if (listLength_surface(*head)>0){
      if (current->layer_num == 1){ 
        double theta_reduction = *bottom_boundary_flux_above_surface_WFs_cm / current->depth_cm;
        double theta_new = current->theta - theta_reduction;
        current->theta = fmax(soil_properties[soil_type[current->layer_num]].theta_r + 1.E-12, fmin(theta_new, soil_properties[soil_type[current->layer_num]].theta_e));
        
        int soil_num_k  = soil_type[current->layer_num];
        double theta_e_k = soil_properties[soil_num_k].theta_e;
        double theta_r_k = soil_properties[soil_num_k].theta_r;
        double vg_a_k    = soil_properties[soil_num_k].vg_alpha_per_cm;
        double vg_m_k    = soil_properties[soil_num_k].vg_m;
        double vg_n_k    = soil_properties[soil_num_k].vg_n;
        double Se = calc_Se_from_theta(current->theta,theta_e_k,theta_r_k);
        current->psi_cm = calc_h_from_Se(Se, vg_a_k, vg_m_k, vg_n_k); 

        if (current->theta==theta_r_k+1.E-12){
          bottom_boundary_flux_cm -= *bottom_boundary_flux_above_surface_WFs_cm;
          bottom_boundary_flux_cm += (mass_at_start_of_extraction - lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
          *bottom_boundary_flux_above_surface_WFs_cm = (mass_at_start_of_extraction - lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
        }

      }
      else{
        int layer_num = current->layer_num;
        int soil_num  = soil_type[layer_num];
        double vg_a_k, vg_m_k, vg_n_k;
        double theta_e_k, theta_r_k;

        double *delta_thetas    = (double *)malloc(sizeof(double)*(layer_num+1));
        double *delta_thickness = (double *)malloc(sizeof(double)*(layer_num+1));

        double psi_cm = current->psi_cm;
        double psi_cm_below = next->psi_cm;

        double new_mass = (current->depth_cm - cum_layer_thickness_cm[layer_num-1]) * (current->theta - next->theta);
        double prior_mass = new_mass; 

        // compute mass in the layers above the current wetting front
        // use the psi of the current wetting front and van Genuchten parameters of
        // the respective layers to get the total mass above the current wetting front
        for (int k=1; k<layer_num; k++) {
          int soil_num_k  = soil_type[k];
          theta_e_k = soil_properties[soil_num_k].theta_e;
          theta_r_k = soil_properties[soil_num_k].theta_r;
          vg_a_k    = soil_properties[soil_num_k].vg_alpha_per_cm;
          vg_m_k    = soil_properties[soil_num_k].vg_m;
          vg_n_k    = soil_properties[soil_num_k].vg_n;

          double layer_thickness = (cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k-1]);

          //-------------------------------------------
          // do the same for the current state
          double theta = calc_theta_from_h(psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);

          double theta_below = calc_theta_from_h(psi_cm_below, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);

            new_mass += layer_thickness * (theta - theta_below);
            prior_mass += layer_thickness * (theta - theta_below);

          delta_thetas[k] = theta_below; 
          delta_thickness[k] = layer_thickness;
        }

        delta_thetas[layer_num] = next->theta;
        delta_thickness[layer_num] = current->depth_cm - cum_layer_thickness_cm[layer_num-1];

        prior_mass -= *bottom_boundary_flux_above_surface_WFs_cm;

        // theta mass balance computes new theta that conserves the mass; new theta is assigned to the current wetting front
        double theta_new = current->theta;
        if (*bottom_boundary_flux_above_surface_WFs_cm!=0.0){
          theta_new = lgar_theta_mass_balance(layer_num, soil_num, current->psi_cm, new_mass, prior_mass, AET_demand_cm, 
                      delta_thetas, delta_thickness, soil_type, soil_properties);
        }

        current->theta = fmax(soil_properties[soil_type[current->layer_num]].theta_r, fmin(theta_new, soil_properties[soil_type[current->layer_num]].theta_e));

        layer_num = current->layer_num;
        soil_num  = soil_type[layer_num];

        theta_e_k   = soil_properties[soil_num].theta_e;
        theta_r_k   = soil_properties[soil_num].theta_r;
        vg_a_k      = soil_properties[soil_num].vg_alpha_per_cm;
        vg_m_k      = soil_properties[soil_num].vg_m;
        vg_n_k      = soil_properties[soil_num].vg_n;

        if (current->theta == theta_r_k){
          current->theta = current->theta + 1.E-12;
          bottom_boundary_flux_cm -= *bottom_boundary_flux_above_surface_WFs_cm;
          bottom_boundary_flux_cm += (mass_at_start_of_extraction - lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
          *bottom_boundary_flux_above_surface_WFs_cm = (mass_at_start_of_extraction - lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
        }

        double Se = calc_Se_from_theta(current->theta,theta_e_k,theta_r_k);

        if (*bottom_boundary_flux_above_surface_WFs_cm!=0.0){
          current->psi_cm = calc_h_from_Se(Se, vg_a_k, vg_m_k, vg_n_k); 
          *bottom_boundary_flux_above_surface_WFs_cm = *bottom_boundary_flux_above_surface_WFs_cm - (mass_at_start_of_extraction - lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
        }

        free(delta_thetas);
        free(delta_thickness);
      
      }

    }
    else{
      bottom_boundary_flux_cm -= *bottom_boundary_flux_above_surface_WFs_cm; //this addresses the case where, in 1 time step, at the time step start, there was a surface WF which had TO WFs above it, and the surface WF merged with a deeper WF. This means that there would have been a demand that the higher TO WFs should have taken from a surface WF, but now no surface WFs exist. So the percolation is reduced for this time step to conserve mass.
    }

    // // In the event that the rightmost surface WF is not in the first layer, surface WFs that are above the WF that had TO ET extracted from it that have to_bottom == TRUE should also be updated 
    listSendToTop(*head);

    current = listFindFront(listLength(*head) - 1 , *head, NULL);
    if (current == NULL){
      current = *head;
    }
    next = current->next;
    if (next!=NULL){
      for (int wf = listLength(*head) - 1; wf != 1; wf--) {
        if ( (current->to_bottom==TRUE)  ){
          current->psi_cm = next->psi_cm;
        }
        current = listFindFront(wf-1, *head, NULL);
        next = current->next;
      }
    }


    current = *head;
    next = current->next;
    bool deepest_surface_WF_became_TO_in_layer_deeper_than_top = false;
    for (int wf = 1; wf != (listLength(*head)); wf++) {
      if ( (current->is_WF_GW==0) & (next->is_WF_GW==1) & (current->layer_num==next->layer_num) & (current->theta<=next->theta) ){
        current->is_WF_GW = 1;
        if (current->layer_num>1){
          deepest_surface_WF_became_TO_in_layer_deeper_than_top = true;
        }
      }
      current = current->next;
      next = current->next;
    }
    //and then if this happened in a layer below the top one, the to_bottom wetting fronts above the converted one with the same psi value have to also have their is_WF_GW set to 1
    current = listFindFront(listLength(*head),*head,NULL);
    if (deepest_surface_WF_became_TO_in_layer_deeper_than_top){
      for (int wf = listLength(*head); wf != 0; wf--) {
        if (current->next!=NULL){
          if ( current->to_bottom==TRUE && (current->next->is_WF_GW==TRUE) ){
            current->is_WF_GW = 1;
          }
        }
        current = listFindFront(current->front_num-1,*head,NULL);
      }
    }


    current = *head;
    if (current->next!=NULL){
      if ( (current->depth_cm==0.0) && (current->next->depth_cm!=0.0) ){
        listDeleteFront(1, head, soil_type, soil_properties);
      }
    }

    current = listFindFront(wf_from_which_to_extract, *head, NULL);
    int layer_num = current->layer_num;
    int soil_num  = soil_type[layer_num];
    double theta_r = soil_properties[soil_num].theta_r;
    bool all_WFs_to_GW = false;
    if (current->theta<theta_r){
      if (verbosity.compare("high") == 0){
        printf("rare case where extracting recharge from surface WFs would make their theta<theta_r, so setting bottom_boundary_flux_above_surface_WFs_cm to 0 and accordingly increasing bottom_boundary_flux_cm \n");
        printf("current->theta: %lf \n", current->theta);
        printf("theta_r: %lf \n", theta_r);
        printf("bottom_boundary_flux_above_surface_WFs_cm: %lf \n", *bottom_boundary_flux_above_surface_WFs_cm);

      }

      current->is_WF_GW = TRUE;
      current->theta = original_theta;
      current->psi_cm = original_psi_cm;

      bottom_boundary_flux_cm = bottom_bdy_flux_at_start_of_extraction + *bottom_boundary_flux_above_surface_WFs_cm;
      *bottom_boundary_flux_above_surface_WFs_cm = 0.0;
      all_WFs_to_GW = true;
    }
    if (all_WFs_to_GW){//make all WFs GW because the wettest WF could not yield the recharge demand 
      current = *head;
      next = current->next;
      for (int wf = 1; wf != (listLength(*head)); wf++) {
        current->is_WF_GW = TRUE;
        current = current->next;
      }
    }

    if (verbosity.compare("high") == 0) {
      printf("After TO demand subtracted from rightmost surface WF ... \n");
      printf("bottom_boundary_flux_above_surface_WFs_cm: %.17lf \n", *bottom_boundary_flux_above_surface_WFs_cm);
      listPrint(*head); 
      printf("associated mass: %lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
    }

  }

  return(bottom_boundary_flux_cm);

}


extern void lgar_global_theta_update(double bottom_boundary_flux_above_surface_WFs_cm, int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head){
  if (verbosity.compare("high") == 0) {
    printf("before lgar_global_theta_update \n");
    listPrint(*head);
  }
  double vg_a_k, vg_m_k, vg_n_k;
  double theta_e_k, theta_r_k;
  int wf_free_drainage_demand = wetting_front_free_drainage(*head);
  for (int wf = wf_free_drainage_demand; wf > 1; wf--) {
    struct wetting_front *current = listFindFront(wf, *head, NULL);
    struct wetting_front *previous = listFindFront(wf-1, *head, NULL);

    if ((previous!=NULL) && (current!=NULL)){
      if ( (bottom_boundary_flux_above_surface_WFs_cm!=0) & (current->layer_num!=previous->layer_num) ){ 
        previous->psi_cm = current->psi_cm;
        int soil_num_k  = soil_type[previous->layer_num];
        
        theta_e_k = soil_properties[soil_num_k].theta_e;
        theta_r_k = soil_properties[soil_num_k].theta_r;
        vg_a_k    = soil_properties[soil_num_k].vg_alpha_per_cm;
        vg_m_k    = soil_properties[soil_num_k].vg_m;
        vg_n_k    = soil_properties[soil_num_k].vg_n;
        previous->theta = calc_theta_from_h(previous->psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
      }
    }
  }

  for (int wf = wf_free_drainage_demand - 1; wf > 1; wf--) {
  struct wetting_front *current = listFindFront(wf,*head, NULL);
  if ( (current!=NULL) && (current->next!=NULL) ){
    if ( (current->to_bottom==TRUE) && (current->next->is_WF_GW==0) ){ 
      current->psi_cm = current->next->psi_cm;
      int soil_num_k  = soil_type[current->layer_num];
      
      double theta_e_k = soil_properties[soil_num_k].theta_e;
      double theta_r_k = soil_properties[soil_num_k].theta_r;
      double vg_a_k    = soil_properties[soil_num_k].vg_alpha_per_cm;
      double vg_m_k    = soil_properties[soil_num_k].vg_m;
      double vg_n_k    = soil_properties[soil_num_k].vg_n;
      current->theta = calc_theta_from_h(current->psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
      }
    }
  }
  if (verbosity.compare("high") == 0) {
    printf("after lgar_global_theta_update \n");
    listPrint(*head);
  }
}


extern void lgar_global_psi_update(int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head){
  //some of the TO-TO merging fxns update theta via mass balance, so after that, psi must be updated accordingly.
  struct wetting_front *current = *head;
  for (int wf = 1; wf != (listLength(*head)); wf++){
    int soil_num     = soil_type[current->layer_num];
    double theta_e   = soil_properties[soil_num].theta_e;
    double theta_r   = soil_properties[soil_num].theta_r;
    double vg_a      = soil_properties[soil_num].vg_alpha_per_cm;
    double vg_m      = soil_properties[soil_num].vg_m;
    double vg_n      = soil_properties[soil_num].vg_n;

    double Se = calc_Se_from_theta(current->theta,theta_e,theta_r);
    current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);

    current = current->next;

    if (current==NULL){
      break;
    }
  }

  if (verbosity.compare("high") == 0) {
    printf("states after lgar_global_psi_update, which comes after TO merging and layer bdy crossing ... \n");
    listPrint(*head);
    }
}

extern bool lgarto_correct_negative_depths(struct wetting_front** head){
  struct wetting_front *current = *head;
  struct wetting_front *next = current->next;
  bool did_a_WF_have_negative_depth = false;
  for (int wf = 1; wf != (listLength(*head)); wf++){
    if (current->depth_cm<0){
      current->depth_cm = 0.0;
      did_a_WF_have_negative_depth = true;
    }
    current = next;
    next = current->next;
  }
  return(did_a_WF_have_negative_depth);
}

extern void lgarto_cleanup_after_surface_TO_merging_in_layer_below_top(bool merged_in_non_top_layer, int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head){
  if (verbosity.compare("high") == 0) {
    printf("before lgarto_cleanup_after_surface_TO_merging_in_layer_below_top: \n");
    listPrint(*head);
  }
  /////if there was merging of surface and TO WFs where the surface WF was dry enough to become a TO WF, but in a layer deeper than the top one, a TO WF must move to the soil surface and therefore have its theta value corrected based on its psi.
  /////This functionality could eventually be part of some kind of merging or something, and possibly it should be made such that it can work iteratively until it is no longer needed, as opposed to running just once
  if (merged_in_non_top_layer){

    listSendToTop(*head); 

    double vg_a_k, vg_m_k, vg_n_k;
    double theta_e_k, theta_r_k;

    for (int wf = listLength(*head)-1; wf > 0; wf--) {
      struct wetting_front *current = listFindFront(wf, *head, NULL);
      struct wetting_front *next = current->next;

      if (current->to_bottom == TRUE){
        current->is_WF_GW = next->is_WF_GW;
      }

      if ( (next->to_bottom==FALSE) && (current->to_bottom==FALSE) && (current->layer_num!=next->layer_num) ){
        current->layer_num = next->layer_num;
      }

        if ( (current->layer_num>next->layer_num) && (next->to_bottom == TRUE) ){
          current->layer_num = next->layer_num; 
        } 

        int soil_num_k  = soil_type[current->layer_num];
        
        theta_e_k = soil_properties[soil_num_k].theta_e;
        theta_r_k = soil_properties[soil_num_k].theta_r;
        vg_a_k    = soil_properties[soil_num_k].vg_alpha_per_cm;
        vg_m_k    = soil_properties[soil_num_k].vg_m;
        vg_n_k    = soil_properties[soil_num_k].vg_n;
        current->theta = calc_theta_from_h(current->psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);

      if ((current->to_bottom==TRUE) && (next->to_bottom==TRUE) && (current->layer_num!=next->layer_num) && (current->psi_cm!=next->psi_cm)){
        current->psi_cm = next->psi_cm;
        int soil_num_k  = soil_type[current->layer_num];
        
        theta_e_k = soil_properties[soil_num_k].theta_e;
        theta_r_k = soil_properties[soil_num_k].theta_r;
        vg_a_k    = soil_properties[soil_num_k].vg_alpha_per_cm;
        vg_m_k    = soil_properties[soil_num_k].vg_m;
        vg_n_k    = soil_properties[soil_num_k].vg_n;
        current->theta = calc_theta_from_h(current->psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
      }//to make adjacent to_bottom WFs have the correct psi and theta values

    }
  }
  if (verbosity.compare("high") == 0) {
    printf("after lgarto_cleanup_after_surface_TO_merging_in_layer_below_top: \n");
    listPrint(*head);
  }
}



extern int lgarto_correction_type(int num_layers, double* cum_layer_thickness_cm, struct wetting_front** head){
  int correction_type = 0;
  struct wetting_front *current = *head;
  struct wetting_front *next = current->next;
  struct wetting_front *next_to_next = NULL;
  bool TO_layer_cross = FALSE;
  struct wetting_front *top_most_TO_front_below_surfs = NULL;
  double top_most_TO_front_below_surfs_psi_cm = 1.E16;
  if (listLength(*head)>(listLength_surface(*head)+listLength_TO_WFs_above_surface_WFs(*head))){
    top_most_TO_front_below_surfs = listFindFront(listLength_surface(*head)+listLength_TO_WFs_above_surface_WFs(*head) + 1, *head, NULL);
    top_most_TO_front_below_surfs_psi_cm = top_most_TO_front_below_surfs->psi_cm;
  }
  if (next!=NULL){
    next_to_next = current->next->next;
  }

  for (int wf = 1; wf != (listLength(*head)); wf++) {
    struct wetting_front *previous = listFindFront(wf - 1, *head, NULL);

    if (next!=NULL){
      if ( ((current->depth_cm>next->depth_cm) && (current->is_WF_GW==0) && (next->is_WF_GW==1) && (next->to_bottom==0)) ){
        correction_type = 1; //this is surface-TO merging 
      }
      if ( (current->depth_cm>next->depth_cm) && (current->is_WF_GW==1) && (next->is_WF_GW==1) && (next->to_bottom==FALSE) && (current->to_bottom==FALSE) && (current->theta < next->theta) ){ 
        correction_type = 2; //TO WFs merge when one gets too deep
      }
      if ( ((current->theta>next->theta) ) && (current->is_WF_GW==1) && (next->is_WF_GW==1) && (next->layer_num==current->layer_num) ){ 
        correction_type = 3; //TO WFs merge when one gets too dry, added (current->theta<next->theta && fabs(current->psi_cm - next->psi_cm)<1.E-2) to make TO WFs merge when they are very close in theta (avoids problems with enormous changes in depth when thetas are very close)
      }
    }
    if (previous!=NULL){
      if ( current->depth_cm < cum_layer_thickness_cm[current->layer_num - 1] && (previous->to_bottom==TRUE) && (current->layer_num>previous->layer_num) && (current->is_WF_GW==TRUE) ){
        correction_type = 4; //this is a TO WF crossing a layer bdy (up)
        TO_layer_cross = TRUE;
      }
    }

    if ( (current->depth_cm > cum_layer_thickness_cm[current->layer_num]) && (next->to_bottom==TRUE) && (current->is_WF_GW==TRUE) && (current->layer_num!=num_layers) ){
      correction_type = 4; //this is a TO WF crossing a layer bdy (down)
      TO_layer_cross = TRUE;
    }

    if (next!=NULL){
      if ( (current->is_WF_GW==0) && (next->is_WF_GW==0) && (current->theta>next->theta) && (current->depth_cm > next->depth_cm) && (current->layer_num == next->layer_num) && (!next->to_bottom) ){
        correction_type = 5; //this is surface-surface WF merging 
      }
      if ( (current->is_WF_GW==0) && (current->depth_cm > cum_layer_thickness_cm[current->layer_num]) && (next->depth_cm == cum_layer_thickness_cm[current->layer_num]) && (current->theta>next->theta) && (current->layer_num!=num_layers) && (current->psi_cm<top_most_TO_front_below_surfs_psi_cm) ){
        correction_type = 6; //this is surface WF layer bdy crossing 
      }
    }
    if ( (next_to_next == NULL) && (current->depth_cm > cum_layer_thickness_cm[current->layer_num]) && (current->is_WF_GW==0) ){
      correction_type = 7; //this is a surface WF crossing the model lower bdy
    }

    if (current->is_WF_GW==FALSE && next->is_WF_GW==TRUE && current->theta<next->theta && current->layer_num==next->layer_num && correction_type==0){
      correction_type = 8;//fairly confident that correction type 8 is no longer necessary due to new functionality of listDeleteFront, will test removal when time is available
    }

    if (current->to_bottom==TRUE && current->is_WF_GW==FALSE && current->next->is_WF_GW==TRUE){
      correction_type = 8;
    }

    current = next;
    next = current->next;
    if (next!=NULL){
      next_to_next = current->next->next;
    }
  }

  if (TO_layer_cross){
    correction_type = 4;//idea is that in some rare cases, a crash will happen if correction types 2 or 3 are necessary at the same time as type 4, and types 2 or 3 are attempted first
  }
  if (verbosity.compare("high") == 0){
    printf("computed correction type: %d \n", correction_type);
  }
  return correction_type;
}

extern int lgarto_correction_type_surf(int num_layers, double* cum_layer_thickness_cm, struct wetting_front** head){
  int correction_type_surf = 0;
  struct wetting_front *current = *head;
  struct wetting_front *next = current->next;
  struct wetting_front *next_to_next = NULL;
  struct wetting_front *top_most_TO_front_below_surfs = NULL;
  double top_most_TO_front_below_surfs_psi_cm = 1.E16;
  if (listLength(*head)>(listLength_surface(*head)+listLength_TO_WFs_above_surface_WFs(*head))){
    top_most_TO_front_below_surfs = listFindFront(listLength_surface(*head)+listLength_TO_WFs_above_surface_WFs(*head) + 1, *head, NULL);
    top_most_TO_front_below_surfs_psi_cm = top_most_TO_front_below_surfs->psi_cm;
  }

  if (next!=NULL){
    next_to_next = current->next->next;
  }

  for (int wf = 1; wf != (listLength(*head)); wf++) {

    if (next!=NULL){
      if ( (current->is_WF_GW==0) && (next->is_WF_GW==0) && (current->theta>next->theta) && (current->depth_cm > next->depth_cm) && (current->layer_num == next->layer_num) && (!next->to_bottom) ){
        correction_type_surf = 1; //this is surface-surface WF merging 
        break;
      }
      if ( (current->is_WF_GW==0) && (current->depth_cm > cum_layer_thickness_cm[current->layer_num]) && (next->depth_cm == cum_layer_thickness_cm[current->layer_num]) && (current->theta>next->theta) && (current->layer_num!=num_layers) && (current->psi_cm<top_most_TO_front_below_surfs_psi_cm) ){
        correction_type_surf = 2; //this is surface WF layer bdy crossing 
        break;
      }
    }
    if ( (next_to_next == NULL) && (current->depth_cm > cum_layer_thickness_cm[current->layer_num]) && (current->is_WF_GW==0) ){
      correction_type_surf = 3; //this is a surface WF crossing the model lower bdy
      break;
    }
    if (lgar_check_dry_over_wet_wetting_fronts(*head)){
      correction_type_surf = 4;
      break;
    }

    current = next;
    next = current->next;
    if (next!=NULL){
      next_to_next = current->next->next;
    }
  }

  if (verbosity.compare("high") == 0){
    printf("computed correction type for surface WFs: %d \n", correction_type_surf);
  }
  return correction_type_surf;
}


extern double lgarto_TO_WFs_merge_via_depth(double initial_mass, double column_depth, double* cum_layer_thickness_cm, struct wetting_front** head, int *soil_type, struct soil_properties_ *soil_properties){
  //in the event that a TO WF passes another one. 

  struct wetting_front *current = *head;
  struct wetting_front *next = current->next;

  for (int wf = 1; wf != (listLength(*head)); wf++){

    if ( (current->depth_cm>next->depth_cm) && (current->is_WF_GW==1) && (next->is_WF_GW==1) && (next->to_bottom==FALSE) && (current->to_bottom==FALSE) && (current->theta < next->theta) ){
      next = listDeleteFront(next->front_num, head, soil_type, soil_properties);
      next = current->next;
      //possible current should be made to next (but in this case next be made to current->next)
        
        double temp_tol = 1e-10;
        double factor;

        if (current->depth_cm>column_depth){
          factor = 1000.0;
        }
        else{
          factor = 1.0; //seems to control how many iterations are needed
        }
        bool switched = false;

        int iter = 0;

        while ( fabs(lgar_calc_mass_bal(cum_layer_thickness_cm, *head)-initial_mass)>temp_tol ){
          iter++;
          if (iter>10000){/////should have more sophisticated loop breaking as in LGAR while loops
            break;
          }
          if (lgar_calc_mass_bal(cum_layer_thickness_cm, *head)>=initial_mass) {
            current->depth_cm += 0.1 * factor;
            switched = false;
          }
          else {
            if (!switched) {
              switched = true;
              factor = factor * 0.1;
            }
            current->depth_cm -= 0.1 * factor;
          }

          if ((next->depth_cm>column_depth)){
            break;
          }
        }
    }

    current = current->next;
    next = current->next;
    if (next==NULL){
      break;
    }
  }

  double mass_diff = initial_mass - lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  //might need to factor in some fluxes too, in case you're seeing a mass balance here
  // if (current->depth_cm > column_depth){
  //   mass_diff = 0.0;
  // }

  if (verbosity.compare("high") == 0) {
    printf("after lgarto_TO_WFs_merge_via_depth ... \n");
    listPrint(*head);
  }

  return mass_diff;
}


extern double lgarto_TO_WFs_merge_via_theta(double initial_mass, double column_depth, double* cum_layer_thickness_cm, struct wetting_front** head, int *soil_type, struct soil_properties_ *soil_properties){
  if (verbosity.compare("high") == 0) {
      printf("before lgarto_TO_WFs_merge_via_theta ... \n");
      listPrint(*head);
    }

  //in the event that a TO WF becomes drier than the one to its left.
  struct wetting_front *current = *head;
  struct wetting_front *next = current->next;
  struct wetting_front *previous = *head;
  double mass_diff = 0.0;

  for (int wf = 1; wf != (listLength(*head)); wf++){

    if (current->next==NULL){
      break;
    }

    if (  ((current->theta>next->theta)) && (current->is_WF_GW==1) && (next->is_WF_GW==1) && (next->layer_num==current->layer_num)  ){ 

      if (current->depth_cm == 0.0){
        current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
        //possible current should be made to next , although this one has been handled by this point 
        double new_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        mass_diff = initial_mass - new_mass; 
      }
      else {
        if (next->to_bottom==FALSE && previous->to_bottom==FALSE){
          current = listDeleteFront(current->front_num, head, soil_type, soil_properties);

          double temp_tol = 1e-12;
          double factor = 1.0; 
          bool switched = false;
          while ( fabs(lgar_calc_mass_bal(cum_layer_thickness_cm, *head)-initial_mass)>temp_tol ){
            if (lgar_calc_mass_bal(cum_layer_thickness_cm, *head)>=initial_mass) {
              current->depth_cm += 0.1 * factor;
              switched = false;
            }
            else {
              if (!switched) {
                switched = true;
                factor = factor * 0.1;
              }
              current->depth_cm -= 0.1 * factor;
            }

            if (current->front_num>1){
              if (listFindFront(current->front_num - 1, *head, NULL)->depth_cm>current->depth_cm){ //the idea here is that in this case, should let the WFs merge based on depth and not theta
                break;
              }
            }

            if ((next->depth_cm>column_depth)){
              break;
            }

            if (current->depth_cm<0.0){
              current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
              double new_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
              mass_diff = initial_mass - new_mass; 
              break;
            }

          }
        }
        else{
          current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
          double new_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
          mass_diff = initial_mass - new_mass; 
          //the idea here is that when the conditions for TO-TO merging based on theta are correct by next->to_bottom==FALSE, you can't use the above code because you can't delete a to_bottom == TRUE WF. 
          //you could also solve this problem by iteratively updating the moisture of the to_bottom WF below and all WFs below it that had the same psi value. 
        }
      }
    }
    previous = current;
    if (previous->front_num==listLength(*head)){
      break;
    }
    current = current->next;
    next = current->next;
    }

    if (verbosity.compare("high") == 0) {
      printf("after lgarto_TO_WFs_merge_via_theta ... \n");
      listPrint(*head);
    }
  return mass_diff;
}


// ############################################################################################
/* The function allows TO wetting fronts to cross layer boundaries. */
// ############################################################################################

extern void lgar_TO_wetting_fronts_cross_layer_boundary(int *front_num_with_negative_depth, int num_layers, double *cum_layer_thickness_cm, int *soil_type, double *frozen_factor, struct soil_properties_ *soil_properties, struct wetting_front** head){
      //here TO WFs cross layer boundaries if they got to a sufficient depth
      
    for (int wf = (listLength(*head)); wf != (listLength_surface(*head)); wf--) { 
      struct wetting_front *current = listFindFront(wf, *head, NULL);
      struct wetting_front *previous = listFindFront(wf-1, *head, NULL);
      struct wetting_front *next = listFindFront(wf+1, *head, NULL);
      int layer_num = current->layer_num;
      // // first the case where a TO WF moves from a higher layer to the one below
      if (current->depth_cm > cum_layer_thickness_cm[layer_num] && (next->to_bottom==TRUE) && (current->is_WF_GW==TRUE) && (layer_num!=num_layers) ){
        if (verbosity.compare("high") == 0) {
          printf("Inside layer boundary crossing for TO wetting front, moving down case ... \n");
          listPrint(*head);
        }

        // local variables
        double theta_e,theta_r;
        double vg_a, vg_m, vg_n;
        int layer_num, soil_num;//

        layer_num   = current->layer_num;
        soil_num    = soil_type[layer_num];
        theta_e     = soil_properties[soil_num].theta_e;
        theta_r     = soil_properties[soil_num].theta_r;
        vg_a        = soil_properties[soil_num].vg_alpha_per_cm;
        vg_m        = soil_properties[soil_num].vg_m;
        vg_n        = soil_properties[soil_num].vg_n;
        double Ksat_cm_per_h  = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num]; 

        struct wetting_front *next_to_next = current->next->next; 
        if (layer_num==num_layers){
          next_to_next=next;//in LGARTO for surface WFs, it was next_to_next==next, which is not correct 
        }


        double current_theta = fmin(theta_e, current->theta);
        double overshot_depth = current->depth_cm - next->depth_cm;
        int soil_num_next;
        if (layer_num==num_layers){
           soil_num_next = soil_type[layer_num];
        }
        else{
          soil_num_next = soil_type[layer_num+1];
        }

        double next_theta_e   = soil_properties[soil_num_next].theta_e;
        double next_theta_r   = soil_properties[soil_num_next].theta_r;
        double next_vg_a      = soil_properties[soil_num_next].vg_alpha_per_cm;
        double next_vg_m      = soil_properties[soil_num_next].vg_m;
        double next_vg_n      = soil_properties[soil_num_next].vg_n;

        double Se = calc_Se_from_theta(current->theta,theta_e,theta_r);
        current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);

        current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m); // AJ check Ksat_cm_per_h here current or next? // PTL: mixed and matched: one layer vg_m and another layer K_s value, fixed

        // current psi with van Genuchten properties of the next layer to get new theta
        double theta_new = calc_theta_from_h(current->psi_cm, next_vg_a, next_vg_m, next_vg_n, next_theta_e, next_theta_r);

        double mbal_correction = overshot_depth * (current_theta - next->theta);
        double mbal_Z_correction = mbal_correction / (theta_new - next_to_next->theta); // this is the new wetting front depth

        double depth_new = cum_layer_thickness_cm[layer_num] + mbal_Z_correction; // this is the new wetting front absolute depth

        if (isinf(depth_new)){
          theta_new = theta_new - 1.E-14; 
          mbal_correction = overshot_depth * (current_theta - next->theta);
          mbal_Z_correction = mbal_correction / (theta_new - next_to_next->theta); // this is the new wetting front depth
          depth_new = cum_layer_thickness_cm[layer_num] + mbal_Z_correction; // this is the new wetting front absolute depth
        }

        current->depth_cm = cum_layer_thickness_cm[layer_num];

        next->theta = theta_new;
        next->psi_cm = current->psi_cm;
        next->depth_cm = depth_new;
        
        next->layer_num = layer_num + 1;
        next->dzdt_cm_per_h = current->dzdt_cm_per_h;
        current->dzdt_cm_per_h = 0;
        current->to_bottom = 1;
        next->to_bottom = 0;

        if (isnan(next->depth_cm)){
          next = listDeleteFront(next->front_num, head, soil_type, soil_properties);
          next = current->next;
        }

        if (verbosity.compare("high") == 0) {
          printf("After layer boundary crossing for TO wetting front, moving down case ... \n");
          listPrint(*head);
        }

      }

      // // next the case where a TO WF moves from a lower layer to the one above
      if (previous!=NULL){
        if (current->depth_cm < cum_layer_thickness_cm[layer_num - 1] && (previous->to_bottom==TRUE) && (current->layer_num>previous->layer_num) && (current->is_WF_GW==TRUE) ){
          if (verbosity.compare("high") == 0) {
            printf("this is the TO WF moving to a higher layer code \n");
            printf("front_num: %d \n",current->front_num);
            listPrint(*head);
          }

          double theta_e,theta_r;
          double vg_a, vg_m, vg_n;
          int layer_num, soil_num;

          layer_num   = current->layer_num;
          soil_num    = soil_type[layer_num];
          theta_e     = soil_properties[soil_num].theta_e;
          theta_r     = soil_properties[soil_num].theta_r;
          vg_a        = soil_properties[soil_num].vg_alpha_per_cm;
          vg_m        = soil_properties[soil_num].vg_m;
          vg_n        = soil_properties[soil_num].vg_n;
          double Ksat_cm_per_h  = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num]; 

          double current_theta = fmin(theta_e, current->theta);
          double overshot_depth = previous->depth_cm - current->depth_cm;
          int soil_num_above = soil_type[layer_num - 1];

          double above_theta_e   = soil_properties[soil_num_above].theta_e;
          double above_theta_r   = soil_properties[soil_num_above].theta_r;
          double above_vg_a      = soil_properties[soil_num_above].vg_alpha_per_cm;
          double above_vg_m      = soil_properties[soil_num_above].vg_m;
          double above_vg_n      = soil_properties[soil_num_above].vg_n;

          double Se = calc_Se_from_theta(current->theta,theta_e,theta_r);
          current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);

          current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m); // AJ check Ksat_cm_per_h here current or next? // PTL: mixed and matched: one layer vg_m and another layer K_s value

          // current psi with van Genuchten properties of the next layer to get new theta
          double theta_new = calc_theta_from_h(current->psi_cm, above_vg_a, above_vg_m, above_vg_n, above_theta_e, above_theta_r);
          double new_boundary_theta = calc_theta_from_h(next->psi_cm, above_vg_a, above_vg_m, above_vg_n, above_theta_e, above_theta_r);

          double mbal_correction = overshot_depth * (next->theta - current_theta);
          double mbal_Z_correction = mbal_correction / (theta_new - new_boundary_theta); 

          double depth_new = cum_layer_thickness_cm[layer_num - 1] + mbal_Z_correction; // this is the new wetting front absolute depth

          current->depth_cm = cum_layer_thickness_cm[layer_num - 1];

          current->theta = new_boundary_theta;

          current->psi_cm = next->psi_cm;
          previous->depth_cm = depth_new;

          current->layer_num = previous->layer_num;
          previous->dzdt_cm_per_h = current->dzdt_cm_per_h;
          current->dzdt_cm_per_h = 0;
          current->to_bottom = 1;
          previous->to_bottom = 0;

          if (isnan(previous->depth_cm)){
            listPrint(*head);
            printf("previous->depth_cm was not a number \n");
            abort();
          }

          if (previous->depth_cm<0.0){
            *front_num_with_negative_depth = previous->front_num;
            if (verbosity.compare("high") == 0) {
              printf("negative WF depth. front_num_with_negative_depth: %d \n", *front_num_with_negative_depth);
              listPrint(*head);
            }
          }
          else {
            *front_num_with_negative_depth = -1;
          }

          if (verbosity.compare("high") == 0) {
            printf("this is after the TO WF moving to a higher layer code \n");
            printf("front_num: %d \n",current->front_num);
            listPrint(*head);
          }

        }
      }
    }
}

extern void lgarto_consolidate_excessive_fronts(double* cum_layer_thickness_cm, struct wetting_front** head, int *soil_type, struct soil_properties_ *soil_properties){
  double start_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  if (verbosity.compare("high") == 0) {
    printf("this is how many WFs are in the zone for consolidation: %d \n", lgarto_count_fronts_for_excessive_calc(head));
    printf("states before lgarto_consolidate_excessive_fronts: \n");
    listPrint(*head);
  }
  while (lgarto_count_fronts_for_excessive_calc(head) > 6){
    listDeleteFront(1, head, soil_type, soil_properties);
  }
  double end_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
  struct wetting_front *current = *head;
  while (start_mass>end_mass){
    current->depth_cm -= 0.001;//one directional mass balance correction loop should probably be iterative, but less important because it only goes in 1 direction
  }

  if (verbosity.compare("high") == 0) {
    printf("states after lgarto_consolidate_excessive_fronts: \n");
    listPrint(*head);
  }
}


extern double adjust_new_theta(int new_wf_num, double target_mass, double *cum_layer_thickness, int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head)
{

  double bottom_boundary_flux_cm_temp = 0.0;
  struct wetting_front *current;
  current = *head;
  while (current->front_num!=new_wf_num){
    current = current->next;
  }
  current = listFindFront(listLength_TO_WFs_above_surface_WFs(*head)+1, *head, NULL);

  double current_mass = lgar_calc_mass_bal(cum_layer_thickness, *head);
  double delta_mass = fabs(current_mass - target_mass);
  double tolerance = 1e-12;
  bool switched = false;
  double factor = 0.1;

  while (delta_mass>tolerance){
    if (current_mass<target_mass) {
      current->theta += 0.1 * factor;
      switched = false;
    }
    else {
      if (!switched) {
        switched = true;
        factor = factor * 0.1;
      }
      current->theta -= 0.1 * factor;
    }
    current_mass = lgar_calc_mass_bal(cum_layer_thickness, *head);
    delta_mass = fabs(current_mass - target_mass);
  }

  int soil_num = soil_type[current->layer_num];

  double theta_e     = soil_properties[soil_num].theta_e;
  double theta_r     = soil_properties[soil_num].theta_r;
  double vg_a        = soil_properties[soil_num].vg_alpha_per_cm;
  double vg_m        = soil_properties[soil_num].vg_m;
  double vg_n        = soil_properties[soil_num].vg_n;

  double theta_e_temp = soil_properties[soil_num].theta_e; 

  double Se = 0.0;
  if (current->theta < theta_e_temp){
    Se              = calc_Se_from_theta(current->theta, theta_e, theta_r);
    current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n); //have to update psi because you just updated theta 
  }

  if (current->theta >= theta_e_temp){
    current->theta = theta_e_temp;
    current->psi_cm = 0.0;
    struct wetting_front *next;
    next = current->next;
    if (next->theta==theta_e_temp){
      int end_loop = current->front_num;
      for (int temp_wf = 1; temp_wf <= end_loop; temp_wf++){
        current = listFindFront(temp_wf, *head, NULL);
        current->depth_cm = 0;
        current->is_WF_GW = 1;
      }
    }

  }


  current = *head;
  if (listLength_surface(*head)>0){
    while(current->is_WF_GW==TRUE){
      current = current->next;
    }
  }
  struct wetting_front *next;
  next = current->next;

  if ((current->theta<next->theta) && (listLength_surface(*head)>0)){

    
    soil_num = soil_type[current->layer_num];
    theta_e_temp = soil_properties[soil_num].theta_e;
    if (current->theta<next->theta){
      current->theta = next->theta + 0.01;
      if (current->theta >= theta_e_temp){
        current->theta = theta_e_temp;
        current->psi_cm = 0.0;
      }
    }
    if (current->theta != theta_e_temp){
      int layer_num   = current->layer_num;
      soil_num    = soil_type[layer_num];
      theta_e     = soil_properties[soil_num].theta_e;
      theta_r     = soil_properties[soil_num].theta_r;
      vg_a        = soil_properties[soil_num].vg_alpha_per_cm;
      vg_m        = soil_properties[soil_num].vg_m;
      vg_n        = soil_properties[soil_num].vg_n;
      double theta       = current->theta;
      double new_Se          = calc_Se_from_theta(theta, theta_e, theta_r);
      current->psi_cm      = calc_h_from_Se(new_Se, vg_a, vg_m, vg_n);
    }

      double temp_mass = lgar_calc_mass_bal(cum_layer_thickness, *head);
      while (temp_mass>(target_mass + 0*bottom_boundary_flux_cm_temp)){
        current->depth_cm = current->depth_cm - 0.001; //again, one directional mass bal loop is not ideal for now but at least we know it's going to mnly move in 1 direction
        if (current->depth_cm<0.0){
          break;
        }
      }
      if (current->depth_cm<0.0){
        current->depth_cm = 1e-4;
      }
  }
  bottom_boundary_flux_cm_temp = target_mass - lgar_calc_mass_bal(cum_layer_thickness, *head);

  if (isnan(bottom_boundary_flux_cm_temp)){
      printf("bottom_boundary_flux_cm_temp is not a number \n");
      listPrint(*head);
      abort();
  }

  return bottom_boundary_flux_cm_temp;
}

extern void lgar_clean_redundant_fronts(struct wetting_front** head, int *soil_type, struct soil_properties_ *soil_properties){
  if (verbosity.compare("high") == 0) {
    printf("before lgar_clean_redundant_fronts: \n");
    listPrint(*head);
  }
  struct wetting_front *current;
  struct wetting_front *next;
  current = *head;
  next = current->next;
  for (int wf = 1; wf != (listLength(*head)); wf++) {
    if ( ((current->layer_num==next->layer_num) && (fabs(current->theta - next->theta)<1.E-10)) ){
      current = listDeleteFront(current->front_num, head, soil_type, soil_properties); 
      break;
    }

    current = next;
    next = current->next;
  }

  if (verbosity.compare("high") == 0) {
    printf("after lgar_clean_redundant_fronts: \n");
    listPrint(*head);
  }
}

extern bool correct_close_psis(int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head){
  if (verbosity.compare("high") == 0) {
    printf("before lgar_clean_redundant_fronts: \n");
    listPrint(*head);
  }
  bool close_psis = false;
  struct wetting_front *current;
  struct wetting_front *next;
  current = *head;
  next = current->next;
  for (int wf = 1; wf != (listLength(*head)); wf++) {
    if ( (current->layer_num==next->layer_num) && (current->is_WF_GW==FALSE && next->is_WF_GW==TRUE) && ((fabs(current->psi_cm - next->psi_cm)<1.E-3)  )){
      current = listDeleteFront(current->front_num, head, soil_type, soil_properties); 
      close_psis = true;
      break;
    }
    current = next;
    next = current->next;
  }
  if (verbosity.compare("high") == 0) {
    printf("after correct_close_psis: \n");
    listPrint(*head);
  }
  return(close_psis);
}

extern double lgarto_calc_largest_mag_TO_dzdt(struct wetting_front* head){
  //edited to return mag not just of TO WFs but of all WFs, idea is that if a surf WF is too fast it can still be bad for recharge 
  struct wetting_front *current;
  current = head;

  double largest_mag_TO_dzdt = 0.0;

  while (current!=NULL){
    double temp_dzdt = 0.0;
    temp_dzdt = fabs(current->dzdt_cm_per_h);
    if (temp_dzdt>largest_mag_TO_dzdt){
      largest_mag_TO_dzdt = temp_dzdt;
    }
    current = current->next;
  }

  return(largest_mag_TO_dzdt);

}


#endif