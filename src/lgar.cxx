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

  InitFromConfigFile(config_file, state);
  state->lgar_bmi_params.shape[0] = state->lgar_bmi_params.num_layers;
  state->lgar_bmi_params.shape[1] = state->lgar_bmi_params.num_wetting_fronts;

  // initial number of wetting fronts are same are number of layers
  state->lgar_bmi_params.num_wetting_fronts = state->lgar_bmi_params.num_layers;
  state->lgar_bmi_params.soil_depth_wetting_fronts = new double[state->lgar_bmi_params.num_wetting_fronts];
  state->lgar_bmi_params.soil_moisture_wetting_fronts = new double[state->lgar_bmi_params.num_wetting_fronts];

  // initialize thickness/depth and soil moisture of wetting fronts (used for model coupling)
  struct wetting_front *current = head;
  for (int i=0; i<state->lgar_bmi_params.num_wetting_fronts; i++) {
    assert (current != NULL);
    state->lgar_bmi_params.soil_moisture_wetting_fronts[i] = current->theta;
    state->lgar_bmi_params.soil_depth_wetting_fronts[i] = current->depth_cm * state->units.cm_to_m;
    current = current->next;
  }


  /* initialize bmi input variables to -1.0 (on purpose), this should be assigned (non-negative) and if not,
     the code will throw an error in the Update method */
  state->lgar_bmi_input_params->precipitation_mm_per_h = -1.0;
  state->lgar_bmi_input_params->PET_mm_per_h = -1.0;

  // initialize all global mass balance variables to zero
  state->lgar_mass_balance.volprecip_cm = 0.0;
  state->lgar_mass_balance.volin_cm = 0.0;
  state->lgar_mass_balance.volend_cm = 0.0;
  state->lgar_mass_balance.volAET_cm = 0.0;
  state->lgar_mass_balance.volrech_cm = 0.0;
  state->lgar_mass_balance.volrunoff_cm = 0.0;
  state->lgar_mass_balance.volrunoff_giuh_cm = 0.0;
  state->lgar_mass_balance.volQ_cm = 0.0;
  state->lgar_mass_balance.volPET_cm = 0.0;
  state->lgar_mass_balance.volon_cm = 0.0;

  // setting volon and precip at the initial time to 0.0 as they determine the creation of surficail wetting front
  state->lgar_mass_balance.volon_timestep_cm = 0.0;
  state->lgar_bmi_params.precip_previous_timestep_cm = 0.0;

  // setting flux from groundwater_reservoir_to_stream to zero, will be non-zero when groundwater reservoir is added/simulated
  state->lgar_mass_balance.volQ_gw_timestep_cm = 0.0;
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
  @param ponded_depth_cm        : amount of water on the surface not available for surface drainage (initialized to zero)
  @param ponded_depth_max cm    : maximum amount of water on the surface not available for surface drainage (default is zero)
  @param nint                   : number of trapezoids used in integrating the Geff function (set to 120)
  @param time_s                 : current time [s] (initially set to zero)
  @param sft_coupled            : model coupling flag. if true, lasam is coupled to soil freeze thaw model; default is uncoupled version
  @param giuh_ordinates         : geomorphological instantaneous unit hydrograph
  @param num_giuh_ordinates     : number of giuh ordinates
*/

// #############################################################################################################################
extern void InitFromConfigFile(string config_file, struct model_state *state)
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
	std::cerr<<"Verbosity is set to \' "<<verbosity<<"\' \n";
	std::cerr<<"          *****         \n";
      }

      fp.clear();
      break;
    }
  }


  if (verbosity.compare("none") != 0) {
    std::cerr<<"------------- Initialization from config file ---------------------- \n";
  }

  // setting these options to false (defualt) 
  state->lgar_bmi_params.sft_coupled       = false;
  state->lgar_bmi_params.use_closed_form_G = false;
  
  bool is_layer_thickness_set      = false;
  bool is_initial_psi_set          = false;
  bool is_timestep_set             = false;
  bool is_endtime_set              = false;
  bool is_forcing_resolution_set   = false;
  bool is_layer_soil_type_set      = false;
  bool is_wilting_point_psi_cm_set = false;
  bool is_soil_params_file_set     = false;
  bool is_max_soil_types_set       = false;
  bool is_giuh_ordinates_set       = false;
  bool is_soil_z_set               = false;
  bool is_ponded_depth_max_cm_set  = false;

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

      // calculate the cumulative (absolute) depth from land surface to bottom of each soil layer
      state->lgar_bmi_params.cum_layer_thickness_cm[0] = 0;

      for (unsigned int layer=1; layer <= vec.size(); layer++) {
      	state->lgar_bmi_params.layer_soil_type[layer] = vec[layer-1];
      }

      is_layer_soil_type_set = true;

      continue;
    }
    else if (param_key == "giuh_ordinates") {
      vector<double> vec = ReadVectorData(param_value);

      state->lgar_bmi_params.giuh_ordinates = new double[vec.size()+1];

      for (unsigned int layer=1; layer <= vec.size(); layer++) {
      	state->lgar_bmi_params.giuh_ordinates[layer] = vec[layer-1];
      }

      state->lgar_bmi_params.num_giuh_ordinates = vec.size();

      is_giuh_ordinates_set = true;

      if (verbosity.compare("high") == 0) {
	for (int i=1; i<=state->lgar_bmi_params.num_giuh_ordinates; i++)
	  std::cerr<<"GIUH ordinates : "<<state->lgar_bmi_params.giuh_ordinates[i]<<"\n";

	std::cerr<<"          *****         \n";
      }

      continue;
    }
    else if (param_key == "soil_z") {
      vector<double> vec = ReadVectorData(param_value);

      state->lgar_bmi_params.soil_temperature_z = new double[vec.size()];

      for (unsigned int i=0; i < vec.size(); i++) {
      	state->lgar_bmi_params.soil_temperature_z[i] = vec[i];
      }

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
    else if (param_key == "max_soil_types") {
      state->lgar_bmi_params.num_soil_types = stod(param_value);
      is_max_soil_types_set = true;
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
    else if (param_key == "use_closed_form_G") { 
      if (param_value == "false") {
        state->lgar_bmi_params.use_closed_form_G = false;
      }
      else if (param_value == "true") {
        state->lgar_bmi_params.use_closed_form_G = true;
      }
      else {
	std::cerr<<"Invalid option: use_closed_form_G must be true or false. \n";
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

      assert (state->lgar_bmi_params.timestep_h > 0);
      is_timestep_set = true;

      if (verbosity.compare("high") == 0) {
	std::cerr<<"Model timestep [hours,seconds]: "<<state->lgar_bmi_params.timestep_h<<" , "<<state->lgar_bmi_params.timestep_h*3600<<"\n";
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
	state->lgar_bmi_params.forcing_resolution_h /= 3600; // convert to hours
      else if (param_unit == "[min]" || param_unit == "[minute]")
	state->lgar_bmi_params.forcing_resolution_h /= 60; // convert to hours
      else if (param_unit == "[h]" || param_unit == "[hr]")
	state->lgar_bmi_params.forcing_resolution_h /= 1.0; // convert to hours

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
  
  if(!is_max_soil_types_set)
     state->lgar_bmi_params.num_soil_types = 15;          // maximum number of soil types defaults to 15

  if (verbosity.compare("high") == 0) {
    std::cerr<<"Maximum number of soil types: "<<state->lgar_bmi_params.num_soil_types<<"\n";
    std::cerr<<"          *****         \n";
  }

  if(is_soil_params_file_set) {
    //allocate memory to create an array of structures to hold the soils properties data.
    //state->soil_properties = (struct soil_properties_*) malloc((state->lgar_bmi_params.num_layers+1)*sizeof(struct soil_properties_));


    state->soil_properties = new soil_properties_[state->lgar_bmi_params.num_soil_types+1];
    int num_soil_types = state->lgar_bmi_params.num_soil_types;
    double wilting_point_psi_cm = state->lgar_bmi_params.wilting_point_psi_cm;
    int max_num_soil_in_file = lgar_read_vG_param_file(soil_params_file.c_str(), num_soil_types, wilting_point_psi_cm, state->soil_properties);

    // check if soil layers provided are within the range
    for (int layer=1; layer <= state->lgar_bmi_params.num_layers; layer++) {
      assert (state->lgar_bmi_params.layer_soil_type[layer] <= state->lgar_bmi_params.num_soil_types);
      assert (state->lgar_bmi_params.layer_soil_type[layer] <= max_num_soil_in_file);
    }

    if (verbosity.compare("high") == 0) {
      for (int layer=1; layer<=state->lgar_bmi_params.num_layers; layer++) {
	int soil = state->lgar_bmi_params.layer_soil_type[layer];// layer_soil_type[layer];
	std::cerr<<"Soil type/name : "<<state->lgar_bmi_params.layer_soil_type[layer]<<" "<<state->soil_properties[soil].soil_name<<"\n";
      }
      std::cerr<<"          *****         \n";
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

  if (!is_endtime_set) {
    stringstream errMsg;
    errMsg << "end time not set in the config file "<< config_file << "\n";
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


  if (!is_giuh_ordinates_set) {
    stringstream errMsg;
    errMsg << "giuh ordinates not set in the config file "<< config_file << "\n";
    throw runtime_error(errMsg.str());
  }

  if (state->lgar_bmi_params.sft_coupled) {
    state->lgar_bmi_params.soil_temperature = new double[state->lgar_bmi_params.num_cells_temp]();
    if (!is_soil_z_set) {
      stringstream errMsg;
      errMsg << "soil_z not set in the config file "<< config_file << "\n";
      throw runtime_error(errMsg.str());
    }
  }
  else {
    state->lgar_bmi_params.soil_temperature = new double[1]();
    state->lgar_bmi_params.soil_temperature_z = new double[1]();
    state->lgar_bmi_params.num_cells_temp = 1;
  }

  if (!is_ponded_depth_max_cm_set) {
    state->lgar_bmi_params.ponded_depth_max_cm = 0.0; // default maximum ponded depth is set to zero (i.e. no surface ponding)
  }


  state->lgar_bmi_params.forcing_interval = int(state->lgar_bmi_params.forcing_resolution_h/state->lgar_bmi_params.timestep_h+1.0e-08); // add 1.0e-08 to prevent truncation error

  // initialize frozen factor array to 1.
  state->lgar_bmi_params.frozen_factor = new double[state->lgar_bmi_params.num_layers+1];
  for (int i=0; i <= state->lgar_bmi_params.num_layers; i++)
    state->lgar_bmi_params.frozen_factor[i] = 1.0;

  InitializeWettingFronts(state->lgar_bmi_params.num_layers, state->lgar_bmi_params.initial_psi_cm,
			  state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.cum_layer_thickness_cm,
			  state->lgar_bmi_params.frozen_factor, state->soil_properties);

  if (verbosity.compare("none") != 0) {
    std::cerr<<"--- Initial state/conditions --- \n";
    listPrint();
    std::cerr<<"          *****         \n";
  }

  // initial mass in the system
  state->lgar_mass_balance.volstart_cm = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm);

  state->lgar_bmi_params.ponded_depth_cm = 0.0; // initially we start with a dry surface (no surface ponding)
  state->lgar_bmi_params.nint = 120; // hacked, not needed to be an input option
  state->lgar_bmi_params.num_wetting_fronts = state->lgar_bmi_params.num_layers;

  assert (state->lgar_bmi_params.num_layers == listLength());

  if (verbosity.compare("high") == 0) {
    std::cerr<<"Initial ponded depth is set to zero. \n";
    std::cerr<<"No. of spatial intervals used in trapezoidal integration to compute G : "<<state->lgar_bmi_params.nint<<"\n";
  }

  state->lgar_bmi_input_params = new lgar_bmi_input_parameters;
  state->lgar_bmi_params.time_s = 0.0;
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
extern void InitializeWettingFronts(int num_layers, double initial_psi_cm, int *layer_soil_type, double *cum_layer_thickness_cm,
				    double *frozen_factor, struct soil_properties_ *soil_properties)
{
  int soil;
  int front=0;
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
    bottom_flag=true;  // all initial wetting fronts are in contact with the bottom of the layer they exist in
    // NOTE: The listInsertFront function does lots of stuff.
    current = listInsertFront(cum_layer_thickness_cm[layer],theta_init,front,layer,bottom_flag);
    current->psi_cm = initial_psi_cm;
    Se = calc_Se_from_theta(current->theta,soil_properties[soil].theta_e,soil_properties[soil].theta_r);

    Ksat_cm_per_h = frozen_factor[layer] * soil_properties[soil].Ksat_cm_per_h;
    current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h , soil_properties[soil].vg_m);  // cm/s

  }

}

// ##################################################################################
/*
  Reads 1D data from the config file
  - used for reading soil discretization (1D)
  - used for reading layers depth from the surface if model `layered` is chosen
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


    assert (layer_temp > 100.0); // just a check to ensure the while loop executes at least once
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
  double volstart  = state->lgar_mass_balance.volstart_cm;
  double volprecip = state->lgar_mass_balance.volprecip_cm;
  double volrunoff = state->lgar_mass_balance.volrunoff_cm;
  double volAET    = state->lgar_mass_balance.volAET_cm;
  double volPET    = state->lgar_mass_balance.volPET_cm;
  double volon     = state->lgar_mass_balance.volon_cm;
  double volin     = state->lgar_mass_balance.volin_cm;
  double volrech   = state->lgar_mass_balance.volrech_cm;
  double volend    = state->lgar_mass_balance.volend_cm;
  double volrunoff_giuh = state->lgar_mass_balance.volrunoff_giuh_cm;
  double volend_giuh_cm = 0.0;
  double total_Q_cm     = state->lgar_mass_balance.volQ_cm;

  //check if the giuh queue have some water left at the end of simulaiton; needs to be included in the global mass balance
  // hold on; this is probably not needed as we have volrunoff in the balance; revist AJK
  for(int i=1; i <= state->lgar_bmi_params.num_giuh_ordinates; i++)
    volend_giuh_cm += giuh_runoff_queue_cm[i];


  double global_error_cm = volstart + volprecip - volrunoff - volAET - volon - volrech - volend;


 printf("\n********************************************************* \n");
 printf("-------------------- Simulation Summary ----------------- \n");
 //printf("Time (sec)                 = %6.10f \n", elapsed);
 printf("------------------------ Mass balance ------------------- \n");
 printf("Initial water in soil    = %14.10f cm\n", volstart);
 printf("Total precipitation      = %14.10f cm\n", volprecip);
 printf("Total infiltration       = %14.10f cm\n", volin);
 printf("Final water in soil      = %14.10f cm\n", volend);
 printf("Surface ponded water     = %14.10f cm\n", volon);
 printf("Surface runoff           = %14.10f cm\n", volrunoff);
 printf("GIUH runoff              = %14.10f cm\n", volrunoff_giuh);
 printf("Total percolation        = %14.10f cm\n", volrech);
 printf("Total AET                = %14.10f cm\n", volAET);
 printf("Total PET                = %14.10f cm\n", volPET);
 printf("Total discharge (Q)      = %14.10f cm\n", total_Q_cm);
 printf("Global balance           =   %.6e cm\n", global_error_cm);

}

// ############################################################################################
/*
 finds the wetting front that corresponds to psi (head) value closest to zero
 (i.e., saturation in terms of psi). This is the wetting front that experiences infiltration
 and actual ET based on precipitatona and PET, respectively. For example, the actual ET
 is extracted from this wetting front plus the wetting fronts above it.
 Note: the free_drainage name came from its python version, which is probably not the correct name.
 */
// ############################################################################################
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
extern void lgar_move_wetting_fronts(double timestep_h, double *volin_cm, int wf_free_drainage_demand,
				     double old_mass, int num_layers, double *AET_demand_cm, double *cum_layer_thickness_cm,
				     int *soil_type, double *frozen_factor, struct soil_properties_ *soil_properties)
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

  double column_depth = cum_layer_thickness_cm[num_layers];

  previous = head;
  double theta_e,theta_r;
  double vg_a, vg_m, vg_n;
  int layer_num, soil_num;

  int number_of_wetting_fronts = listLength();

  current = head;

  int last_wetting_front_index = number_of_wetting_fronts;
  int layer_num_above, layer_num_below;

  double precip_mass_to_add = (*volin_cm); // water to be added to the soil

  double bottom_boundary_flux_cm = 0.0; // water leaving the system through the bottom boundary

  *volin_cm = 0.0; // assuming that all the water can fit in, if not then re-assign the left over water at the end

  /* ************************************************************ */
  // main loop advancing all wetting fronts and doing the mass balance
  // loop goes over deepest to top most wetting front
  // wf denotes wetting front

  for (int wf = number_of_wetting_fronts; wf != 0; wf--) {

    if (verbosity.compare("high") == 0) {
      printf("Moving |******** Wetting Front = %d *********| \n", wf);
    }

    if (wf == 1 && number_of_wetting_fronts >0) {
      current = listFindFront(wf, NULL);
      next = listFindFront(wf+1, NULL);
      previous = NULL;

      current_old = listFindFront(wf, state_previous);
      next_old = listFindFront(wf+1, state_previous);
    }
    else if (wf < number_of_wetting_fronts) {
      current = listFindFront(wf, NULL);
      next = listFindFront(wf+1, NULL);
      previous = listFindFront(wf-1, NULL);

      current_old = listFindFront(wf, state_previous);
      next_old = listFindFront(wf+1, state_previous);
    }
    else if (wf == number_of_wetting_fronts) {
      current = listFindFront(wf, NULL);
      next = NULL;
      previous = listFindFront(wf-1, NULL);

      current_old = listFindFront(wf, state_previous);
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
       listPrint();
    }

    double free_drainage_demand = 0.0;
    double actual_ET_demand = *AET_demand_cm;

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
      //double psi_cm_below_old = 0.0;

      double psi_cm = current->psi_cm;
      //double psi_cm_below = 0.0;

      // mass = delta(depth) * delta(theta)
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
	delta_thetas[k] = theta_below;
	delta_thickness[k] = layer_thickness;
      }

      delta_thetas[layer_num] = 0.0;
      delta_thickness[layer_num] = current->depth_cm - cum_layer_thickness_cm[layer_num-1];

      double free_drainage_demand = 0;

      if (wf_free_drainage_demand == wf)
	prior_mass += precip_mass_to_add - (free_drainage_demand + actual_ET_demand);

      // theta mass balance computes new theta that conserves the mass; new theta is assigned to the current wetting front
      double theta_new = lgar_theta_mass_balance(layer_num, soil_num, psi_cm, new_mass, prior_mass,
						 delta_thetas, delta_thickness, soil_type, soil_properties);

      current->theta = fmin(theta_new, theta_e);

      double Se = calc_Se_from_theta(current->theta,theta_e,theta_r);
      current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);

      /* note: theta and psi of the current wetting front are updated here based on the wetting front's mass balance,
	 upper wetting fronts will be updated later in the lgar_merge_ module (the place where all state variables are updated
	 before proceeding to the next timestep. */

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

	double free_drainage_demand = 0;
	// prior mass = mass contained in the current old wetting front
	double prior_mass = current_old->depth_cm * (current_old->theta -  next_old->theta);

	if (wf_free_drainage_demand == wf)
	  prior_mass += precip_mass_to_add - (free_drainage_demand + actual_ET_demand);

	current->depth_cm += current->dzdt_cm_per_h * timestep_h;

	/* condition to bound the wetting front depth, if depth of a wf, at this timestep,
	   gets greater than the domain depth, it will be merge anyway as it is passing
	   the layer depth */
	current->depth_cm = fmin(current->depth_cm, column_depth);

	if (current->dzdt_cm_per_h == 0.0 && current->to_bottom == FALSE) // a new front was just created, so don't update it.
	  current->theta = current->theta;
	else
	  current->theta = fmin(theta_e, prior_mass/current->depth_cm + next->theta);

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

	  - LGAR paper (currently under review) has a better description, using diagrams, of the mass balance of wetting fronts
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

	double free_drainage_demand = 0;

	if (wf_free_drainage_demand == wf)
	  prior_mass += precip_mass_to_add - (free_drainage_demand + actual_ET_demand);


	// theta mass balance computes new theta that conserves the mass; new theta is assigned to the current wetting front
	double theta_new = lgar_theta_mass_balance(layer_num, soil_num, psi_cm, new_mass, prior_mass,
						   delta_thetas, delta_thickness, soil_type, soil_properties);

	current->theta = fmin(theta_new, theta_e);

      }


      double Se = calc_Se_from_theta(current->theta,theta_e,theta_r);
      current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
    }


    // if f_p (predicted infiltration) causes theta > theta_e, mass correction is needed.
    // depth of the wetting front is increased to close the mass balance when theta > theta_e.
    // l == 1 is the last iteration (top most wetting front), so do a check on the mass balance)
    // this part should be moved out of here to a subroutine; add a call to that subroutine
    if (wf == 1) {

      int soil_num_k1  = soil_type[wf_free_drainage_demand];
      double theta_e_k1 = soil_properties[soil_num_k1].theta_e;

      struct wetting_front *wf_free_drainage = listFindFront(wf_free_drainage_demand,NULL);

      double mass_timestep = (old_mass + precip_mass_to_add) - (actual_ET_demand + free_drainage_demand);

      assert (old_mass > 0.0);

      if (wf_free_drainage->theta == theta_e_k1) {

	double current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm);

	double mass_balance_error = fabs(current_mass - mass_timestep); // mass error

	double factor = 1.0;
	bool switched = false;
	double tolerance = 1e-12;

	// check if the difference is less than the tolerance
	if (mass_balance_error <= tolerance) {
	  // return current_mass;
	}

	double depth_new = wf_free_drainage->depth_cm;

	// loop to adjust the depth for mass balance
	while (mass_balance_error > tolerance) {

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

	  wf_free_drainage->depth_cm = depth_new;

	  current_mass = lgar_calc_mass_bal(cum_layer_thickness_cm);
	  mass_balance_error = fabs(current_mass - mass_timestep);

	}

      }
    }

  
  }
  
  /*******************************************************************/
  // end of the for loop
  /*******************************************************************/

  if (verbosity.compare("high") == 0) {
    printf("State after moving but before merging wetting fronts...\n");
    listPrint();
  }

  // **************************** MERGING WETTING FRONT ********************************
  
  /* In the python version of LGAR, wetting front merging, layer boundary crossing, and lower boundary crossing
     all occur in a loop that happens after wetting fronts have been moved. This prevents the model from crashing,
     as there are rare but possible cases where multiple merging / layer boundary crossing events will happen in
     the same time step. For example, if two wetting fronts cross a layer boundary in the same time step, it will
     be necessary for merging to occur before layer boundary crossing. So, LGAR-C now approaches merging in the
     same way as in LGAR-Py, where wetting fronts are moved, then a new loop does merging for all wetting fronts,
     then a new loop does layer boundary corssing for all wetting fronts, then a new loop does merging again for
     all wetting fronts, and then a new loop does lower boundary crossing for all wetting fronts. this is a thorough
     way to deal with these scenarios. */

  //merge
  lgar_merge_wetting_fronts(soil_type, frozen_factor, soil_properties);
  
  //cross
  lgar_wetting_fronts_cross_layer_boundary(num_layers, cum_layer_thickness_cm, soil_type, frozen_factor,
					   soil_properties);
  
  //merge
  lgar_merge_wetting_fronts(soil_type, frozen_factor, soil_properties);
  
  //lower bound
  bottom_boundary_flux_cm += lgar_wetting_front_cross_domain_boundary(cum_layer_thickness_cm, soil_type,
								      frozen_factor, soil_properties);


  *volin_cm = bottom_boundary_flux_cm;

  // reset current to head to fix any mass balance issues and dry-over-wet wetting fronts conditions
  double mass_change = 0.0;

  lgar_fix_dry_over_wet_fronts(&mass_change, cum_layer_thickness_cm, soil_type, soil_properties);

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


  /***********************************************/
  // make sure all psi values are updated
  current = head;

  for (int wf=1; wf != listLength(); wf++) {

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
    current = current->next;

  }

  if (verbosity.compare("high") == 0) {
    printf("Moving/merging wetting fronts done... \n");
  }


  //Just a check to make sure that, when there is only 1 layer, than the existing wetting front is at the correct depth.
  //This might have been fixed with other debugging related to scenarios with just 1 layer where the wetting front is completely satruated. Not sure this is necessary.
  if (listLength()==1){
    if (current->depth_cm != cum_layer_thickness_cm[1]){
      current->depth_cm = cum_layer_thickness_cm[1];
    }
  }


}


// ############################################################################################
/*
  the function merges wetting fronts; called from lgar_move_wetting_fronts.
*/
// ############################################################################################

extern void lgar_merge_wetting_fronts(int *soil_type, double *frozen_factor, struct soil_properties_ *soil_properties)
{
  struct wetting_front *current;
  struct wetting_front *next;
  struct wetting_front *next_to_next;
  current = head; 

  if (verbosity.compare("high") == 0) {
    printf("Merging wetting fronts... \n");
  }

  // local variables
  double theta_e,theta_r;
  double vg_a, vg_m, vg_n;
  double Se, Ksat_cm_per_h;
  int layer_num, soil_num;
    
  for (int wf=1; wf != listLength(); wf++) {
    
    if (verbosity.compare("high") == 0) {
      printf("Merge | ********* Wetting Front = %d *********\n", wf);
    }


    next = current->next;
    next_to_next = current->next->next;

    // case : wetting front passing another wetting front within a layer
    /**********************************************************/
    // 'current->depth_cm > next->depth_cm' ensures that merging is needed
    // 'current->layer_num == next->layer_num' ensures wetting fronts are in the same layer
    // '!next->to_bottom' ensures that the next wetting front is not the deepest wetting front in the layer
    if ( (current->depth_cm > next->depth_cm) && (current->layer_num == next->layer_num) && !next->to_bottom) {
      
      double current_mass_this_layer = current->depth_cm * (current->theta - next->theta) + next->depth_cm*(next->theta - next_to_next->theta);
      current->depth_cm = current_mass_this_layer / (current->theta - next_to_next->theta);
      
      layer_num = current->layer_num;
      soil_num  = soil_type[layer_num];
      theta_e   = soil_properties[soil_num].theta_e;
      theta_r   = soil_properties[soil_num].theta_r;
      vg_a      = soil_properties[soil_num].vg_alpha_per_cm;
      vg_m      = soil_properties[soil_num].vg_m;
      vg_n      = soil_properties[soil_num].vg_n;
      Se        = calc_Se_from_theta(current->theta,theta_e,theta_r);

      Ksat_cm_per_h  = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num];

      current->psi_cm     = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
      current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);
      
      if (verbosity.compare("high") == 0) {
        printf ("Deleting wetting front (before)... \n");
        listPrint();
      }
      
      listDeleteFront(next->front_num);
      
      if (verbosity.compare("high") == 0) {
        printf ("Deleting wetting front (after) ... \n");
        listPrint();
      }
    }
    
    current = current->next;
  }

  if (verbosity.compare("high") == 0) {
    printf("State after merging wetting fronts...\n");
    listPrint();
  }

}



// ############################################################################################
/*
  the function lets wetting fronts of a sufficient depth cross layer boundaries; called from lgar_move_wetting_fronts.
*/
// ############################################################################################

extern void lgar_wetting_fronts_cross_layer_boundary(int num_layers,
					double* cum_layer_thickness_cm, int *soil_type,
					double *frozen_factor, struct soil_properties_ *soil_properties)
{
  struct wetting_front *current;
  struct wetting_front *next;
  struct wetting_front *next_to_next;
  current = head; 

  if (verbosity.compare("high") == 0) {
    printf("Layer boundary crossing... \n");
  }
  for (int wf=1; wf != listLength(); wf++) {
    
    if (verbosity.compare("high") == 0) {
      printf("Boundary Crossing | ******* Wetting Front = %d ****** \n", wf);
    }
    
    // local variables
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
    double Ksat_cm_per_h  = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num]; //PTL addition to make K_cm_per_h for this conditon to be correct
    
    next = current->next;
    next_to_next = current->next->next;

    if (current->depth_cm > cum_layer_thickness_cm[layer_num] && (next->depth_cm == cum_layer_thickness_cm[layer_num])
	&& (layer_num!=num_layers) ) {
      
      double current_theta = fmin(theta_e, current->theta);
      double overshot_depth = current->depth_cm - next->depth_cm;
      double next_theta_e   = soil_properties[soil_num+1].theta_e;
      double next_theta_r   = soil_properties[soil_num+1].theta_r;
      double next_vg_a      = soil_properties[soil_num+1].vg_alpha_per_cm;
      double next_vg_m      = soil_properties[soil_num+1].vg_m;
      double next_vg_n      = soil_properties[soil_num+1].vg_n;
      //double next_Ksat_cm_per_h  = soil_properties[soil_num+1].Ksat_cm_per_h * frozen_factor[current->layer_num]; 

      double Se = calc_Se_from_theta(current->theta,theta_e, theta_r);
      current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);

      current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m);
      
      // current psi with van Gunechton properties of the next layer to get new theta
      double theta_new = calc_theta_from_h(current->psi_cm, next_vg_a, next_vg_m, next_vg_n, next_theta_e, next_theta_r);

      double mbal_correction = overshot_depth * (current_theta - next->theta);
      double mbal_Z_correction = mbal_correction / (theta_new - next_to_next->theta); // this is the new wetting front depth

      double depth_new = cum_layer_thickness_cm[layer_num] + mbal_Z_correction; // this is the new wetting front absolute depth


      current->depth_cm = cum_layer_thickness_cm[layer_num];
      
      next->theta = theta_new;
      next->psi_cm = current->psi_cm;
      next->depth_cm = depth_new;
      next->layer_num = layer_num + 1;
      next->dzdt_cm_per_h = current->dzdt_cm_per_h;
      current->dzdt_cm_per_h = 0;
      current->to_bottom = TRUE;
      next->to_bottom = FALSE;
      
    }
    
    if (verbosity.compare("high") == 0) {
      printf("State after wetting fronts cross layer boundary...\n");
      listPrint();
    }
    current = current->next;
  }
  
}


// ############################################################################################
/*
  the function lets wetting fronts of a sufficient depth interact with the lower boundary; called from lgar_move_wetting_fronts.
*/
// ############################################################################################

extern double lgar_wetting_front_cross_domain_boundary(double* cum_layer_thickness_cm, int *soil_type,
						       double *frozen_factor, struct soil_properties_ *soil_properties)
{
  struct wetting_front *current;
  struct wetting_front *next;
  struct wetting_front *next_to_next;
  current = head; 
  double bottom_flux_cm = 0.0;

  if (verbosity.compare("high") == 0) {
    printf("Domain boundary crossing (bottom flux calc.) \n");
  }

  // local variables
  double theta_e,theta_r;
  double vg_a, vg_m, vg_n;
  double bottom_flux_cm_temp;
  int layer_num, soil_num;
    
  for (int wf=1; wf != listLength(); wf++) {

    if (verbosity.compare("high") == 0) {
      printf("Domain boundary crossing | ***** Wetting Front = %d ****** \n", wf);
    }
    
    bottom_flux_cm_temp = 0.0;
    
    
    layer_num   = current->layer_num;
    soil_num    = soil_type[layer_num];
    
    next = current->next;
    next_to_next = current->next->next;

    // case : wetting front is the deepest one in the last layer (most deepested wetting front in the domain)
    /**********************************************************/
    if (next_to_next == NULL && current->depth_cm > cum_layer_thickness_cm[layer_num]) {
      
      //  this is the water leaving the system through the bottom of the soil
      bottom_flux_cm_temp = (current->theta - next->theta) *  (current->depth_cm - next->depth_cm);
      theta_e   = soil_properties[soil_num].theta_e;
      theta_r   = soil_properties[soil_num].theta_r;
      vg_a      = soil_properties[soil_num].vg_alpha_per_cm;
      vg_m      = soil_properties[soil_num].vg_m;
      vg_n      = soil_properties[soil_num].vg_n;
      double Ksat_cm_per_h  = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[current->layer_num];

      next->theta = current->theta;
      double Se_k = calc_Se_from_theta(current->theta,theta_e,theta_r);
      next->psi_cm = calc_h_from_Se(Se_k, vg_a, vg_m, vg_n);
      next->K_cm_per_h = calc_K_from_Se(Se_k, Ksat_cm_per_h, vg_m);
      listDeleteFront(current->front_num);
      
    }
    
    if (verbosity.compare("high") == 0) {
      printf("State after lowest wetting front contributes to flux through the bottom boundary...\n");
      listPrint();
    }
    
    bottom_flux_cm += bottom_flux_cm_temp;
    current = current->next;
    
    if (verbosity.compare("high") == 0){
      printf("Bottom boundary flux = %lf \n",bottom_flux_cm); }
    if (bottom_flux_cm_temp!=0){//PTL: perhaps due to potential problems with precision error, this should just be replaces with a flag that gets toggled when bottom_boundary_flux_cm is changed from 0 to nonzero
      break;
    }
  }

  return bottom_flux_cm;
  
}


// ############################################################################################
/* The function handles situation of dry over wet wetting fronts
  mainly happen when AET extracts more water from the upper wetting front
  and the front gets drier than the lower wetting front */
// ############################################################################################
extern void lgar_fix_dry_over_wet_fronts(double *mass_change, double* cum_layer_thickness_cm, int *soil_type,
					 struct soil_properties_ *soil_properties)
{
  if (verbosity.compare("high") == 0) {
    printf("Fix Dry over Wet Wetting Front... \n");
  }

  struct wetting_front *current;
  struct wetting_front *next;
  current = head;
  next = current->next;

  for (int l=1; l <= listLength(); l++) {

    if (next != NULL) {
      // this part fixes case of upper theta less than lower theta due to AET extraction
      // also handles the case when the current and next wetting fronts have the same theta
      // and are within the same layer
      /***************************************************/

      if ( (current->theta <= next->theta) && (current->layer_num == next->layer_num) ) {
	int layer_num_k = current->layer_num;
	double mass_before = lgar_calc_mass_bal(cum_layer_thickness_cm);

	current = listDeleteFront(current->front_num);
	
	// if the dry wetting front is the most surficial then simply track the mass change
	// due to the deletion of the wetting front;
	// this needs to be revised
	if (layer_num_k > 1) {
	  int soil_num_k    = soil_type[current->layer_num];
	  double theta_e_k  = soil_properties[soil_num_k].theta_e;
	  double theta_r_k  = soil_properties[soil_num_k].theta_r;
	  double vg_a_k     = soil_properties[soil_num_k].vg_alpha_per_cm;
	  double vg_m_k     = soil_properties[soil_num_k].vg_m;
	  double vg_n_k     = soil_properties[soil_num_k].vg_n;
	  double Se_k       = calc_Se_from_theta(current->theta,theta_e_k,theta_r_k);

	  // now this is the wetting front that was below the dry wetting front
	  current->psi_cm = calc_h_from_Se(Se_k, vg_a_k, vg_m_k, vg_n_k);

	  struct wetting_front *current_local = head;

	  // update psi and theta for all wetting fronts above the current wetting front
	  while (current_local->layer_num < layer_num_k) {
	    int soil_num_k1 = soil_type[current_local->layer_num];
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

	double mass_after = lgar_calc_mass_bal(cum_layer_thickness_cm);
	*mass_change += fabs(mass_after - mass_before);

	/* note: mass_before is less when we have wetter front over drier front condition,
	   however, lgar_calc_mass_bal returns mass_before > mass_after due to fabs(theta_current - theta_next);
	   for mass_before the functions compuates more than the actual mass; removing fabs in that function
	   might be one option, but for now we are adding fabs to mass_change to make sure we added extra water
	   back to AET after deleting the drier front */

      }

      current = current->next;

      if (current == NULL)
	next = NULL;
      else
	next = current->next;

    }
  }

}

// ############################################################################################
/* The module computes the potential infiltration capacity, fp (in the lgar manuscript),
   potential infiltration capacity = the maximum amount of water that can be inserted into
   the soil depending on the availability of water.
   this module is called when a new superficial wetting front is not created
   in the current timestep, that is precipitation in the current and previous
   timesteps was greater than zero */
// ############################################################################################
extern double lgar_insert_water(bool use_closed_form_G, int nint, double timestep_h, double *ponded_depth_cm, double *volin_this_timestep,
				double precip_timestep_cm, int wf_free_drainage_demand,
			        int num_layers, double ponded_depth_max_cm, int *soil_type,
				double *cum_layer_thickness_cm, double *frozen_factor,
				struct soil_properties_ *soil_properties)
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
  current_free_drainage      = listFindFront(wf_that_supplies_free_drainage_demand, NULL);
  current_free_drainage_next = listFindFront(wf_that_supplies_free_drainage_demand+1, NULL);

  int number_of_wetting_fronts = listLength();

  //int last_wetting_front_index = number_of_wetting_fronts;
  int layer_num_fp = current_free_drainage->layer_num;


  double Geff;

  if (number_of_wetting_fronts == num_layers) {
    Geff = 0.0; // i.e., case of no capillary suction, dz/dt is also zero for all wetting fronts
  }
  else {

    //double theta = current_free_drainage->theta;
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

    //Se = calc_Se_from_theta(theta,theta_e,theta_r);
    //psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);

    Geff = calc_Geff(use_closed_form_G, theta_below, theta_e, theta_e, theta_r, vg_a, vg_n, vg_m, h_min_cm, Ksat_cm_per_h, nint, lambda, bc_psib_cm); 

  }

  // if the free_drainage wetting front is the top most, then the potential infiltration capacity has the following simple form
  if (layer_num_fp == 1) {
      f_p = Ksat_cm_per_h * (1 + (Geff + h_p)/current_free_drainage->depth_cm);
  }
  else {
    // point here to the equation in lgar paper once published
    double bottom_sum = (current_free_drainage->depth_cm - cum_layer_thickness_cm[layer_num_fp-1])/Ksat_cm_per_h;
    //printf("##################################################################Inside f_p calc. bottom_sum %.17g\n", bottom_sum);

    for (int k = 1; k < layer_num_fp; k++) {
      int soil_num_k = soil_type[layer_num_fp-k];
      double Ksat_cm_per_h_k = soil_properties[soil_num_k].Ksat_cm_per_h * frozen_factor[layer_num_fp - k];

      bottom_sum += (cum_layer_thickness_cm[layer_num_fp - k] - cum_layer_thickness_cm[layer_num_fp - (k+1)])/ Ksat_cm_per_h_k;
    }

    f_p = (current_free_drainage->depth_cm / bottom_sum) + ((Geff + h_p)*Ksat_cm_per_h/(current_free_drainage->depth_cm)); //Geff + h_p

  }

  // checkpoint # AJ
  int soil_num_k = soil_type[head->layer_num];
  double theta_e1 = soil_properties[soil_num_k].theta_e; // saturated theta of top layer

  // if free drainge has to be included, which currently we don't, then the following will be set to hydraulic conductivity
  // of the deeepest layer
  if ((layer_num_fp == num_layers) && (current_free_drainage->theta == theta_e1) && (num_layers == number_of_wetting_fronts))
    f_p = 0.0;

  double ponded_depth_temp = *ponded_depth_cm;

  double free_drainage_demand = 0;

  // 'if' condition is not needed ... AJ
  if ((layer_num_fp==num_layers) && (num_layers == number_of_wetting_fronts))
    ponded_depth_temp = *ponded_depth_cm - f_p * timestep_h - free_drainage_demand*0;
  else
    ponded_depth_temp = *ponded_depth_cm - f_p * timestep_h - free_drainage_demand*0;

  ponded_depth_temp   = fmax(ponded_depth_temp, 0.0);

  double fp_cm = f_p * timestep_h + free_drainage_demand/timestep_h; // infiltration in cm

  if (ponded_depth_max_cm > 0.0 ) {

    if (ponded_depth_temp < ponded_depth_max_cm) {
      runoff = 0.0;
      *volin_this_timestep = fmin(*ponded_depth_cm, fp_cm); //PTL: does this code account for the case where volin_this_timestep can not all infiltrate?
      *ponded_depth_cm     = *ponded_depth_cm - *volin_this_timestep;
      return runoff;
    }
    else if (ponded_depth_temp > ponded_depth_max_cm ) {
      runoff = ponded_depth_temp - ponded_depth_max_cm;
      *ponded_depth_cm     = ponded_depth_max_cm;
      *volin_this_timestep = fp_cm;

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

  return runoff;
}

// ######################################################################################
/* This subroutine is called iff there is no surfacial front, it creates a new front and
   inserts ponded depth, and will return some amount if can't fit all water */
// ######################################################################################
extern void lgar_create_surfacial_front(double *ponded_depth_cm, double *volin, double dry_depth,
					double theta1, int *soil_type, double *cum_layer_thickness_cm,
					double *frozen_factor, struct soil_properties_ *soil_properties)
{
  // into the soil.  Note ponded_depth_cm is a pointer.   Access its value as (*ponded_depth_cm).

  // local vars
  double theta_e,Se,theta_r;
  double delta_theta;
  double vg_alpha_per_cm, vg_m, vg_n, Ksat_cm_per_h;

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

  double theta_new;

  if(dry_depth * delta_theta > (*ponded_depth_cm))  // all the ponded depth enters the soil
    {
      *volin = *ponded_depth_cm;
      theta_new = fmin(theta1 + (*ponded_depth_cm) /dry_depth, theta_e);
      listInsertFirst(dry_depth, theta_new, front_num, layer_num, to_bottom);
      *ponded_depth_cm = 0.0;
      //hp_cm =0.0;
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
      //hp_cm = *ponded_depth_cm;
    }

  current = head;  // must do this again because listInsertFirst() created a new *head
  vg_alpha_per_cm    = soil_properties[soil_num].vg_alpha_per_cm;
  vg_m               = soil_properties[soil_num].vg_m;
  vg_n               = soil_properties[soil_num].vg_n;
  Ksat_cm_per_h      = soil_properties[soil_num].Ksat_cm_per_h * frozen_factor[layer_num];

  Se = calc_Se_from_theta(theta_new,theta_e,theta_r);
  current->psi_cm = calc_h_from_Se(Se, vg_alpha_per_cm , vg_m, vg_n);

  current->K_cm_per_h = calc_K_from_Se(Se, Ksat_cm_per_h, vg_m) * frozen_factor[layer_num]; // AJ - K_temp in python version for 1st layer

  current->dzdt_cm_per_h = 0.0; //for now assign 0 to dzdt as it will be computed/updated in lgar_dzdt_calc function

  return;

}

// ############################################################################################
/* This routine calculates the "dry depth" of a newly created wetting front in the top soil layer after
   a non-rainy period or a big increase in rainrate  on an unsaturated first layer.
   Note: Calculation of the initial depth of a new wetting front in the first layer uses the concept of "dry depth",
   described in the 2015 GARTO paper (Lai et al., An efficient and guaranteed stable numerical method ffor
   continuous modeling of infiltration and redistribution with a shallow dynamic water table). */
// ############################################################################################
extern double lgar_calc_dry_depth(bool use_closed_form_G, int nint, double timestep_h, double *delta_theta, int *soil_type,
				  double *cum_layer_thickness_cm, double *frozen_factor,
				  struct soil_properties_ *soil_properties)
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

  return dry_depth;

}

// ###########################################################################
/* function to calculate the amount of soil moisture (total mass of water)
   in the profile (cm) */
// ###########################################################################
double lgar_calc_mass_bal(double *cum_layer_thickness)
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
	if(next->layer_num == current->layer_num)
	  sum += (current->depth_cm - base_depth) * (current->theta - next->theta); // note no need for fabs() here otherwise we get more mass for the case dry-over-wet front within a layer
	else
	  sum += (current->depth_cm - base_depth) * current->theta;
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
   Open file to read in the van Genuchton parameters for standard soil types*/
// ############################################################################################
extern int lgar_read_vG_param_file(char const* vG_param_file_name, int num_soil_types, double wilting_point_psi_cm,
				    struct soil_properties_ *soil_properties)
{

  if (verbosity.compare("high") == 0) {
    std::cerr<<"Reading van Genuchton parameters files...\n";
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

  //for(soil=1;soil<=num_soil_types;soil++) {// read the num_soil_types lines with data
  while (fgets(jstr,255,in_vG_params_fptr) != NULL) {

    sscanf(jstr,"%s %lf %lf %lf %lf %lf %lf",soil_name,&theta_r,&theta_e,&vg_alpha_per_cm,&vg_n,&vg_m,&Ksat_cm_per_h);
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
    soil_properties[soil].vg_m            = vg_m;
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
extern void lgar_dzdt_calc(bool use_closed_form_G, int nint, double h_p, int *soil_type, double *cum_layer_thickness_cm, double *frozen_factor,
			  struct soil_properties_ *soil_properties)
{
  if (verbosity.compare("high") == 0) {
    std::cerr<<"Calculating dz/dt .... \n";
  }

  struct wetting_front* current;
  struct wetting_front* next;

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

  // make sure to use previous state values as current state is updated during the timestep (that's how it is done is Peter's python version)

  current = head;

  do {  // loop through the wetting fronts
    dzdt = 0.0;

    // copy structure elements into shorter variables names to increase readability
    // WETTING FRONT PROPERTIES
    layer_num    = current->layer_num;    // what layer the front is in
    K_cm_per_h   = current->K_cm_per_h;   // K(theta)

    if (K_cm_per_h <=0) {
      printf("K is zero: %d %lf \n", layer_num, K_cm_per_h);
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

    if(current->to_bottom == TRUE) { // checkpoint # AJ
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

    current->dzdt_cm_per_h = dzdt;

    current = current->next;  // point to the next link

  } while(current != NULL );   // putting conditional at end of do looop makes sure it executes at least once

}

// ############################################################################################
/* The function does mass balance for a wetting front to get an updated theta.
   The head (psi) value is iteratively altered until the error between prior mass and new mass
   is within a tolerance. */
// ############################################################################################
extern double lgar_theta_mass_balance(int layer_num, int soil_num, double psi_cm, double new_mass,
				      double prior_mass, double *delta_theta, double *delta_thickness,
				      int *soil_type, struct soil_properties_ *soil_properties)
{

  double psi_cm_loc = psi_cm; // location psi
  double delta_mass = fabs(new_mass - prior_mass); // mass different between the new and prior
  double tolerance = 1e-12;

  double factor = 1.0;
  bool switched = false; // flag that determines capillary head to be incremented or decremented

  double theta = 0; // this will be updated and returned


  // check if the difference is less than the tolerance
  if (delta_mass <= tolerance) {
    theta = calc_theta_from_h(psi_cm_loc, soil_properties[soil_num].vg_alpha_per_cm, soil_properties[soil_num].vg_m,
			      soil_properties[soil_num].vg_n,soil_properties[soil_num].theta_e,soil_properties[soil_num].theta_r);
    return theta;
  }

  // the loop increments/decrements the capillary head until mass difference between
  // the new and prior is within the tolerance
  while (delta_mass > tolerance) {

    if (new_mass>prior_mass) {
      psi_cm_loc += 0.1 * factor;
      switched = false;
    }
    else {
      if (!switched) {
	switched = true;
	factor = factor * 0.1;
  //factor = factor * 1;
      }
      double psi_cm_loc_temp = psi_cm_loc;
      psi_cm_loc -= 0.1 * factor;

      if (psi_cm_loc<0 && psi_cm_loc_temp!=0){// this is for the extremely rare case when iterative psi_cm_loc calculation temporarily yields a negative value
                                              // // and the actual answer for psi_cm_loc is nonzero. For example when a completely saturated wetting front with a tiny
                                              // // amount of ET should yield a resulting theta that is slightly below saturation.
        //abort();
        psi_cm_loc = psi_cm_loc_temp * 0.1;
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

      theta_layer = calc_theta_from_h(psi_cm_loc, soil_properties[soil_num_loc].vg_alpha_per_cm, soil_properties[soil_num_loc].vg_m,
				      soil_properties[soil_num_loc].vg_n,soil_properties[soil_num_loc].theta_e,
				      soil_properties[soil_num_loc].theta_r);

      mass_layers += delta_thickness[k] * (theta_layer - delta_theta[k]);
    }

    new_mass = mass_layers;
    delta_mass = fabs(new_mass - prior_mass);
    //delta_mass = new_mass - prior_mass; 
  }
  return theta;

}

#endif
