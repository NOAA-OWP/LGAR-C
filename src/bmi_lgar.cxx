#ifndef BMI_LGAR_CXX_INCLUDED
#define BMI_LGAR_CXX_INCLUDED


#include <stdio.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>
#include "../bmi/bmi.hxx"
#include "../include/bmi_lgar.hxx"
#include "../include/all.hxx"


// default verbosity is set to 'none' other option 'high' or 'low' needs to be specified in the config file
string verbosity="none";

// some macros for flux caching, these are not relevant if flux caching is disabled in the config file.
// precip intensity above which flux caching does not occur
#ifndef PRECIP_THRESHOLD_MM_PER_H
#define PRECIP_THRESHOLD_MM_PER_H ((double)1.0e-6)
#endif

// Small epsilon added to PET denominator to prevent division by 0, and in another case to check if there is significant PET
#ifndef PET_EPSILON
#define PET_EPSILON ((double)1.0e-8)
#endif

// Threshold for AET/PET ratio to indicate low enough AET for flux caching to start
#ifndef AET_PET_RATIO_THRESHOLD
#define AET_PET_RATIO_THRESHOLD ((double)0.75)
#endif

// while ponded head is greater than this, flux caching will not occur
#ifndef VOLON_TIMESTEP_THRESHOLD_CM
#define VOLON_TIMESTEP_THRESHOLD_CM ((double)1.0e-16)
#endif

// While the vadose zone lower boundary flux is greater than this, flux caching will not happen
#ifndef BOTTOM_BDY_FLUX_THRESHOLD_CM
#define BOTTOM_BDY_FLUX_THRESHOLD_CM ((double)1.0e-4)
#endif

// While the top most wetting front has a dzdt value greater than this, flux caching will not happen
#ifndef THRESHOLD_DZDT_CM_PER_H
#define THRESHOLD_DZDT_CM_PER_H ((double)1.0e0)
#endif

// Number of timesteps before flux cache reset
#ifndef NUM_TIMESTEPS_BEFORE_RESET_CACHE
#define NUM_TIMESTEPS_BEFORE_RESET_CACHE 24
#endif

// small epsillon that is used to determine if the difference between two quantities is 0 while avoiding machine precision errors
#define SMALL_EPS 1.E-12


/**
 * @brief Delete dynamic arrays allocated in Initialize() and held by this object
 * 
 */
BmiLGAR::~BmiLGAR(){
  delete [] giuh_ordinates;
  delete [] giuh_runoff_queue;
}

/* The `head` pointer stores the address in memory of the first member of the linked list containing
   all the wetting fronts. The contents of struct wetting_front are defined in "all.h" */

void BmiLGAR::
Initialize (std::string config_file)
{
  if (config_file.compare("") != 0 ) {
    this->state = new model_state;
    state->head = NULL;
    state->state_previous = NULL;
    lgar_initialize(config_file, state);
  }

  num_giuh_ordinates = state->lgar_bmi_params.num_giuh_ordinates;

  /* giuh ordinates are static and read in the lgar.cxx, and we need to have a copy of it to pass to
     giuh.cxx, so allocating/copying here*/

  giuh_ordinates = new double[num_giuh_ordinates];
  giuh_runoff_queue = new double[num_giuh_ordinates+1];

  for (int i=0; i<num_giuh_ordinates;i++){
    giuh_ordinates[i] = state->lgar_bmi_params.giuh_ordinates[i+1]; // note lgar uses 1-indexing
  }

  for (int i=0; i<=num_giuh_ordinates;i++){
    giuh_runoff_queue[i] = 0.0;
  }

}

/**
 * @brief Allocate (or reallocate) storage for soil parameters
 * 
 */
void BmiLGAR::realloc_soil(){
  
  delete [] state->lgar_bmi_params.soil_depth_wetting_fronts;
  delete [] state->lgar_bmi_params.soil_moisture_wetting_fronts;

  state->lgar_bmi_params.soil_depth_wetting_fronts = new double[state->lgar_bmi_params.num_wetting_fronts];
  state->lgar_bmi_params.soil_moisture_wetting_fronts = new double[state->lgar_bmi_params.num_wetting_fronts];
}

/*
  This is the main function calling lgar subroutines for creating, moving, and merging wetting fronts.
  Calls to AET and mass balance module are also happening here
  If the model's timestep is smaller than the forcing's timestep then we take subtimesteps inside the subcycling loop
*/
void BmiLGAR::
Update()
{
  if (verbosity.compare("none") != 0) {
    std::cerr<<"---------------------------------------------------------\n";
    std::cerr<<"|****************** LASAM BMI Update... ******************|\n";
    std::cerr<<"---------------------------------------------------------\n";
  }

  double mm_to_cm = 0.1; // unit conversion
  double mm_to_m = 0.001;
  
  if (state->lgar_bmi_params.is_invalid_soil_type) {
    // add to mass balance accumulated variables
    state->lgar_mass_balance.volprecip_cm  += state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm;
    state->lgar_mass_balance.volin_cm       = 0.0;
    state->lgar_mass_balance.volon_cm       = 0.0;
    state->lgar_mass_balance.volend_cm      = state->lgar_mass_balance.volstart_cm;
    state->lgar_mass_balance.volCRend_cm    = state->lgar_mass_balance.volCRstart_cm;
    state->lgar_mass_balance.volAET_cm      = 0.0;
    state->lgar_mass_balance.volrech_cm     = 0.0;
    state->lgar_mass_balance.volrunoff_cm  += state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm;
    state->lgar_mass_balance.volQ_cm       += state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm;
    state->lgar_mass_balance.volQ_CR_cm     = 0.0;
    state->lgar_mass_balance.volPET_cm      = 0.0;
    state->lgar_mass_balance.volrunoff_giuh_cm  = 0.0;
    state->lgar_mass_balance.volchange_calib_cm = 0.0;

    // converted values, a struct local to the BMI and has bmi output variables
    bmi_unit_conv.mass_balance_m        = 0.0;
    bmi_unit_conv.volprecip_timestep_m  = state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_m;
    bmi_unit_conv.volin_timestep_m      = 0.0;
    bmi_unit_conv.volend_timestep_m     = 0.0;
    bmi_unit_conv.volCRend_timestep_m   = 0.0;
    bmi_unit_conv.volAET_timestep_m     = 0.0;
    bmi_unit_conv.volrech_timestep_m    = 0.0;
    bmi_unit_conv.volrunoff_timestep_m  = state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_m;
    bmi_unit_conv.volQ_timestep_m       = state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_m;
    bmi_unit_conv.volQ_CR_timestep_m    = 0.0;
    bmi_unit_conv.volPET_timestep_m     = 0.0;
    bmi_unit_conv.volrunoff_giuh_timestep_m = 0.0;

    return;
  }
  
  // if lasam is coupled to soil freeze-thaw, frozen fraction module is called
  if (state->lgar_bmi_params.sft_coupled)
    frozen_factor_hydraulic_conductivity(state->lgar_bmi_params);

  double volchange_calib_cm = 0.0;

  if(state->lgar_bmi_params.calib_params_flag) {
    volchange_calib_cm = update_calibratable_parameters(); // change in soil water volume due to calibratable parameters
    state->lgar_bmi_params.calib_params_flag = false;
  }


  // local variables for readibility
  int subcycles;
  int num_layers = state->lgar_bmi_params.num_layers;

  // local variables for a full timestep (i.e., timestep of the forcing data)
  // see 'struct lgar_mass_balance_variables' in all.hxx for full description of the variables
  double precip_timestep_cm   = 0.0;
  double PET_timestep_cm      = 0.0;
  double AET_timestep_cm      = 0.0;
  double volend_timestep_cm   = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head); // this should not be reset to 0.0 in the for loop
  double volCRend_timestep_cm = state->lgar_mass_balance.CR_fast_storage_cm + state->lgar_mass_balance.CR_slow_storage_cm;
  double volin_timestep_cm    = 0.0;
  double volon_timestep_cm    = state->lgar_mass_balance.volon_timestep_cm;
  double volrunoff_timestep_cm      = 0.0;
  double volrech_timestep_cm        = 0.0;
  double volrunoff_giuh_timestep_cm = 0.0;
  double volQ_timestep_cm           = 0.0;
  double volQ_CR_timestep_cm        = 0.0;
  
  // local variables for a subtimestep (i.e., timestep of the model)
  double precip_subtimestep_cm;
  double precip_subtimestep_cm_per_h;
  double PET_subtimestep_cm;
  double PET_subtimestep_cm_per_h;
  double ponded_depth_subtimestep_cm;
  double AET_subtimestep_cm;
  double volstart_subtimestep_cm;
  double volend_subtimestep_cm = volend_timestep_cm; // this should not be reset to 0.0 in the for loop
  double volin_subtimestep_cm;
  double volon_subtimestep_cm;
  double volrunoff_subtimestep_cm;
  double volrech_subtimestep_cm;
  double precip_previous_subtimestep_cm;
  double volCRstart_subtimestep_cm;
  double volCRend_subtimestep_cm = volCRend_timestep_cm;
  
  double subtimestep_h = state->lgar_bmi_params.timestep_h;
  int nint = state->lgar_bmi_params.nint;
  double wilting_point_psi_cm = state->lgar_bmi_params.wilting_point_psi_cm;
  double field_capacity_psi_cm = state->lgar_bmi_params.field_capacity_psi_cm;
  double a = state->lgar_bmi_params.a;
  double b = state->lgar_bmi_params.b;
  double frac_to_CR = state->lgar_bmi_params.frac_to_CR;
  double a_slow = state->lgar_bmi_params.a_slow;
  double b_slow = state->lgar_bmi_params.b_slow;
  double frac_slow = state->lgar_bmi_params.frac_slow;
  double spf_factor = state->lgar_bmi_params.spf_factor;
  bool use_closed_form_G = state->lgar_bmi_params.use_closed_form_G; 
  bool adaptive_timestep = state->lgar_bmi_params.adaptive_timestep;
  bool PET_affects_precip = state->lgar_bmi_params.PET_affects_precip;
  double mbal_tol = state->lgar_bmi_params.mbal_tol;

  // constant value used in the AET function
  double AET_thresh_Theta = 0.85;    // scaled soil moisture (0-1) above which AET=PET (fix later!)
  double AET_expon        = 1.0;     // exponent that allows curvature of the rising portion of the Budyko curve (fix later!)

  double ponded_depth_max_cm = state->lgar_bmi_params.ponded_depth_max_cm;

  if (verbosity.compare("high") == 0) {
    std::cerr<<"Pr  [cm/h] (timestep) = "<<state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm <<"\n";
    std::cerr<<"PET [cm/h] (timestep) = "<<state->lgar_bmi_input_params->PET_mm_per_h * mm_to_cm <<"\n"; 
  }

  assert (state->lgar_bmi_input_params->precipitation_mm_per_h >= 0.0);
  assert(state->lgar_bmi_input_params->PET_mm_per_h >=0.0);

  // adaptive time step is set 
  if (adaptive_timestep) {
    subtimestep_h = state->lgar_bmi_params.forcing_resolution_h;
    if (state->lgar_bmi_input_params->precipitation_mm_per_h > 10.0 || volon_timestep_cm > 0.0 ) {
      subtimestep_h = state->lgar_bmi_params.minimum_timestep_h;  //case where precip > 1 cm/h, or there is ponded head from the last time step
    }
    else if (state->lgar_bmi_input_params->precipitation_mm_per_h > 0.0) {
      subtimestep_h = state->lgar_bmi_params.minimum_timestep_h * 2.0;  //case where precip is less than 1 cm/h but greater than 0, and there is no ponded head 
    }
    subtimestep_h = fmin(subtimestep_h, state->lgar_bmi_params.forcing_resolution_h);  //just in case the user has specified a minimum time step that would make the subtimestep_h greater than the forcing resolution 
    state->lgar_bmi_params.timestep_h = subtimestep_h;
  }

  bool caching_at_start = state->lgar_mass_balance.cache_fluxes;
  bool switch_caching = FALSE;

  double wf_free_drain_dzdt = 0.0;

  if (!caching_at_start){
    struct wetting_front *wf_free_drainage_for_cache = listFindFront(wetting_front_free_drainage(state->head), state->head, NULL);
    wf_free_drain_dzdt = wf_free_drainage_for_cache->dzdt_cm_per_h;
  }

  if (state->lgar_bmi_params.allow_flux_caching){
    //The idea here is that, during dry periods, AET will become a small fraction of PET and wetting fronts will be very slow moving. In these cases, it is not necessary to compute fluxes for every time step.
    //To save on runtime, and if allow_flux_caching is set to true in the config file, we simply cache computed fluxes to be used for subsequent time steps rather than recomputing them.
    //The current implementation is to not move the wetting fronts under these conditions, but then move them more rapidly once it is time to calculate fluxes again. Also, during these periods, PET will be 0 but made to be larger to conserve mass when it is time to recalculate fluxes.
    //Recharge is copied from the last time step. If free drainage is enabled, this will technically result in a small mass balance error. 
    //There is very little change to the simulation when this is enabled.
    //Note that for NextGen models, it is ultimately desirable that there is technically output for every hour, so simply relaxing the adaptive time step to be coarser than 1 hour isn't the best solution
    if (subtimestep_h == state->lgar_bmi_params.forcing_resolution_h){
      if ( (state->lgar_bmi_input_params->precipitation_mm_per_h < PRECIP_THRESHOLD_MM_PER_H) && ( (state->lgar_mass_balance.previous_AET / (state->lgar_mass_balance.previous_PET + PET_EPSILON )) < AET_PET_RATIO_THRESHOLD) && (state->lgar_bmi_params.cache_count!=NUM_TIMESTEPS_BEFORE_RESET_CACHE) && 
           (volon_timestep_cm<VOLON_TIMESTEP_THRESHOLD_CM) && (state->lgar_mass_balance.previous_recharge<BOTTOM_BDY_FLUX_THRESHOLD_CM) ){
        state->lgar_mass_balance.cache_fluxes = TRUE;
      }
      if ( wf_free_drain_dzdt>THRESHOLD_DZDT_CM_PER_H ){
        state->lgar_mass_balance.cache_fluxes = FALSE;
      }
      if (volon_timestep_cm>=VOLON_TIMESTEP_THRESHOLD_CM){
        state->lgar_mass_balance.cache_fluxes = FALSE;
      }
      if ( (state->lgar_bmi_params.cache_count==NUM_TIMESTEPS_BEFORE_RESET_CACHE ) || (state->lgar_bmi_input_params->precipitation_mm_per_h >= PRECIP_THRESHOLD_MM_PER_H) ){
        state->lgar_mass_balance.cache_fluxes = FALSE;
      }
    }
    else {
      state->lgar_mass_balance.cache_fluxes = FALSE;
    }
  }

  if (caching_at_start && !state->lgar_mass_balance.cache_fluxes){
    switch_caching = TRUE;//if you switch from cached to not, you need to add the "missing" PET back into the mass balance and AET calculation 
  }

  if (state->lgar_mass_balance.cache_fluxes){
    state->lgar_bmi_params.cache_count ++;
  }

  state->lgar_bmi_params.forcing_interval = int(state->lgar_bmi_params.forcing_resolution_h/state->lgar_bmi_params.timestep_h+1.0e-08); // add 1.0e-08 to prevent truncation error
  subcycles = state->lgar_bmi_params.forcing_interval;

  if (verbosity.compare("high") == 0) {
    printf("time step size in hours: %lf \n", state->lgar_bmi_params.timestep_h);
  }

  // ensure precip and PET are non-negative
  if (state->lgar_bmi_input_params->precipitation_mm_per_h < 0.0) {
      std::cerr << "Warning: Pr [mm/h] (timestep) is negative ("
                << state->lgar_bmi_input_params->precipitation_mm_per_h
                << "), setting to 0.\n";
      state->lgar_bmi_input_params->precipitation_mm_per_h = 0.0;
  }

  if (state->lgar_bmi_input_params->PET_mm_per_h < 0.0) {
      std::cerr << "Warning: PET [mm/h] (timestep) is negative ("
                << state->lgar_bmi_input_params->PET_mm_per_h
                << "), setting to 0.\n";
      state->lgar_bmi_input_params->PET_mm_per_h = 0.0;
  }

  if (PET_affects_precip){ // if the user wants PET subtracted from precip
    if (state->lgar_bmi_input_params->precipitation_mm_per_h > state->lgar_bmi_input_params->PET_mm_per_h){
      state->lgar_bmi_input_params->precipitation_mm_per_h = state->lgar_bmi_input_params->precipitation_mm_per_h - state->lgar_bmi_input_params->PET_mm_per_h;
      state->lgar_bmi_input_params->PET_mm_per_h = 0.0;
    }
    else{
      state->lgar_bmi_input_params->PET_mm_per_h = state->lgar_bmi_input_params->PET_mm_per_h - state->lgar_bmi_input_params->precipitation_mm_per_h;
      state->lgar_bmi_input_params->precipitation_mm_per_h = 0.0;
    }
  }

  if ( (verbosity.compare("high") == 0) && (PET_affects_precip)) {
    std::cerr<<"Pr  [cm/h] (timestep), after PET is subtracted from precip = "<<state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm <<"\n";
    std::cerr<<"PET [cm/h] (timestep), after PET is subtracted from precip = "<<state->lgar_bmi_input_params->PET_mm_per_h * mm_to_cm <<"\n"; 
  }
  
  // subcycling loop (loop over model's timestep)
  for (int cycle=1; cycle <= subcycles; cycle++) {

    bool top_near_sat = false;
    this->state->lgar_bmi_params.time_s    += subtimestep_h * state->units.hr_to_sec;
    this->state->lgar_bmi_params.timesteps ++;

    double precip_for_CR_subtimestep_cm_per_h = 0.0;
    precip_subtimestep_cm_per_h = state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm; // rate [cm/hour]
    
    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      std::cerr<<"BMI Update |---------------------------------------------------------------|\n";
      std::cerr<<"BMI Update |Timesteps = "<< state->lgar_bmi_params.timesteps<<", Time [h] = "<<this->state->lgar_bmi_params.time_s / 3600.<<", Subcycle = "<< cycle <<" of "<<subcycles<<std::endl;
    }

    if( state->state_previous != NULL ){
      listDelete(state->state_previous);
      state->state_previous = NULL;
    }
    state->state_previous = listCopy(state->head);

    double ponded_flux_for_CR = 0.0;

    // allocates some water to conceptual reservoir storage via conditional preferential flow
    if (state->lgar_bmi_params.runoff_in_prev_step){
      double precip_subtimestep_cm_per_h_total = precip_subtimestep_cm_per_h;
      precip_for_CR_subtimestep_cm_per_h = frac_to_CR * precip_subtimestep_cm_per_h_total;
      precip_subtimestep_cm_per_h = (1.0 - frac_to_CR) * precip_subtimestep_cm_per_h_total;

      ponded_flux_for_CR = volon_timestep_cm/subtimestep_h*frac_to_CR;
      volon_timestep_cm = volon_timestep_cm*(1 - frac_to_CR);
    }

    /* Note unit conversion:
       Pr and PET are rates (fluxes) in mm/h
       Pr [mm/h] * 1h/3600sec = Pr [mm/3600sec]
       Model timestep (dt) = 300 sec (5 minutes for example)
       convert rate to amount
       Pr [mm/3600sec] * dt [300 sec] = Pr[mm] * 300/3600.
       in the code below, subtimestep_h is this 300/3600 factor (see initialize from config in lgar.cxx)
    */

    AET_subtimestep_cm            = 0.0;
    volstart_subtimestep_cm       = 0.0;
    volin_subtimestep_cm          = 0.0;
    volrunoff_subtimestep_cm      = 0.0;
    volrech_subtimestep_cm        = 0.0;
    double temp_rch               = 0.0; //handles case when a fraction of a wetting front technically crosses the lower boundary of the vadose zone
    double free_drainage_subtimestep_cm = 0.0;
    double free_drainage_for_CR = 0.0;

    PET_subtimestep_cm_per_h = state->lgar_bmi_input_params->PET_mm_per_h * mm_to_cm;

    ponded_depth_subtimestep_cm = precip_subtimestep_cm_per_h * subtimestep_h; // the amount of water on the surface before any infiltration and runoff

    ponded_depth_subtimestep_cm += volon_timestep_cm; // add volume of water on the surface (from the last timestep) to ponded depth as well

    precip_subtimestep_cm = precip_subtimestep_cm_per_h * subtimestep_h; // rate x dt = amount (portion of the water on the suface for model's timestep [cm])
    PET_subtimestep_cm = PET_subtimestep_cm_per_h * subtimestep_h;      // potential ET for this subtimestep [cm]

    volstart_subtimestep_cm = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);

    if (!state->lgar_mass_balance.cache_fluxes){

      //this code makes sure that AET or free drainage will not be extracted in a way that would result in an impossible storage
      double min_storage = 0.0;
      double mass_used_to_check_impossible_storages = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);
      for (int k = 1; k < num_layers+1; k++) {
        int layer_num_min_check = k;
        int soil_num_min_check = state->lgar_bmi_params.layer_soil_type[layer_num_min_check];
        min_storage += state->soil_properties[soil_num_min_check].theta_r * (state->lgar_bmi_params.cum_layer_thickness_cm[k]-state->lgar_bmi_params.cum_layer_thickness_cm[k-1]);
      }

      int wf_free_drainage_demand = wetting_front_free_drainage(state->head);

      double min_water_possible_for_FD_WF = calc_min_water_possible_for_free_drainage_wetting_front(wf_free_drainage_demand,  &state->head, state->lgar_bmi_params.layer_soil_type, state->soil_properties);
      double storage_in_FD_WF = calc_storage_in_free_drainage_wetting_front(wf_free_drainage_demand, &state->head);

      double mass_correction_for_cached_free_drainage_fluxes = 0.0;

      if (switch_caching){
        PET_subtimestep_cm_per_h += state->lgar_mass_balance.accumulated_PET;
        state->lgar_mass_balance.accumulated_PET = 0.0;
        mass_correction_for_cached_free_drainage_fluxes += state->lgar_mass_balance.accumulated_free_drainage;
        state->lgar_mass_balance.accumulated_free_drainage = 0.0;
      }

      if (state->lgar_bmi_params.free_drainage_enabled){
        struct wetting_front *front = listFindFront(listLength(state->head), state->head, NULL);
        free_drainage_subtimestep_cm += subtimestep_h*front->K_cm_per_h;
        int iter_mass_check_FD = 0;
        while ( (mass_used_to_check_impossible_storages - free_drainage_subtimestep_cm < min_storage) || (storage_in_FD_WF - free_drainage_subtimestep_cm < min_water_possible_for_FD_WF) ){
          free_drainage_subtimestep_cm *= 0.5; //give it a chance to merely become smaller before setting to 0
          if (iter_mass_check_FD > 5){
            free_drainage_subtimestep_cm = 0.0;
            break;
          }
          iter_mass_check_FD ++;
        }
        if (free_drainage_subtimestep_cm<1.E-7){
          free_drainage_subtimestep_cm = 0.0;
        }
        if (front->psi_cm>1.E6){
          free_drainage_subtimestep_cm = 0.0;
        }
      }

      if ((state->lgar_bmi_params.free_drainage_enabled) && verbosity.compare("high") == 0){
        printf("free_drainage_subtimestep_cm: %.10lf \n", free_drainage_subtimestep_cm);
      }

      //using cerr instead of cout due to some cout buffering issues when running in the ngen framework, cerr doesn't buffer so it prints immediately to the sreeen.
      if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {

        std::cerr<<"Pr [cm/h], Pr [cm] (subtimestep), subtimestep [h] = "<<state->lgar_bmi_input_params->precipitation_mm_per_h * mm_to_cm <<", "<< precip_subtimestep_cm <<", "<< subtimestep_h<<" ("<<subtimestep_h*3600<<" sec)"<<"\n";
        std::cerr<<"PET [cm/h], PET [cm] (subtimestep) = "<<state->lgar_bmi_input_params->PET_mm_per_h * mm_to_cm <<", "<< PET_subtimestep_cm<<"\n";
      }

      precip_previous_subtimestep_cm = state->lgar_bmi_params.precip_previous_timestep_cm; // creation of a new wetting front depends on previous timestep's rainfall

      num_layers = state->lgar_bmi_params.num_layers;
      double delta_theta;   // the width of a front, such that its volume=depth*delta_theta
      double dry_depth;

      // Calculate AET from PET if PET is non-zero
      if (PET_subtimestep_cm_per_h > 0.0) {
        AET_subtimestep_cm = calc_aet(PET_subtimestep_cm_per_h, subtimestep_h, wilting_point_psi_cm, field_capacity_psi_cm,
                                      state->lgar_bmi_params.layer_soil_type, AET_thresh_Theta, AET_expon,
                                      state->head, state->soil_properties);
      }

      int iter_mass_check_AET = 0;
      while ( (mass_used_to_check_impossible_storages - AET_subtimestep_cm < min_storage) || (storage_in_FD_WF - AET_subtimestep_cm < min_water_possible_for_FD_WF) ){
        AET_subtimestep_cm *= 0.5; //give it a chance to merely become smaller before setting to 0
        if (iter_mass_check_AET > 5){
          AET_subtimestep_cm = 0.0;
          break;
        }
        iter_mass_check_AET ++;
      }

      int iter_mass_check_AET_and_FD = 0;
      while ( ( (mass_used_to_check_impossible_storages - AET_subtimestep_cm - free_drainage_subtimestep_cm - mass_correction_for_cached_free_drainage_fluxes) < min_storage) || ( (storage_in_FD_WF - AET_subtimestep_cm - free_drainage_subtimestep_cm - mass_correction_for_cached_free_drainage_fluxes) < min_water_possible_for_FD_WF) ){ 
        // both should also be checked at the same because while individually these might not make an impossible storage, together they might
        AET_subtimestep_cm *= 0.5;
        free_drainage_subtimestep_cm *= 0.5;
        mass_correction_for_cached_free_drainage_fluxes *=0.5;
        if (iter_mass_check_AET_and_FD > 5){
          AET_subtimestep_cm = 0.0;
          free_drainage_subtimestep_cm = 0.0;
          mass_correction_for_cached_free_drainage_fluxes = 0.0;
          break;
        }
        iter_mass_check_AET_and_FD ++;
      }

      // precip_timestep_cm += precip_subtimestep_cm;
      precip_timestep_cm += precip_subtimestep_cm + precip_for_CR_subtimestep_cm_per_h*subtimestep_h;
      PET_timestep_cm += fmax(PET_subtimestep_cm,0.0); // ensures non-negative PET

      //addressed machine precision issues where volon_timestep_error could be for example -1E-17 or 1.E-20 or smaller
      volon_timestep_cm = fmax(volon_timestep_cm,0.0);
      volon_timestep_cm = volon_timestep_cm > SMALL_EPS ? volon_timestep_cm : 0.0;

      /*----------------------------------------------------------------------*/
      // Should a new wetting front be created?
      int soil_num = state->lgar_bmi_params.layer_soil_type[state->head->layer_num];
      double theta_e = state->soil_properties[soil_num].theta_e;
      bool is_top_wf_saturated = (state->head->theta+SMALL_EPS) >= theta_e ? true : false; //sometimes a machine precision error would erroneously create a new wetting front during saturated conditions. The epsillon seems to prevent this.
      double theta_above_which_precip_contribs_to_GW = theta_e * spf_factor;
      top_near_sat = state->head->theta > theta_above_which_precip_contribs_to_GW ? true : false; //is the top WF near saturation, thus triggering simple preferential flow if enabled

      // checks on creatign a new surficial front
      // 1. check current and previous timestep precipitation
      // bool create_surficial_front = (precip_previous_subtimestep_cm == 0.0 && precip_subtimestep_cm > 0.0);
      bool create_surficial_front = (precip_previous_subtimestep_cm == 0.0 && precip_subtimestep_cm > 0.0 && volon_timestep_cm == 0) || ( (precip_subtimestep_cm > 0.0 || volon_timestep_cm > 0) && (listLength(state->head)==num_layers) );
      
      // 2. check soil top wetting front condition (saturated/unsaturated), and surface ponded water
      if (is_top_wf_saturated || volon_timestep_cm > 0.0)
        create_surficial_front = false;

      if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
        std::string flag        = (create_surficial_front && !is_top_wf_saturated) == true ? "Yes" : "No";
        std::string flag_top_wf = is_top_wf_saturated == true ? "Yes" : "No";
        std::cerr<<"Is top wetting front saturated? "<< flag_top_wf  << "\n";
        std::cerr<<"Create superficial wetting front? "<< flag << "\n";
      }

      /*----------------------------------------------------------------------*/
      /* create a new wetting front if the following is true. Meaning there is no
        wetting front in the top layer to accept the water, must create one. */
      if(create_surficial_front) {

        double temp_pd = 0.0; // necessary to assign zero precip due to the creation of new wetting front; AET will still be taken out of the layers

        // move the wetting fronts without adding any water; this is done to close the mass balance
        // and also to merge / cross if necessary 
        temp_rch = lgar_move_wetting_fronts(subtimestep_h, &free_drainage_subtimestep_cm, &temp_pd, wf_free_drainage_demand, volend_subtimestep_cm, mass_correction_for_cached_free_drainage_fluxes,
              num_layers, &AET_subtimestep_cm, state->lgar_bmi_params.cum_layer_thickness_cm,
              state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.frozen_factor,
              &state->head, state->state_previous, state->soil_properties);

        // if (temp_pd != 0.0){ //if temp_pd != 0.0, that means that some water left the model through the lower model bdy. For LGARTO preparation, this has been refactored such that temp_rch handles this now.
        //   // volrech_subtimestep_cm = temp_pd;
        //   // volrech_timestep_cm += volrech_subtimestep_cm;
        //   // temp_pd = 0.0;
        // }
        
        // depth of the surficial front to be created
        dry_depth = lgar_calc_dry_depth(use_closed_form_G, nint, subtimestep_h, &delta_theta, state->lgar_bmi_params.layer_soil_type,
                state->lgar_bmi_params.cum_layer_thickness_cm, state->lgar_bmi_params.frozen_factor,
                state->head, state->soil_properties);

        if (verbosity.compare("high") == 0) {
          printf("State before moving creating new WF...\n");
          listPrint(state->head);
        }
        
        lgar_create_surficial_front(num_layers, &ponded_depth_subtimestep_cm, &volin_subtimestep_cm, dry_depth, state->head->theta,
            state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.cum_layer_thickness_cm,
            state->lgar_bmi_params.frozen_factor, &state->head, state->soil_properties);

        if (verbosity.compare("high") == 0) {
          printf("State after moving creating new WF...\n");
          listPrint(state->head);
        }

        if(state->state_previous != NULL ){
          listDelete(state->state_previous);
          state->state_previous = NULL;
        }
        state->state_previous = listCopy(state->head);

        // volin_timestep_cm += volin_subtimestep_cm;

        if (verbosity.compare("high") == 0) {
    std::cerr<<"New wetting front created...\n";
    listPrint(state->head);
        }
      }

      /*----------------------------------------------------------------------*/
      /* infiltrate water based on the infiltration capacity given no new wetting front
        is created and that there is water on the surface (or raining). */

      if (ponded_depth_subtimestep_cm > 0 && !create_surficial_front) {
        volrunoff_subtimestep_cm = lgar_insert_water(use_closed_form_G, nint, subtimestep_h, AET_subtimestep_cm, free_drainage_subtimestep_cm, &ponded_depth_subtimestep_cm,
                &volin_subtimestep_cm, precip_subtimestep_cm_per_h,
                wf_free_drainage_demand, num_layers,
                ponded_depth_max_cm, state->lgar_bmi_params.layer_soil_type,
                state->lgar_bmi_params.cum_layer_thickness_cm,
                state->lgar_bmi_params.frozen_factor, state->head,
                state->soil_properties); 
        // volin_timestep_cm += volin_subtimestep_cm;
        // volrunoff_timestep_cm += volrunoff_subtimestep_cm;
        
        volrech_subtimestep_cm = volin_subtimestep_cm;
        volon_subtimestep_cm = ponded_depth_subtimestep_cm;
        if (volrunoff_subtimestep_cm < 0) {
          std::cerr<<"Runoff is less than 0, which should not happen.\n";
          abort();
        }
      }
      else {

        if (ponded_depth_subtimestep_cm < ponded_depth_max_cm) {
          volrunoff_timestep_cm += 0.0;
          volon_subtimestep_cm = ponded_depth_subtimestep_cm;
          ponded_depth_subtimestep_cm = 0.0;
          volrunoff_subtimestep_cm = 0.0;
        }
        else {
          volrunoff_subtimestep_cm = (ponded_depth_subtimestep_cm - ponded_depth_max_cm);
          volon_subtimestep_cm = ponded_depth_max_cm;
          ponded_depth_subtimestep_cm = ponded_depth_max_cm;
        }
      }
      /*----------------------------------------------------------------------*/

      /* move wetting fronts if no new wetting front is created. Otherwise, movement
        of wetting fronts has already happened at the time of creating surficial front,
        so no need to move them here. */
      if (!create_surficial_front) {
        double volin_subtimestep_cm_temp = volin_subtimestep_cm;  /* passing this for mass balance only, the method modifies it
                    and returns percolated value, so we need to keep its original
                    value stored to copy it back*/
        temp_rch = lgar_move_wetting_fronts(subtimestep_h, &free_drainage_subtimestep_cm, &volin_subtimestep_cm, wf_free_drainage_demand, volend_subtimestep_cm, mass_correction_for_cached_free_drainage_fluxes,
              num_layers, &AET_subtimestep_cm, state->lgar_bmi_params.cum_layer_thickness_cm,
              state->lgar_bmi_params.layer_soil_type, state->lgar_bmi_params.frozen_factor,
              &state->head, state->state_previous, state->soil_properties);

        // this is the volume of water leaving through the bottom
        volrech_subtimestep_cm = volin_subtimestep_cm;
        // volrech_timestep_cm += volrech_subtimestep_cm;

        volin_subtimestep_cm = volin_subtimestep_cm_temp;
      }

      lgar_clean_redundant_fronts(&state->head, state->lgar_bmi_params.layer_soil_type, state->soil_properties); //deletes WFs that are very close in capillary head value

      /*----------------------------------------------------------------------*/
      // calculate derivative (dz/dt) for all wetting fronts
      int new_front = -1;
      if (switch_caching){
        if (create_surficial_front){
          new_front = state->head->front_num;
        }
      }
      lgar_dzdt_calc(use_closed_form_G, nint, num_layers, ponded_depth_subtimestep_cm, subtimestep_h, state->lgar_bmi_params.layer_soil_type,
        state->lgar_bmi_params.cum_layer_thickness_cm, state->lgar_bmi_params.frozen_factor,
        state->head, state->soil_properties, switch_caching, state->lgar_bmi_params.cache_count, new_front);

      if (switch_caching){
        state->lgar_bmi_params.cache_count = 1;
      }

      volrech_subtimestep_cm = temp_rch;
      if (!state->lgar_bmi_params.free_drainage_to_CR){
        volrech_subtimestep_cm += free_drainage_subtimestep_cm;
      }
      else {
        free_drainage_for_CR += free_drainage_subtimestep_cm + volrech_subtimestep_cm;
        volrech_subtimestep_cm = 0.0;
      }

    }

    else {//in this case, we just use cached fluxes in order to save time. No direct computation of fluxes are necessary, and the ones from the last time step are used.
      //also in this case, there will be no runoff due to precipitation partitioning because there is no precipitation 
      PET_timestep_cm += fmax(PET_subtimestep_cm,0.0); // ensures non-negative PET
      state->lgar_mass_balance.accumulated_PET += PET_timestep_cm;

      if (!state->lgar_bmi_params.free_drainage_to_CR){
        volrech_subtimestep_cm += state->lgar_mass_balance.previous_recharge;
      }
      else {
        free_drainage_for_CR += state->lgar_mass_balance.previous_recharge;
      }
      state->lgar_mass_balance.accumulated_free_drainage += state->lgar_mass_balance.previous_recharge;
    }

    volend_subtimestep_cm = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);
    volend_timestep_cm = volend_subtimestep_cm;
    state->lgar_bmi_params.precip_previous_timestep_cm = precip_subtimestep_cm;

    volCRstart_subtimestep_cm = state->lgar_mass_balance.CR_fast_storage_cm + state->lgar_mass_balance.CR_slow_storage_cm;
    double volin_CR_subtimestep_cm = (precip_for_CR_subtimestep_cm_per_h + ponded_flux_for_CR)*subtimestep_h + free_drainage_for_CR;
    double volQ_CR_subtimestep_cm = calc_CR_Q(subtimestep_h, a, a_slow, b, b_slow, frac_slow, precip_for_CR_subtimestep_cm_per_h + ponded_flux_for_CR + free_drainage_for_CR/subtimestep_h, &state->lgar_mass_balance.CR_fast_storage_cm, &state->lgar_mass_balance.CR_slow_storage_cm);
    state->lgar_mass_balance.volrunoff_CR_cm += volQ_CR_subtimestep_cm;
    volQ_CR_timestep_cm += volQ_CR_subtimestep_cm;
    volCRend_subtimestep_cm = state->lgar_mass_balance.CR_fast_storage_cm + state->lgar_mass_balance.CR_slow_storage_cm;
    volCRend_timestep_cm = volCRend_subtimestep_cm;

    // set runoff_in_prev_step for next step
    if ((volrunoff_subtimestep_cm > SMALL_EPS) || (top_near_sat)){
      state->lgar_bmi_params.runoff_in_prev_step = true;
    }
    else {
      state->lgar_bmi_params.runoff_in_prev_step = false;
    }

    //add precip_for_CR_subtimestep_cm_per_h back into precip for mass balance, note that this means volin_CR_timestep is not necessary for mass balance 
    precip_subtimestep_cm += precip_for_CR_subtimestep_cm_per_h * subtimestep_h;

    /*----------------------------------------------------------------------*/
    // mass balance at the subtimestep (local mass balance)

    double local_mb = volstart_subtimestep_cm + precip_subtimestep_cm + volon_timestep_cm - volrunoff_subtimestep_cm - volQ_CR_subtimestep_cm - volCRend_subtimestep_cm + volCRstart_subtimestep_cm 
                      - AET_subtimestep_cm - volon_subtimestep_cm - volrech_subtimestep_cm - volend_subtimestep_cm;


    /*----------------------------------------------------------------------*/

    ///////
    //separating code such that most non substep vars (so xxx_timestep and not xxx_subtimestep) are updated in just one place. not all, because some must be set before substepping.
    volin_timestep_cm += volin_subtimestep_cm;
    volrech_timestep_cm += volrech_subtimestep_cm;

    volrunoff_timestep_cm += volrunoff_subtimestep_cm;

    AET_timestep_cm += AET_subtimestep_cm;
    volon_timestep_cm = volon_subtimestep_cm; // surface ponded water at the end of the timestep 
    ///////
    
    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      printf("Printing wetting fronts at this subtimestep... \n");
      listPrint(state->head);
    }

    bool unexpected_local_error = fabs(local_mb) > mbal_tol ? true : false; //1.0E-4 was the default for initial stability testing 
    if (isinf(local_mb)){
      unexpected_local_error = true;
    }
    
    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0 || unexpected_local_error) {
      if (!state->lgar_bmi_params.frac_to_CR){
        printf("\nLocal mass balance at this timestep... \n\
        Error         = %14.10f \n\
        Initial water = %14.10f \n\
        Water added   = %14.10f \n\
        Ponded water  = %14.10f \n\
        Infiltration  = %14.10f \n\
        Runoff        = %14.10f \n\
        AET           = %14.10f \n\
        Percolation   = %14.10f \n\
        Final water   = %14.10f \n", local_mb, volstart_subtimestep_cm, precip_subtimestep_cm, volon_subtimestep_cm,
        volin_subtimestep_cm, volrunoff_subtimestep_cm, AET_subtimestep_cm, volrech_subtimestep_cm,
        volend_subtimestep_cm);
      }
      else {
        printf("\nLocal mass balance at this timestep... \n\
        Error                   = %14.10f \n\
        Initial water (LGAR)    = %14.10f \n\
        Water added (total)     = %14.10f \n\
        Ponded water            = %14.10f \n\
        Infiltration (LGAR)     = %14.10f \n\
        Runoff (total)          = %14.10f \n\
        AET                     = %14.10f \n\
        Percolation             = %14.10f \n\
        Final water (LGAR)      = %14.10f \n\
        Water added (con res)   = %14.10f \n\
        Initial water (con res) = %14.10f \n\
        Final water (con res)   = %14.10f \n\
        Runoff (con res)        = %14.10f \n", local_mb, volstart_subtimestep_cm, precip_subtimestep_cm, volon_subtimestep_cm,
        volin_subtimestep_cm, volrunoff_subtimestep_cm, AET_subtimestep_cm, volrech_subtimestep_cm,
        volend_subtimestep_cm, volin_CR_subtimestep_cm, volCRstart_subtimestep_cm,
        volCRend_subtimestep_cm, volQ_CR_subtimestep_cm);
      }

      if (unexpected_local_error) {
	printf("Local mass balance (in this timestep) is %14.10f, larger than expected, needs some debugging...\n ",local_mb);
	abort();
      }

    }

    // store local mass balance error to the struct
    state->lgar_mass_balance.local_mass_balance = local_mb;

    assert (state->head->depth_cm > 0.0); // check on negative layer depth --> move this to somewhere else AJ (later)

    bool lasam_standalone = true;
#ifdef NGEN
    lasam_standalone = false;
#endif
    // simuation time can't exceed the endtime when running standalone
    if ( (this->state->lgar_bmi_params.time_s >= this->state->lgar_bmi_params.endtime_s) && lasam_standalone)
      break;

  } // end of subcycling

  //update giuh at the time step level (was previously updated at the sub time step level)
  volrunoff_giuh_timestep_cm = giuh_convolution_integral(volrunoff_timestep_cm + volQ_CR_timestep_cm , num_giuh_ordinates, giuh_ordinates, giuh_runoff_queue);

  // total mass of water leaving the system, at this time it is the giuh-only, but later will add groundwater component as well.
  // when groundwater component is added, it should probably happen inside of the subcycling loop.
  volQ_timestep_cm = volrunoff_giuh_timestep_cm; //note that volQ_CR_timestep_cm was added to volQ_timestep_cm in the input for giuh_convolution_integral

  /*----------------------------------------------------------------------*/
  // Everything related to lgar state is done at this point, now time to update some dynamic variables

  // update number of wetting fronts
  state->lgar_bmi_params.num_wetting_fronts = listLength(state->head);

  // allocate new memory based on updated wetting fronts; we could make it conditional i.e. create only if no. of wf are changed
  realloc_soil();

  // update thickness/depth and soil moisture of wetting fronts (used for state coupling)
  struct wetting_front *current = state->head;
  int to_bottom_count = 0;
  for (int i=0; i<state->lgar_bmi_params.num_wetting_fronts; i++) {
    if (current->to_bottom){
      to_bottom_count ++;
    }
    if (to_bottom_count>num_layers){
      std::cerr << "Error: too many to_bottom WFs! This should be equal to the number of layers.\n";
      listPrint(state->head);
	    abort();
    }
    assert (current != NULL);
    state->lgar_bmi_params.soil_moisture_wetting_fronts[i] = current->theta;
    state->lgar_bmi_params.soil_depth_wetting_fronts[i] = current->depth_cm * state->units.cm_to_m;
    current = current->next;
    if (verbosity.compare("high") == 0)
      std::cerr<<"Wetting fronts (bmi outputs) (depth in meters, theta)= "
	       <<state->lgar_bmi_params.soil_depth_wetting_fronts[i]
	       <<" "<<state->lgar_bmi_params.soil_moisture_wetting_fronts[i]<<"\n";
  }
  
  // add to mass balance timestep variables
  state->lgar_mass_balance.volprecip_timestep_cm  = precip_timestep_cm;
  state->lgar_mass_balance.volin_timestep_cm      = volin_timestep_cm;
  state->lgar_mass_balance.volon_timestep_cm      = volon_timestep_cm;
  state->lgar_mass_balance.volend_timestep_cm     = volend_timestep_cm;
  state->lgar_mass_balance.volCRend_timestep_cm   = volCRend_timestep_cm;
  state->lgar_mass_balance.volAET_timestep_cm     = AET_timestep_cm;
  state->lgar_mass_balance.volrech_timestep_cm    = volrech_timestep_cm;
  state->lgar_mass_balance.volrunoff_timestep_cm  = volrunoff_timestep_cm;
  state->lgar_mass_balance.volQ_timestep_cm       = volQ_timestep_cm;
  state->lgar_mass_balance.volQ_CR_timestep_cm    = volQ_CR_timestep_cm;
  state->lgar_mass_balance.volPET_timestep_cm     = PET_timestep_cm;
  state->lgar_mass_balance.volrunoff_giuh_timestep_cm = volrunoff_giuh_timestep_cm;

  //for caching
  state->lgar_mass_balance.previous_AET         = AET_subtimestep_cm;
  state->lgar_mass_balance.previous_PET         = PET_subtimestep_cm;
  if (!state->lgar_mass_balance.cache_fluxes){
    state->lgar_mass_balance.previous_recharge    = volrech_subtimestep_cm;
  }

  // add to mass balance accumulated variables
  state->lgar_mass_balance.volprecip_cm  += precip_timestep_cm;
  state->lgar_mass_balance.volin_cm      += volin_timestep_cm;
  state->lgar_mass_balance.volon_cm       = volon_timestep_cm;
  state->lgar_mass_balance.volend_cm      = volend_timestep_cm;
  state->lgar_mass_balance.volCRend_cm    = volCRend_timestep_cm;
  state->lgar_mass_balance.volAET_cm     += AET_timestep_cm;
  state->lgar_mass_balance.volrech_cm    += volrech_timestep_cm;
  state->lgar_mass_balance.volrunoff_cm  += volrunoff_timestep_cm;
  state->lgar_mass_balance.volQ_cm       += volQ_timestep_cm;
  state->lgar_mass_balance.volQ_CR_cm    += volQ_CR_timestep_cm;
  state->lgar_mass_balance.volPET_cm     += PET_timestep_cm;
  state->lgar_mass_balance.volrunoff_giuh_cm  += volrunoff_giuh_timestep_cm;
  state->lgar_mass_balance.volchange_calib_cm += volchange_calib_cm ;
 
  // converted values, a struct local to the BMI and has bmi output variables
  bmi_unit_conv.mass_balance_m        = state->lgar_mass_balance.local_mass_balance * state->units.cm_to_m;
  bmi_unit_conv.volprecip_timestep_m  = precip_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volin_timestep_m      = volin_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volend_timestep_m     = volend_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volCRend_timestep_m   = volCRend_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volAET_timestep_m     = AET_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volrech_timestep_m    = volrech_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volrunoff_timestep_m  = volrunoff_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volQ_timestep_m       = volQ_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volQ_CR_timestep_m    = volQ_CR_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volPET_timestep_m     = PET_timestep_cm * state->units.cm_to_m;
  bmi_unit_conv.volrunoff_giuh_timestep_m = volrunoff_giuh_timestep_cm * state->units.cm_to_m;
  
}


void BmiLGAR::
UpdateUntil(double t)
{
  assert (t > 0.0);
  this->Update();
}

struct model_state* BmiLGAR::get_model()
{
  return state;
}

void BmiLGAR::
global_mass_balance()
{
  lgar_global_mass_balance(this->state, giuh_runoff_queue);
}

double BmiLGAR::
update_calibratable_parameters()
{
  int soil, layer_num;
  struct wetting_front *current = state->head;

  if (verbosity.compare("high") == 0)
    listPrint(state->head);
  
  double volstart_before = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);

  // first we update the parameters that depend on soil layer, for each layer. 
  // This no longer relies on arrays (for compatability with ngen-cal) and instead supports calibration of up to 3 soil layers with scalar parameters.
  for (int i=0; i<state->lgar_bmi_params.num_wetting_fronts; i++) {
    layer_num  = current->layer_num;
    soil = state->lgar_bmi_params.layer_soil_type[layer_num];
    
    assert (current != NULL);

    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      std::cerr<<"----------- Calibratable parameters depending on soil layer (initial values) ----------- \n";
      std::cerr<<"| soil_type = "<< soil <<", layer = "<<layer_num
	       <<", smcmax = "   << state->soil_properties[soil].theta_e
	       <<", smcmin = "   << state->soil_properties[soil].theta_r
	       <<", vg_n = "     << state->soil_properties[soil].vg_n
	       <<", vg_alpha = " << state->soil_properties[soil].vg_alpha_per_cm
	       <<", Ksat = "     << state->soil_properties[soil].Ksat_cm_per_h
	       <<", theta = "    << current->theta <<"\n";
    }
    
    // state->soil_properties[soil].theta_e = state->lgar_calib_params.theta_e[layer_num-1];
    // state->soil_properties[soil].theta_r = state->lgar_calib_params.theta_r[layer_num-1];
    // state->soil_properties[soil].vg_n    = state->lgar_calib_params.vg_n[layer_num-1];
    // state->soil_properties[soil].vg_m    = 1.0 - 1.0/state->soil_properties[soil].vg_n;
    // state->soil_properties[soil].vg_alpha_per_cm = state->lgar_calib_params.vg_alpha[layer_num-1];
    // state->soil_properties[soil].Ksat_cm_per_h   = state->lgar_calib_params.Ksat[layer_num-1];

    if (layer_num==1){
      state->soil_properties[soil].theta_e = state->lgar_calib_params.theta_e_1;
      state->soil_properties[soil].theta_r = state->lgar_calib_params.theta_r_1;
      state->soil_properties[soil].vg_n    = state->lgar_calib_params.vg_n_1;
      state->soil_properties[soil].vg_m    = 1.0 - 1.0/state->soil_properties[soil].vg_n;
      state->soil_properties[soil].vg_alpha_per_cm = state->lgar_calib_params.vg_alpha_1;
      state->soil_properties[soil].Ksat_cm_per_h   = state->lgar_calib_params.Ksat_1;
      if (state->lgar_bmi_params.log_mode){
        state->soil_properties[soil].vg_alpha_per_cm = pow(10.0, state->lgar_calib_params.vg_alpha_1);
        state->soil_properties[soil].Ksat_cm_per_h   = pow(10.0, state->lgar_calib_params.Ksat_1);
      }
    }

    if (layer_num==2){
      state->soil_properties[soil].theta_e = state->lgar_calib_params.theta_e_2;
      state->soil_properties[soil].theta_r = state->lgar_calib_params.theta_r_2;
      state->soil_properties[soil].vg_n    = state->lgar_calib_params.vg_n_2;
      state->soil_properties[soil].vg_m    = 1.0 - 1.0/state->soil_properties[soil].vg_n;
      state->soil_properties[soil].vg_alpha_per_cm = state->lgar_calib_params.vg_alpha_2;
      state->soil_properties[soil].Ksat_cm_per_h   = state->lgar_calib_params.Ksat_2;
      if (state->lgar_bmi_params.log_mode){
        state->soil_properties[soil].vg_alpha_per_cm = pow(10.0, state->lgar_calib_params.vg_alpha_2);
        state->soil_properties[soil].Ksat_cm_per_h   = pow(10.0, state->lgar_calib_params.Ksat_2);
      }
    }

    if (layer_num==3){
      state->soil_properties[soil].theta_e = state->lgar_calib_params.theta_e_3;
      state->soil_properties[soil].theta_r = state->lgar_calib_params.theta_r_3;
      state->soil_properties[soil].vg_n    = state->lgar_calib_params.vg_n_3;
      state->soil_properties[soil].vg_m    = 1.0 - 1.0/state->soil_properties[soil].vg_n;
      state->soil_properties[soil].vg_alpha_per_cm = state->lgar_calib_params.vg_alpha_3;
      state->soil_properties[soil].Ksat_cm_per_h   = state->lgar_calib_params.Ksat_3;
      if (state->lgar_bmi_params.log_mode){
        state->soil_properties[soil].vg_alpha_per_cm = pow(10.0, state->lgar_calib_params.vg_alpha_3);
        state->soil_properties[soil].Ksat_cm_per_h   = pow(10.0, state->lgar_calib_params.Ksat_3);
      }
    }
    
    current->theta = calc_theta_from_h(current->psi_cm, state->soil_properties[soil].vg_alpha_per_cm,
				       state->soil_properties[soil].vg_m, state->soil_properties[soil].vg_n,
				       state->soil_properties[soil].theta_e, state->soil_properties[soil].theta_r);

    if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
      std::cerr<<"----------- Calibratable parameters depending on soil layer (updated values) ----------- \n";
      std::cerr<<"| soil_type = "<< soil <<", layer = "<<layer_num
	       <<", smcmax = "   << state->soil_properties[soil].theta_e
	       <<", smcmin = "   << state->soil_properties[soil].theta_r
	       <<", vg_n = "     << state->soil_properties[soil].vg_n
	       <<", vg_alpha = " << state->soil_properties[soil].vg_alpha_per_cm
	       <<", Ksat = "     << state->soil_properties[soil].Ksat_cm_per_h
	       <<", theta = "    << current->theta <<"\n";
    }
    
    current = current->next;
  }

  //next we update the parameters that apply to the whole model domain and do not depend on soil layer
  if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
    std::cerr<<"----------- Calibratable parameters independent of soil layer (initial values) ----------- \n";
    std::cerr<<"field_capacity_psi = "   << state->lgar_bmi_params.field_capacity_psi_cm
      <<", ponded_depth_max = "     << state->lgar_bmi_params.ponded_depth_max_cm
      <<", a = "     << state->lgar_bmi_params.a
      <<", b = "     << state->lgar_bmi_params.b
      <<", frac_to_CR = "     << state->lgar_bmi_params.frac_to_CR
      <<", spf_factor = "     << state->lgar_bmi_params.spf_factor <<
      "\n";
    if (state->lgar_bmi_params.frac_slow){
      std::cerr
      <<", a_slow = "     << state->lgar_bmi_params.a_slow
      <<", b_slow = "     << state->lgar_bmi_params.b_slow
      <<", frac__slow = "     << state->lgar_bmi_params.frac_slow <<
      "\n";
    }
  }



  state->lgar_bmi_params.field_capacity_psi_cm = state->lgar_calib_params.field_capacity_psi;
  state->lgar_bmi_params.ponded_depth_max_cm   = state->lgar_calib_params.ponded_depth_max;
  state->lgar_bmi_params.a                     = state->lgar_calib_params.a;
  state->lgar_bmi_params.b                     = state->lgar_calib_params.b;
  state->lgar_bmi_params.frac_to_CR            = state->lgar_calib_params.frac_to_CR;
  state->lgar_bmi_params.spf_factor            = state->lgar_calib_params.spf_factor;

  if (state->lgar_bmi_params.frac_slow){
    state->lgar_bmi_params.a_slow              = state->lgar_calib_params.a_slow;
    state->lgar_bmi_params.b_slow              = state->lgar_calib_params.b_slow;
    state->lgar_bmi_params.frac_slow           = state->lgar_calib_params.frac_slow;
  }

  if (state->lgar_bmi_params.log_mode){
    state->lgar_bmi_params.a = pow(10.0, state->lgar_calib_params.a);
    if (state->lgar_bmi_params.frac_slow){
      state->lgar_bmi_params.a_slow = pow(10.0, state->lgar_calib_params.a_slow);
    }
  }

  if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0) {
    std::cerr<<"----------- Calibratable parameters independent of soil layer (updated values) ----------- \n";
    std::cerr<<"field_capacity_psi = "   << state->lgar_bmi_params.field_capacity_psi_cm
      <<", ponded_depth_max = "     << state->lgar_bmi_params.ponded_depth_max_cm
      <<", a = "     << state->lgar_bmi_params.a
      <<", b = "     << state->lgar_bmi_params.b
      <<", frac_to_CR = "     << state->lgar_bmi_params.frac_to_CR
      <<", spf_factor = "     << state->lgar_bmi_params.spf_factor <<
      "\n";
    if (state->lgar_bmi_params.frac_slow){
      std::cerr
      <<", a_slow = "     << state->lgar_bmi_params.a_slow
      <<", b_slow = "     << state->lgar_bmi_params.b_slow
      <<", frac__slow = "     << state->lgar_bmi_params.frac_slow << 
      "\n";
    }
  }
  
  if (verbosity.compare("high") == 0)
    listPrint(state->head);
  
  double volstart_after = lgar_calc_mass_bal(state->lgar_bmi_params.cum_layer_thickness_cm, state->head);

  if (verbosity.compare("high") == 0 || verbosity.compare("low") == 0)
    std::cerr<<"Mass of water (before and after) = "<< volstart_before<<", "<< volstart_after <<"\n";
  
  return volstart_after - volstart_before;
}

void BmiLGAR::
Finalize()
{
  global_mass_balance();
  listDelete(state->head);
  listDelete(state->state_previous);

  delete [] state->soil_properties;

  delete [] state->lgar_bmi_params.soil_depth_wetting_fronts;
  delete [] state->lgar_bmi_params.soil_moisture_wetting_fronts;

  delete [] state->lgar_bmi_params.soil_temperature;
  delete [] state->lgar_bmi_params.soil_temperature_z;
  delete [] state->lgar_bmi_params.layer_soil_type;

  // no longer needed as per layer calibratable parameters are no longer handled with arrays
  // delete [] state->lgar_calib_params.theta_e;
  // delete [] state->lgar_calib_params.theta_r;
  // delete [] state->lgar_calib_params.vg_n;
  // delete [] state->lgar_calib_params.vg_alpha;
  // delete [] state->lgar_calib_params.Ksat;

  delete [] state->lgar_bmi_params.layer_thickness_cm;
  delete [] state->lgar_bmi_params.cum_layer_thickness_cm;
  delete [] state->lgar_bmi_params.giuh_ordinates;
  delete [] state->lgar_bmi_params.frozen_factor;
  delete state->lgar_bmi_input_params;
  delete state;
}


int BmiLGAR::
GetVarGrid(std::string name)
{
  if (name.compare("soil_storage_model") == 0 || name.compare("soil_num_wetting_fronts") == 0)   // int
    return 0;
  else if (name.compare("precipitation_rate") == 0 || name.compare("precipitation") == 0)
    return 1;
  else if (name.compare("potential_evapotranspiration_rate") == 0
	   || name.compare("potential_evapotranspiration") == 0
	   || name.compare("actual_evapotranspiration") == 0) // double
    return 1;
  else if (name.compare("surface_runoff") == 0 || name.compare("giuh_runoff") == 0 || name.compare("a") == 0 || name.compare("b") == 0 || name.compare("frac_to_CR") == 0 || name.compare("spf_factor") == 0
	   || name.compare("a_slow") == 0 || name.compare("b_slow") == 0 || name.compare("frac_slow") == 0 || name.compare("soil_storage") == 0 || name.compare("field_capacity") == 0 || name.compare("ponded_depth_max") == 0)// double
    return 1;
  else if (name.compare("total_discharge") == 0 || name.compare("infiltration") == 0
	   || name.compare("percolation") == 0  || name.compare("conceptual_reservoir_to_stream_discharge") == 0) // double
    return 1;
  else if (name.compare("mass_balance") == 0)
    return 1;
  else if (name.compare("smcmax_1") == 0 || name.compare("smcmin_1") == 0 // per layer parameters are now handled with scalars and not arrays
	   || name.compare("van_genuchten_m_1") == 0 || name.compare("van_genuchten_alpha_1") == 0 || name.compare("van_genuchten_n_1") == 0 
	   || name.compare("hydraulic_conductivity_1") == 0)
    return 1;
  else if (name.compare("smcmax_2") == 0 || name.compare("smcmin_2") == 0 // per layer parameters are now handled with scalars and not arrays
	   || name.compare("van_genuchten_m_2") == 0 || name.compare("van_genuchten_alpha_2") == 0 || name.compare("van_genuchten_n_2") == 0 
	   || name.compare("hydraulic_conductivity_2") == 0)
    return 1;
  else if (name.compare("smcmax_2") == 0 || name.compare("smcmin_2") == 0 // per layer parameters are now handled with scalars and not arrays
	   || name.compare("van_genuchten_m_2") == 0 || name.compare("van_genuchten_alpha_2") == 0 || name.compare("van_genuchten_n_2") == 0 
	   || name.compare("hydraulic_conductivity_2") == 0)
    return 1;
  // else if (name.compare("soil_depth_layers") == 0  || name.compare("smcmax") == 0 || name.compare("smcmin") == 0
	//    || name.compare("van_genuchten_m") == 0 || name.compare("van_genuchten_alpha") == 0 || name.compare("van_genuchten_n") == 0 
	//    || name.compare("hydraulic_conductivity") == 0) // array of doubles (fixed length)
  else if (name.compare("soil_depth_layers") == 0 )// array of doubles (fixed length)
    return 2;
  else if (name.compare("soil_moisture_wetting_fronts") == 0 || name.compare("soil_depth_wetting_fronts") == 0) // array of doubles (dynamic length)
    return 3;
  else if (name.compare("soil_temperature_profile") == 0) // array of doubles (fixed and of the size of soil temperature profile)
    return 4;
  else
    return -1;
}


std::string BmiLGAR::
GetVarType(std::string name)
{
  int var_grid = GetVarGrid(name);

  if (var_grid == 0)
    return "int";
  else if (var_grid == 1 || var_grid == 2 || var_grid == 3 || var_grid == 4)
    return "double";
  else
    return "none";
}


int BmiLGAR::
GetVarItemsize(std::string name)
{
  int var_grid = GetVarGrid(name);

   if (var_grid == 0)
    return sizeof(int);
  else if (var_grid == 1 || var_grid == 2 || var_grid == 3 || var_grid == 4)
    return sizeof(double);
  else
    return 0;
}


std::string BmiLGAR::
GetVarUnits(std::string name)
{
  if (name.compare("precipitation_rate") == 0 || name.compare("potential_evapotranspiration_rate") == 0)
    return "mm h^-1";
  else if (name.compare("precipitation") == 0 || name.compare("potential_evapotranspiration") == 0
	   || name.compare("actual_evapotranspiration") == 0) // double
    return "m";
  else if (name.compare("surface_runoff") == 0 || name.compare("giuh_runoff") == 0 
	   || name.compare("soil_storage") == 0) // double
    return "m";
  else if (name.compare("total_discharge") == 0 || name.compare("infiltration") == 0
	   || name.compare("percolation") == 0) // double
    return "m";
  else if (name.compare("mass_balance") == 0 || name.compare("conceptual_reservoir_to_stream_discharge") == 0)
    return "m";
  else if (name.compare("soil_moisture_wetting_fronts") == 0) // array of doubles
    return "none";
  else if (name.compare("soil_depth_layers") == 0 || name.compare("soil_depth_wetting_fronts") == 0) // array of doubles
    return "m";
  else if (name.compare("soil_temperature_profile") == 0)
    return "K";
  else
    return "none";

}


int BmiLGAR::
GetVarNbytes(std::string name)
{
  int itemsize;
  int gridsize;

  itemsize = this->GetVarItemsize(name);
  gridsize = this->GetGridSize(this->GetVarGrid(name));
  return itemsize * gridsize;
}


std::string BmiLGAR::
GetVarLocation(std::string name)
{
  if (name.compare("precipitation_rate") == 0 || name.compare("precipitation") == 0 ||
      name.compare("potential_evapotranspiration") == 0 || name.compare("potential_evapotranspiration_rate") == 0
      || name.compare("actual_evapotranspiration") == 0) // double
    return "node";
  else if (name.compare("surface_runoff") == 0 || name.compare("giuh_runoff") == 0 || name.compare("a") == 0 || name.compare("b") == 0 || name.compare("frac_to_CR") == 0 || name.compare("spf_factor") == 0
	   || name.compare("a_slow") == 0 || name.compare("b_slow") == 0 || name.compare("frac_slow") == 0 || name.compare("soil_storage") == 0) // double
    return "node";
   else if (name.compare("total_discharge") == 0 || name.compare("infiltration") == 0
	    || name.compare("percolation") == 0 || name.compare("conceptual_reservoir_to_stream_discharge") == 0) // double
    return "node";
  else if (name.compare("soil_moisture_wetting_fronts") == 0) // array of doubles
    return "node";
  else if (name.compare("mass_balance") == 0)
    return "node";
  else if (name.compare("soil_depth_layers") == 0 || name.compare("soil_depth_wetting_fronts") == 0
	   || name.compare("soil_num_wetting_fronts") == 0) // array of doubles
    return "node";
  else if (name.compare("soil_temperature_profile") == 0)
    return "node";
  else
    return "none";
}


void BmiLGAR::
GetGridShape(const int grid, int *shape)
{
  if (grid == 2)
    shape[0] = this->state->lgar_bmi_params.num_layers;
  else if (grid == 3) // number of wetting fronts (dynamic)
    shape[1] = this->state->lgar_bmi_params.num_wetting_fronts;
}


void BmiLGAR::
GetGridSpacing (const int grid, double * spacing)
{
  if (grid == 0) {
    spacing[0] = this->state->lgar_bmi_params.spacing[0];
  }
}


void BmiLGAR::
GetGridOrigin (const int grid, double *origin)
{
  if (grid == 0) {
    origin[0] = this->state->lgar_bmi_params.origin[0];
  }
}


int BmiLGAR::
GetGridRank(const int grid)
{
  if (grid == 0 || grid == 1 || grid == 2 || grid == 3 || grid == 4)
    return 1;
  else
    return -1;
}


int BmiLGAR::
GetGridSize(const int grid)
{
  if (grid == 0 || grid == 1)
    return 1;
  else if (grid == 2) // number of layers (fixed)
    return this->state->lgar_bmi_params.num_layers;
  else if (grid == 3) // number of wetting fronts (dynamic)
    return this->state->lgar_bmi_params.num_wetting_fronts;
  else if (grid == 4) // number of cells (discretized temperature profile, input from SFT)
    return this->state->lgar_bmi_params.num_cells_temp;
  else
    return -1;
}



void BmiLGAR::
GetValue (std::string name, void *dest)
{
  void * src = NULL;
  int nbytes = 0;

  src = this->GetValuePtr(name);
  nbytes = this->GetVarNbytes(name);
  memcpy (dest, src, nbytes);
}


void *BmiLGAR::
GetValuePtr (std::string name)
{
  if (name.compare("precipitation_rate") == 0)
    return (void*)(&this->state->lgar_bmi_input_params->precipitation_mm_per_h);
  else if (name.compare("precipitation") == 0)
    return (void*)(&bmi_unit_conv.volprecip_timestep_m);
  else if (name.compare("potential_evapotranspiration_rate") == 0)
    return (void*)(&this->state->lgar_bmi_input_params->PET_mm_per_h);
  else if (name.compare("potential_evapotranspiration") == 0)
    return (void*)(&bmi_unit_conv.volPET_timestep_m);
  else if (name.compare("actual_evapotranspiration") == 0)
    return (void*)(&bmi_unit_conv.volAET_timestep_m);
  else if (name.compare("surface_runoff") == 0)
    return (void*)(&bmi_unit_conv.volrunoff_timestep_m);
  else if (name.compare("giuh_runoff") == 0)
    return (void*)(&bmi_unit_conv.volrunoff_giuh_timestep_m);
  else if (name.compare("soil_storage") == 0)
    return (void*)(&bmi_unit_conv.volend_timestep_m);
  else if (name.compare("conceptual_reservoir_storage") == 0)
    return (void*)(&bmi_unit_conv.volCRend_timestep_m);
  else if (name.compare("total_discharge") == 0)
    return (void*)(&bmi_unit_conv.volQ_timestep_m);
  else if (name.compare("infiltration") == 0)
    return (void*)(&bmi_unit_conv.volin_timestep_m);
  else if (name.compare("percolation") == 0)
    return (void*)(&bmi_unit_conv.volrech_timestep_m);
  else if (name.compare("conceptual_reservoir_to_stream_discharge") == 0)
    return (void*)(&bmi_unit_conv.volQ_CR_timestep_m);
  else if (name.compare("mass_balance") == 0)
    return (void*)(&bmi_unit_conv.mass_balance_m);
  else if (name.compare("soil_depth_layers") == 0)
    return (void*)this->state->lgar_bmi_params.cum_layer_thickness_cm;  // this too and, if needed, change soil_moisture_layers to soil_thickness_layers
  else if (name.compare("soil_moisture_wetting_fronts") == 0)
    return (void*)this->state->lgar_bmi_params.soil_moisture_wetting_fronts;
  else if (name.compare("soil_depth_wetting_fronts") == 0)
    return (void*)this->state->lgar_bmi_params.soil_depth_wetting_fronts;
  else if (name.compare("soil_num_wetting_fronts") == 0)
    return (void*)(&state->lgar_bmi_params.num_wetting_fronts);
  else if (name.compare("soil_temperature_profile") == 0)
    return (void*)this->state->lgar_bmi_params.soil_temperature;
  // else if (name.compare("smcmax") == 0)
  //   return (void*)this->state->lgar_calib_params.theta_e;
  // else if (name.compare("smcmin") == 0)
  //   return (void*)this->state->lgar_calib_params.theta_r;
  // else if (name.compare("van_genuchten_n") == 0)
  //   return (void*)this->state->lgar_calib_params.vg_n;
  // else if (name.compare("van_genuchten_alpha") == 0)
  //   return (void*)this->state->lgar_calib_params.vg_alpha;
  // else if (name.compare("hydraulic_conductivity") == 0)
  //   return (void*)this->state->lgar_calib_params.Ksat;

  // per layer calibratable params are now scalars and not arrays
  else if (name.compare("smcmax_1") == 0)
    return (void*)&this->state->lgar_calib_params.theta_e_1;
  else if (name.compare("smcmin_1") == 0)
    return (void*)&this->state->lgar_calib_params.theta_r_1;
  else if (name.compare("van_genuchten_n_1") == 0)
    return (void*)&this->state->lgar_calib_params.vg_n_1;
  else if (name.compare("van_genuchten_alpha_1") == 0)
    return (void*)&this->state->lgar_calib_params.vg_alpha_1;
  else if (name.compare("hydraulic_conductivity_1") == 0)
    return (void*)&this->state->lgar_calib_params.Ksat_1;
  else if (name.compare("smcmax_2") == 0)
    return (void*)&this->state->lgar_calib_params.theta_e_2;
  else if (name.compare("smcmin_2") == 0)
    return (void*)&this->state->lgar_calib_params.theta_r_2;
  else if (name.compare("van_genuchten_n_2") == 0)
    return (void*)&this->state->lgar_calib_params.vg_n_2;
  else if (name.compare("van_genuchten_alpha_2") == 0)
    return (void*)&this->state->lgar_calib_params.vg_alpha_2;
  else if (name.compare("hydraulic_conductivity_2") == 0)
    return (void*)&this->state->lgar_calib_params.Ksat_2;
  else if (name.compare("smcmax_3") == 0)
    return (void*)&this->state->lgar_calib_params.theta_e_3;
  else if (name.compare("smcmin_3") == 0)
    return (void*)&this->state->lgar_calib_params.theta_r_3;
  else if (name.compare("van_genuchten_n_3") == 0)
    return (void*)&this->state->lgar_calib_params.vg_n_3;
  else if (name.compare("van_genuchten_alpha_3") == 0)
    return (void*)&this->state->lgar_calib_params.vg_alpha_3;
  else if (name.compare("hydraulic_conductivity_3") == 0)
    return (void*)&this->state->lgar_calib_params.Ksat_3;

  else if (name.compare("ponded_depth_max") == 0)
    return (void*)&this->state->lgar_calib_params.ponded_depth_max;
  else if (name.compare("field_capacity") == 0)
    return (void*)&this->state->lgar_calib_params.field_capacity_psi;
  else if (name.compare("a") == 0)
    return (void*)&this->state->lgar_calib_params.a;
  else if (name.compare("b") == 0)
    return (void*)&this->state->lgar_calib_params.b;
  else if (name.compare("frac_to_CR") == 0)
    return (void*)&this->state->lgar_calib_params.frac_to_CR;
  else if (name.compare("a_slow") == 0)
    return (void*)&this->state->lgar_calib_params.a_slow;
  else if (name.compare("b_slow") == 0)
    return (void*)&this->state->lgar_calib_params.b_slow;
  else if (name.compare("frac_slow") == 0)
    return (void*)&this->state->lgar_calib_params.frac_slow;
  else if (name.compare("spf_factor") == 0)
    return (void*)&this->state->lgar_calib_params.spf_factor;
  else {
    std::stringstream errMsg;
    errMsg << "variable "<< name << " does not exist";
    throw std::runtime_error(errMsg.str());
    return NULL;
  }
  
  // delete it later
  return NULL;
}

void BmiLGAR::
GetValueAtIndices (std::string name, void *dest, int *inds, int len)
{
  void * src = NULL;

  src = this->GetValuePtr(name);

  if (src) {
    int i;
    int itemsize = 0;
    int offset;
    char *ptr;

    itemsize = this->GetVarItemsize(name);

    for (i=0, ptr=(char *)dest; i<len; i++, ptr+=itemsize) {
      offset = inds[i] * itemsize;
      memcpy(ptr, (char *)src + offset, itemsize);
    }
  }
}


void BmiLGAR::
SetValue (std::string name, void *src)
{
  void * dest = NULL;
  dest = this->GetValuePtr(name);

  if (dest) {
    int nbytes = 0;
    nbytes = this->GetVarNbytes(name);
    memcpy(dest, src, nbytes);
  }

}


void BmiLGAR::
SetValueAtIndices (std::string name, int * inds, int len, void *src)
{
  void * dest = NULL;

  dest = this->GetValuePtr(name);

  if (dest) {
    int i;
    int itemsize = 0;
    int offset;
    char *ptr;

    itemsize = this->GetVarItemsize(name);

    for (i=0, ptr=(char *)src; i<len; i++, ptr+=itemsize) {
      offset = inds[i] * itemsize;
      memcpy((char *)dest + offset, ptr, itemsize);
    }
  }
}


std::string BmiLGAR::
GetComponentName()
{
  return "LASAM (Lumped Arid/Semi-arid Model)";
}


int BmiLGAR::
GetInputItemCount()
{
  return this->input_var_name_count;
}


int BmiLGAR::
GetOutputItemCount()
{
  return this->output_var_name_count;
}


std::vector<std::string> BmiLGAR::
GetInputVarNames()
{
  std::vector<std::string> names;

  for (int i=0; i<this->input_var_name_count; i++)
    names.push_back(this->input_var_names[i]);

  return names;
}


std::vector<std::string> BmiLGAR::
GetOutputVarNames()
{
  std::vector<std::string> names;

  for (int i=0; i<this->output_var_name_count; i++)
    names.push_back(this->output_var_names[i]);

  return names;
}


double BmiLGAR::
GetStartTime () {
  return 0.0;
}


double BmiLGAR::
GetEndTime () {
  return this->state->lgar_bmi_params.endtime_s;
}


double BmiLGAR::
GetCurrentTime () {
  return this->state->lgar_bmi_params.time_s;
}


std::string BmiLGAR::
GetTimeUnits() {
  return "s";
}


double BmiLGAR::
GetTimeStep () {
  return this->state->lgar_bmi_params.forcing_resolution_h * 3600.; // convert hours to seconds
}

std::string BmiLGAR::
GetGridType(const int grid)
{
  if (grid == 0)
    return "uniform_rectilinear";
  else
    return "";
}


void BmiLGAR::
GetGridX(const int grid, double *x)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridX: "<<grid<<" "<<x[0]<<"\n";
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridY(const int grid, double *y)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridY: "<<grid<<" "<<y[0]<<"\n";
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridZ(const int grid, double *z)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridZ: "<<grid<<" "<<z[0]<<"\n";
  throw bmi_lgar::NotImplemented();
}


int BmiLGAR::
GetGridNodeCount(const int grid)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridNodeCount: "<<grid<<"\n";
  throw bmi_lgar::NotImplemented();
}


int BmiLGAR::
GetGridEdgeCount(const int grid)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridEdgeCount: "<<grid<<"\n";
  throw bmi_lgar::NotImplemented();
}


int BmiLGAR::
GetGridFaceCount(const int grid)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridFaceCount: "<<grid<<"\n";
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridEdgeNodes(const int grid, int *edge_nodes)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridEdgeNodes: "<<grid<<" "<<edge_nodes[0]<<"\n";
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridFaceEdges(const int grid, int *face_edges)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridFaceNodes: "<<grid<<" "<<face_edges[0]<<"\n";
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridFaceNodes(const int grid, int *face_nodes)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridFaceNodes: "<<grid<<" "<<face_nodes[0]<<"\n";
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridNodesPerFace(const int grid, int *nodes_per_face)
{
  // this is not needed but printing here to avoid compiler warnings
  std::cerr<<"GetGridNodesPerFace: "<<grid<<" "<<nodes_per_face[0]<<"\n";
  throw bmi_lgar::NotImplemented();
}

#endif
