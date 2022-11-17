#ifndef AET_CXX_INCLUDE
#define AET_CXX_INCLUDE

#include "../include/all.hxx"

//################################################################################
/* authors : Fred Ogden and Ahmad Jan
   year    : 2022
   the code computes actual evapotranspiration given PET.
   It uses an S-shaped function used in HYDRUS-1D (Simunek & Sejna, 2018).
   AET = PET * 1/(1 + (h/h_50) )^3
   h is the capillary head at the surface and
   h_50 is the capillary head at which AET = 0.5 * PET. */
//################################################################################


extern double calc_aet(double PET_timestep_cm, double time_step_h, double wilting_point_psi_cm,
		       int *soil_type, double AET_thresh_Theta, double AET_expon,
		       struct soil_properties_ *soil_properties)
{

  if (verbosity.compare("high") == 0) {
    printf("Computing AET... \n");
  }
  
  double actual_ET_demand=0.0;
  struct wetting_front *current;
  double wettest_theta_in_root_zone = 0.0;  // initial value of max() trap function
  
  double wilting_point = wilting_point_psi_cm;
  double Se_wp, theta_wp;

  double relative_moisture_at_which_PET_equals_AET = 0.75;
  
  double theta,Se,theta_e,theta_r;
  double vg_a, vg_m, vg_n;
  double psi_cm;
  int layer_num, soil_num;
  
  current = head;

  layer_num = current->layer_num;
  soil_num  = soil_type[layer_num];
  theta_e   = soil_properties[soil_num].theta_e;
  theta_r   = soil_properties[soil_num].theta_r;
  vg_a      = soil_properties[soil_num].vg_alpha_per_cm;
  vg_m      = soil_properties[soil_num].vg_m;
  vg_n      = soil_properties[soil_num].vg_n;

  // compute theta field capacity
  double theta_fc = (theta_e - theta_r) * relative_moisture_at_which_PET_equals_AET + theta_r;
  
  double wp_head_theta = calc_theta_from_h(wilting_point_psi_cm, vg_a,vg_m, vg_n, theta_e, theta_r);
  
  theta_wp = (theta_fc - wp_head_theta)*1/2 + wp_head_theta; // theta_50 in python

  Se = calc_Se_from_theta(theta_wp,theta_e,theta_r);
  double psi_wp_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);

  double h_ratio = 1.0 + pow(current->psi_cm/psi_wp_cm, 3.0);

  actual_ET_demand = PET_timestep_cm * (1/h_ratio) * time_step_h;

  if (actual_ET_demand<0)
    actual_ET_demand=0.0;
  else if (actual_ET_demand>(PET_timestep_cm*time_step_h))
    actual_ET_demand = PET_timestep_cm*time_step_h;

  if (verbosity.compare("high") == 0) {
    printf("AET =  %14.10f \n",actual_ET_demand);
  }
  
  return actual_ET_demand;

}

#endif
