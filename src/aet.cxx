#ifndef AET_CXX_INCLUDE
#define AET_CXX_INCLUDE

#include "../include/all.hxx"

extern double calc_aet(double PET_timestep_cm, double time_step_h, double wilting_point_psi_cm,
                       struct soil_properties_ *soil_properties, int *soil_type, double AET_thresh_Theta, double AET_expon)
{

  double actual_ET_demand=0.0;
  struct wetting_front *current;
  //double PET_timestep_cm;
  double wettest_theta_in_root_zone = 0.0;  // initial value of max() trap function
  int layer = -1;  // the layer in the root zone containing the wettest soil
  int soil = -1;   // the soil type in the root zone layer with the wettest soil
  
  double wilting_point= wilting_point_psi_cm;
  double Se_wp, theta_wp;

  double relative_moisture_at_which_PET_equals_AET = 0.75;
  
  double theta,Se,theta_e,theta_r;
  double vg_a, vg_m, vg_n;
  double psi_cm, K_cm_per_s;
  int layer_num, soil_num;
  
  //PET_timestep_cm = PET_mm_per_15min*time_step_s/forcing_resolution_s/10.0;  // 10 converts from mm to cm
  
  //printf("PET_timestep = %lf %lf \n", PET_timestep_cm, wilting_point_psi_cm);
  

  
  // go through all the wetting fronts to find the wettest theta in the layer with the lowest capillary head
  current = head;

  layer_num           = current->layer_num;
  soil_num            = soil_type[layer_num];
  theta_e             = soil_properties[soil_num].theta_e;
  theta_r             = soil_properties[soil_num].theta_r;
  vg_a                = soil_properties[soil_num].vg_alpha_per_cm;
  vg_m                = soil_properties[soil_num].vg_m;
  vg_n                = soil_properties[soil_num].vg_n;
  theta = current->theta;
    
  double theta_fc = (theta_e-theta_r)*relative_moisture_at_which_PET_equals_AET+theta_r;
  //printf("theta_fc = %lf \n", theta_fc);
  double wp_head_theta = calc_theta_from_h(wilting_point_psi_cm, vg_a,vg_m, vg_n, theta_e, theta_r);
  
  theta_wp = (theta_fc-wp_head_theta)*1/2+wp_head_theta; // theta_50 in python
  //printf("theta_50 = %lf %lf %lf \n", theta_wp, wilting_point_psi_cm,wp_head_theta);

  Se = calc_Se_from_theta(theta_wp,theta_e,theta_r);
  double psi_wp_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
  //printf("psi_50 = %lf %lf \n", psi_wp_cm,current->psi_cm);

  double h_ratio = 1+pow(current->psi_cm/psi_wp_cm, 3.0);

  actual_ET_demand = PET_timestep_cm*(1/h_ratio)*time_step_h;
  
  //printf("AET.. = %lf %lf %lf \n", PET_timestep_cm, actual_ET_demand,(1/h_ratio), time_step_h);

  if (actual_ET_demand<0)
    actual_ET_demand=0.0;
  else if (actual_ET_demand>(PET_timestep_cm*time_step_h))
    actual_ET_demand = PET_timestep_cm*time_step_h;

  //printf("AET.... = %0.6e %lf %lf \n", actual_ET_demand,(1/h_ratio), time_step_h);		  
  return actual_ET_demand;

}

#endif
