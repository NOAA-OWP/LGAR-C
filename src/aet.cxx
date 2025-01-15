#ifndef AET_CXX_INCLUDE
#define AET_CXX_INCLUDE

#include "../include/all.hxx"

//################################################################################
/* authors : Fred Ogden and Ahmad Jan and Peter La Follette
   year    : 2022
   the code computes actual evapotranspiration given PET.
   It uses an S-shaped function used in HYDRUS-1D (Simunek & Sejna, 2018).
   AET = PET * 1/(1 + (h/h_50) )^3
   h is the capillary head at the surface and
   h_50 is the capillary head at which AET = 0.5 * PET. */
//################################################################################

bool correct_AET_out_of_root_zone_via_surf_WFs = FALSE;

extern void calc_aet(bool TO_enabled, double PET_timestep_cm, double time_step_h, double wilting_point_psi_cm, double field_capacity_psi_cm, double root_zone_depth_cm, double *surf_frac_rz, 
                       int *soil_type, double AET_thresh_Theta, double AET_expon, struct wetting_front* head,
            		       struct soil_properties_ *soil_properties, double *surf_AET_vec)
{//calc_aet just calcualtes AET that should be extracted from surface WFs, GW WF contribution to AET is calculated elsewhere

  if (verbosity.compare("high") == 0) {
    printf("Computing AET... \n");
    printf("Note: AET_thresh_theta = %lf and AET_expon = %lf are not used in the computation of the current AET model. \n", AET_thresh_Theta, AET_expon);
  }

  if (!TO_enabled){
    double actual_ET_demand = 0.0;
    struct wetting_front *current;
    
    double theta_wp;
    
    double Se,theta_e,theta_r;
    double vg_a, vg_m, vg_n;
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
    double head_at_which_PET_equals_AET_cm = field_capacity_psi_cm; //340.9 is -0.33 atm, expressed in negative water depth, which is a good field capacity for most soils.
    //Coarser soils like sand will have a field capacity of -0.1 atm or so, which would be 103.3 cm expressed as a negative water pressure.
    double theta_fc = calc_theta_from_h(head_at_which_PET_equals_AET_cm, vg_a,vg_m, vg_n, theta_e, theta_r);
    
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
    
    int wetting_front_free_drainage_temp = wetting_front_free_drainage(head);
    surf_AET_vec[wetting_front_free_drainage_temp] = actual_ET_demand;
    *surf_frac_rz = 0.0;
  }

  else {//TO mode is on
    double actual_ET_demand = 0.0;

    if (listLength_surface(head)==0){
      actual_ET_demand = 0.0;
    }
    else{

      struct wetting_front *current;

      double prior_depth = 0.0;
      double *surf_thicknesses_vec = (double *)malloc(sizeof(double)*(listLength_surface(head)+1));

      current = listFindFront(listLength_TO_WFs_above_surface_WFs(head) + 1, head, NULL);
      double vector_sum = 0.0;
      double frac_of_rz_occupied_by_surf_WFs = 0.0;
      double frac_of_rz_occupied_by_top_surf_WF = 0.0;
      for (int wf = listLength_TO_WFs_above_surface_WFs(head) + 1; wf != (listLength(head)); wf++){
        int listLength_TO_WFs_above_surface_WFs_temp = listLength_TO_WFs_above_surface_WFs(head);
        if (current->is_WF_GW==TRUE){
          break;
        }

        surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] += current->depth_cm - prior_depth; 
        prior_depth = current->depth_cm;

        if (current->to_bottom==TRUE){
          surf_thicknesses_vec[wf + 1 - listLength_TO_WFs_above_surface_WFs_temp] += surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp];
          surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] = 0.0;
        }

        struct wetting_front *lowest_eligible_TO_wf;

        lowest_eligible_TO_wf = listFindFront(listLength_surface(head), head, NULL);

        while (lowest_eligible_TO_wf->next!=NULL){

          lowest_eligible_TO_wf = lowest_eligible_TO_wf->next;

          if (lowest_eligible_TO_wf->depth_cm>root_zone_depth_cm){
            break;
          }
        }

        if (correct_AET_out_of_root_zone_via_surf_WFs && listLength_surface(head)>0){
          if ((listLength_TO_WFs_in_rz_nonzero_depth(root_zone_depth_cm, head)==0) && (wf==(listLength_TO_WFs_above_surface_WFs(head) + listLength_surface(head)))){
            surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] = 0.0;
            double sum = 0.0;
            for (int ii=0; ii<listLength_surface(head); ii++) {
              sum = sum + surf_thicknesses_vec[ii];
            }
            surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] = root_zone_depth_cm - sum;
          }
        }

        if (current->depth_cm > root_zone_depth_cm){
          surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] = root_zone_depth_cm - vector_sum;
          int count = 1;
          struct wetting_front *temp_current = current;
          while(current->to_bottom==TRUE){
            surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp + count] = surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] ;
            surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] = 0.0;
            count ++;
            current = current->next;
            if (current->next==NULL){
              break;
            }
          }
          current = temp_current;
        }

        vector_sum += surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp];

        if (surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp]>root_zone_depth_cm){
          printf("surf_thicknesses_vec[wf]>root_zone_depth_cm, which is not allowed \n");
          printf("wf: %d \n", wf);
          abort();
        }

        double theta_wp;

        double head_at_which_PET_equals_AET_cm = field_capacity_psi_cm;
        
        double Se,theta_e,theta_r;
        double vg_a, vg_m, vg_n;
        int layer_num, soil_num;

        layer_num = 1;//current->layer_num; 
        // HYDRUS uses a single P50 value for whole soil profile, so for now we do that too
        soil_num  = soil_type[layer_num];
        theta_e   = soil_properties[soil_num].theta_e;
        theta_r   = soil_properties[soil_num].theta_r;
        vg_a      = soil_properties[soil_num].vg_alpha_per_cm;
        vg_m      = soil_properties[soil_num].vg_m;
        vg_n      = soil_properties[soil_num].vg_n;

        // compute theta field capacity
        double theta_fc = calc_theta_from_h(head_at_which_PET_equals_AET_cm, vg_a,vg_m, vg_n, theta_e, theta_r);
        
        double wp_head_theta = calc_theta_from_h(wilting_point_psi_cm, vg_a,vg_m, vg_n, theta_e, theta_r);
        
        theta_wp = (theta_fc - wp_head_theta)*1/2 + wp_head_theta; // theta_50 in python

        Se = calc_Se_from_theta(theta_wp,theta_e,theta_r);
        double psi_wp_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
        // if (verbosity.compare("high") == 0) {
        //   printf("this is for running with HYDRUS: this is P50 in HYDRUS \n");
        //   printf("%lf \n",-1*psi_wp_cm); //this is for running with HYDRUS: this is P50 in HYDRUS 
        // }

        double h_ratio = 1.0 + pow(current->psi_cm/psi_wp_cm, 3.0);

        frac_of_rz_occupied_by_surf_WFs += fmin(1, surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] / root_zone_depth_cm);
        
        double frac_of_rz_occupied_by_this_WF = fmin(1, surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] / root_zone_depth_cm);

        actual_ET_demand = frac_of_rz_occupied_by_this_WF * PET_timestep_cm * (1/h_ratio) * time_step_h;

        if (frac_of_rz_occupied_by_top_surf_WF == 0.0){
          frac_of_rz_occupied_by_top_surf_WF = frac_of_rz_occupied_by_this_WF;
        }

        surf_AET_vec[wf] = actual_ET_demand;

        if (verbosity.compare("high") == 0) {
          printf("\n ");
          printf("AET =  %14.10f \n",actual_ET_demand);
          printf("PET_timestep_cm*time_step_h =  %14.10f \n",PET_timestep_cm*time_step_h);
          printf("1/h_ratio =  %14.10f \n",1/h_ratio);
          printf("current->front_num: %d \n",current->front_num);
          printf("\n ");
        }
        
        double front_num_deepest_surface_WF = listLength_surface(head) + listLength_TO_WFs_above_surface_WFs(head);
        if ( front_num_deepest_surface_WF==current->front_num ){
          break;
        }

        if ((current->depth_cm > root_zone_depth_cm) && (current->to_bottom==FALSE)){
          break;
        }

        if (current->front_num==listLength(head)){
          printf("breaking because this is the deepest WF \n");
          break;
        }

        current = current->next;
      }

      *surf_frac_rz = frac_of_rz_occupied_by_surf_WFs;

      free(surf_thicknesses_vec);

    }

  }

}



extern double calc_aet_for_individual_TO_WFs(int WF_num, double WF_thickness_cm, double rooting_zone_depth, double PET_timestep_cm, double time_step_h, double wilting_point_psi_cm,
		       double field_capacity_psi_cm, int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head){
  struct wetting_front *current; 
  current = listFindFront(WF_num, *head, NULL);
  double frac = fmin(WF_thickness_cm / rooting_zone_depth, 1.0);
  
  if (verbosity.compare("high") == 0) {
    printf("frac, used in calc_aet_for_individual_TO_WFs: %lf \n",frac);
  }

  double actual_ET_demand = 0.0;
  
  double theta_wp;

  double head_at_which_PET_equals_AET_cm = field_capacity_psi_cm;
  
  double Se,theta_e,theta_r;
  double vg_a, vg_m, vg_n;
  int layer_num, soil_num;
  
  layer_num = 1;
  soil_num  = soil_type[layer_num];
  theta_e   = soil_properties[soil_num].theta_e;
  theta_r   = soil_properties[soil_num].theta_r;
  vg_a      = soil_properties[soil_num].vg_alpha_per_cm;
  vg_m      = soil_properties[soil_num].vg_m;
  vg_n      = soil_properties[soil_num].vg_n;

  double theta_fc = calc_theta_from_h(head_at_which_PET_equals_AET_cm, vg_a,vg_m, vg_n, theta_e, theta_r);
  
  double wp_head_theta = calc_theta_from_h(wilting_point_psi_cm, vg_a,vg_m, vg_n, theta_e, theta_r);
  
  theta_wp = (theta_fc - wp_head_theta)*1/2 + wp_head_theta; // theta_50 in python

  Se = calc_Se_from_theta(theta_wp,theta_e,theta_r);

  // if (verbosity.compare("high") == 0) {
  //   double p_hydrus;
  //   p_hydrus = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
  //   printf("this is P50 for HYDRUS, but in a TO WF, so possibly in a layer deeper than top: %lf \n",p_hydrus);
  //   printf("theta_fc: %lf \n",theta_fc);
  //   printf("wp_head_theta: %lf \n",wp_head_theta);
  //   printf("theta_wp: %lf \n",theta_wp);
  //   printf("Se: %lf \n",Se);
  //   printf("theta: %lf \n",current->theta);
  // }

  double psi_wp_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n); 
  double h_ratio = 1.0 + pow(current->psi_cm/psi_wp_cm, 3.0); //1 July 2024 experiment

  actual_ET_demand = PET_timestep_cm * (1/h_ratio) * time_step_h;

  if (actual_ET_demand<0.0){
    actual_ET_demand=0.0;
  }
  else if (actual_ET_demand>(PET_timestep_cm*time_step_h))
    actual_ET_demand = PET_timestep_cm*time_step_h;
  
  actual_ET_demand = actual_ET_demand * frac;

  if (current->depth_cm==0.0){//this is just to fix an apparent rounding error 
    actual_ET_demand = 0.0;
  }

  if (current->is_WF_GW==FALSE){
    actual_ET_demand = 0.0;
  }

  if (verbosity.compare("high") == 0) {
    printf("current->front_num: %d \n", current->front_num);
    printf("AET component =  %14.17f \n",actual_ET_demand);
  }

  return actual_ET_demand;

}


extern double lgarto_calc_aet_from_TO_WFs(int num_layers, double deepest_surf_depth_at_start, double root_zone_depth, double PET_timestep_cm, double timestep_h, double surf_frac_rz, 
                                          double wilting_point_psi_cm, double field_capacity_psi_cm, int *soil_type, double *cum_layer_thickness_cm, struct soil_properties_ *soil_properties, struct wetting_front** head){

  struct wetting_front *current;
  struct wetting_front *next;
  current = *head;
  next = current->next;
  for (int wf = 1; wf != (listLength(*head)-1); wf++){
    if (current->to_bottom==FALSE && next->is_WF_GW==TRUE && current->is_WF_GW==TRUE && fabs(current->psi_cm - next->psi_cm)<1.e-3){
      current = listDeleteFront(current->front_num, head, soil_type, soil_properties);
    }
    current = next;
    if (current==NULL){
      break;
    }
    next = current->next;
  }                                          

  if (verbosity.compare("high") == 0) {
    printf("calculated mass before AET extraction: %.17lf \n",lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
    listPrint(*head);
  }

  if (surf_frac_rz >= 1.0){
    return 0;
  }

  double cumulative_ET_from_TO_WFs_cm = 0.0;

  int front_num_deepest_surface_WF = listLength_surface(*head) + listLength_TO_WFs_above_surface_WFs(*head);
  if (front_num_deepest_surface_WF>0){
    surf_frac_rz = listFindFront(front_num_deepest_surface_WF, *head, NULL)->depth_cm / root_zone_depth;
  }
  else{
    surf_frac_rz = 0.0;
  }

  if (verbosity.compare("high") == 0) {
    printf("front_num_deepest_surface_WF: %d \n", front_num_deepest_surface_WF);
    printf("surf_frac_rz: %lf \n", surf_frac_rz);
  }

  double already_computed_TO_RZ_frac = 0.0;

  // find deepest non to_bottom TO WF in the RZ, should probably do this before loop
  struct wetting_front *deepest_AET_accepting_TO_WF_in_RZ = *head;
  struct wetting_front *temp_WF = *head;

  while (temp_WF->to_bottom==TRUE || temp_WF->is_WF_GW==FALSE || temp_WF->depth_cm==0.0 || (temp_WF->to_bottom==FALSE && temp_WF->is_WF_GW==TRUE && temp_WF->depth_cm>0.0) ){
    if ((temp_WF->to_bottom==FALSE && temp_WF->is_WF_GW==TRUE && temp_WF->depth_cm>0.0)){
      deepest_AET_accepting_TO_WF_in_RZ = temp_WF;
    }
    temp_WF = temp_WF->next;
    if (temp_WF==NULL){
      break;
    }
    if (temp_WF->depth_cm>=root_zone_depth){
      break;
    }
  }
  double deepest_AET_accepting_TO_WF_in_RZ_depth_before_move = deepest_AET_accepting_TO_WF_in_RZ->depth_cm;

  if ( (root_zone_depth>0) && (surf_frac_rz<(1-1e-5) ) ){

    int WF_num = 1; 
    double WF_thickness_cm = 0.0; 
    double temp_AET_value = 0.0;
    double thickness_of_WF_crossing_rooting_zone = 0.0;
    double *thicknesses_vec = (double *)malloc(sizeof(double)*(listLength(*head)+1)); //the +1 is to avoid the error "malloc: Heap corruption detected, free list is damaged"
    double frac_vec_sum = 0.0; 

    thicknesses_vec[0] = 0.0;

    /////for AET in TO stuff, first we make a vector thicknesses_vec that will be used to determine the fraction of AET each WF in the rooting zone contributes. 
    for (int wf = 1; wf != 1+(listLength(*head)); wf++) { 

      WF_num = wf;
      struct wetting_front *current = listFindFront(WF_num, *head, NULL);

      if (current->depth_cm >= root_zone_depth){
        break;
      }

      if (current->is_WF_GW==FALSE){
        continue;
      }

      if (current->next!=NULL){
        if (current->next->is_WF_GW==FALSE){
          continue;
        }
      }

      if (current->depth_cm == 0.0){
        continue;
      }

      if (wf==1){
        WF_thickness_cm = current->depth_cm;
      }
      else{
        WF_thickness_cm = current->depth_cm - listFindFront(WF_num-1, *head, NULL)->depth_cm;
      }

      if ( (current->to_bottom==TRUE) && (current->depth_cm <= root_zone_depth)){
        WF_thickness_cm = 0.0;
      }

      struct wetting_front *previous = listFindFront(WF_num - 1, *head, NULL);

      if (previous!=NULL){
        while(previous->to_bottom == TRUE){
          struct wetting_front *previous_previous = listFindFront(previous->front_num - 1, *head, NULL);
          if (previous_previous!=NULL){
            WF_thickness_cm += previous->depth_cm - previous_previous->depth_cm;
          }
          else{
            WF_thickness_cm += previous->depth_cm;
          }
          previous = previous_previous;
          if (previous==NULL){
            break;
          }
        }
      }


      if (deepest_AET_accepting_TO_WF_in_RZ!=NULL){
        if (deepest_AET_accepting_TO_WF_in_RZ->depth_cm<=root_zone_depth && deepest_AET_accepting_TO_WF_in_RZ->depth_cm!=0.0 && current->front_num==deepest_AET_accepting_TO_WF_in_RZ->front_num){
          WF_thickness_cm += root_zone_depth - deepest_AET_accepting_TO_WF_in_RZ->depth_cm;
        }
      }

      if (wf==listLength(*head)){
        WF_thickness_cm = (1 - frac_vec_sum - surf_frac_rz) * root_zone_depth;
      }
      else if ( current->depth_cm >= root_zone_depth ){//should never happen now because we break if depth>root zone depth earlier 
        WF_thickness_cm = (1 - frac_vec_sum - surf_frac_rz) * root_zone_depth;
        if ( (WF_num-1)>0){
          thickness_of_WF_crossing_rooting_zone = current->depth_cm - listFindFront(WF_num-1, *head, NULL)->depth_cm;
        }
        else{
          thickness_of_WF_crossing_rooting_zone = current->depth_cm;
        }
      }

      frac_vec_sum += (WF_thickness_cm / root_zone_depth);

      if (verbosity.compare("high") == 0) {
        printf("frac_vec_sum: %lf \n",frac_vec_sum); 
      }

      thicknesses_vec[WF_num] = WF_thickness_cm;

      if (thicknesses_vec[0]>1e10){
        printf("thicknesses_vec[0]>1e10 \n");
        abort();
      }

      if (verbosity.compare("high") == 0) {
        printf("thicknesses_vec[WF_num]: %.17lf \n",thicknesses_vec[WF_num]);
      }

      if (current->depth_cm > root_zone_depth){
        break;
      }
      
    }


    /////second for the AET in TO WFs process, we calculate the contribution each WF in the rooting zone yields for the AET, and then the depths of these WFs are adjusted, and the AET flux is adjusted too
    int second_to_last_thicknesses_vec_index_yielding_nonzero_result = 0;
    int last_thicknesses_vec_index_yielding_nonzero_result = 0;

    int n = listLength(*head);

    for (int ii = 0; ii<n; ii++){
      if (thicknesses_vec[ii] != 0.0){
        second_to_last_thicknesses_vec_index_yielding_nonzero_result = last_thicknesses_vec_index_yielding_nonzero_result; //not correctly named, 
        last_thicknesses_vec_index_yielding_nonzero_result = ii;
      }
    }

    for (int ii = 0; ii<last_thicknesses_vec_index_yielding_nonzero_result; ii++){
      already_computed_TO_RZ_frac = already_computed_TO_RZ_frac + thicknesses_vec[ii];
    }
    already_computed_TO_RZ_frac = already_computed_TO_RZ_frac / root_zone_depth;

    if (thicknesses_vec[second_to_last_thicknesses_vec_index_yielding_nonzero_result]>1e10){//TODO: fix error where second_to_last_thicknesses_vec_index_yielding_nonzero_result is an index greater than the length of thicknesses_vec
      thicknesses_vec[second_to_last_thicknesses_vec_index_yielding_nonzero_result] = 0.0;//I think this code can be deleted
    }

    if ( (root_zone_depth>0.0) ){
      for (int wf = 1; wf != (listLength(*head)); wf++) { 
        WF_num = wf;
        struct wetting_front *current = listFindFront(WF_num, *head, NULL);

        if (current->is_WF_GW==FALSE){
          continue;
        }

        if (fabs(current->depth_cm)>cum_layer_thickness_cm[num_layers]){
          break;
        }

        double depth_at_start = current->depth_cm;

        if (depth_at_start>root_zone_depth){
          break;
        }

        WF_thickness_cm = thicknesses_vec[WF_num];
        if (verbosity.compare("high") == 0) {
          printf("WF_thickness_cm directly from vec: %lf \n",WF_thickness_cm);
          printf("front_num: %d \n", current->front_num);
        }

        if ( (WF_thickness_cm>0.0) && (current->to_bottom==FALSE) ){
          temp_AET_value = calc_aet_for_individual_TO_WFs(WF_num, WF_thickness_cm, root_zone_depth, PET_timestep_cm, timestep_h, wilting_point_psi_cm, field_capacity_psi_cm,
                soil_type, soil_properties, head);

          if (verbosity.compare("high") == 0) {
            printf("temp_AET_value: %lf \n",temp_AET_value); 
            printf("WF_thickness_cm: %lf \n",WF_thickness_cm);
          }

          struct wetting_front *previous = listFindFront(current->front_num - 1, *head, NULL); //6 September 2023 addition to move seemingly stuck fronts
          int count = 1;
          if (previous!=NULL){
            while (previous->to_bottom==1){
              previous = listFindFront(current->front_num - count, *head, NULL);
              count ++;
              if (previous==NULL){
                break;
              }
            }
          }

          if (previous!=NULL){
            if ( (current->depth_cm>root_zone_depth) ){

                current->depth_cm += temp_AET_value / (current->next->theta - current->theta); 
                if (previous->depth_cm < 0.0){
                  double overshot_height = fabs(previous->depth_cm);
                  previous->depth_cm = 0.0;
                  cumulative_ET_from_TO_WFs_cm -= (previous->next->theta - previous->theta)*overshot_height;
                }
            }
            else{
              if (thickness_of_WF_crossing_rooting_zone==0){
                current->depth_cm += temp_AET_value / (current->next->theta - current->theta); 
              }
              else{
                current->depth_cm += temp_AET_value / (current->next->theta - current->theta); 
              }
            }
          }
          else{
            if (thickness_of_WF_crossing_rooting_zone==0.0){
              current->depth_cm += temp_AET_value / (current->next->theta - current->theta); 
            }
            else{
              current->depth_cm += temp_AET_value / (current->next->theta - current->theta); 
            }
          }

          cumulative_ET_from_TO_WFs_cm += temp_AET_value;
        }

      if ( (depth_at_start > root_zone_depth) && (cumulative_ET_from_TO_WFs_cm!=0.0) ){
        break;
      }

      }
    }

    if (verbosity.compare("high") == 0) {
      if (temp_AET_value>0){
        printf("WFs after TO WF depths updated via AET extraction: \n");
        listPrint(*head);
      }
    }

    free(thicknesses_vec);
  }
  if (verbosity.compare("high") == 0) {
    printf("cumulative_ET_from_TO_WFs_cm: %lf \n", cumulative_ET_from_TO_WFs_cm);
  }

  if ( (cumulative_ET_from_TO_WFs_cm<1e-8) && (root_zone_depth>0.0) && (PET_timestep_cm>0.0) && (surf_frac_rz<1.0) ){
    struct wetting_front *current = *head;
    struct wetting_front *next = current->next;
    for (int wf = 1; wf != (listLength(*head)); wf++) {
      if ( (current->depth_cm<1e-8) && (next->depth_cm>0.0) && (next->to_bottom==TRUE) && (current->is_WF_GW==TRUE) && (next->is_WF_GW==TRUE) ){
        current->depth_cm += 0.0000001;
        cumulative_ET_from_TO_WFs_cm += 0.0000001*(next->theta - current->theta);
        break;
      }
      current = next;
      next = current->next;
    }
  }


  bool run_AET_correction = FALSE;
  if (deepest_AET_accepting_TO_WF_in_RZ==NULL){
    run_AET_correction = TRUE;
  }
  if (deepest_AET_accepting_TO_WF_in_RZ_depth_before_move>root_zone_depth || deepest_AET_accepting_TO_WF_in_RZ_depth_before_move==0.0 || deepest_AET_accepting_TO_WF_in_RZ->to_bottom==TRUE || deepest_AET_accepting_TO_WF_in_RZ->is_WF_GW==FALSE){
    run_AET_correction = TRUE;
  }
  if (correct_AET_out_of_root_zone_via_surf_WFs){
    run_AET_correction = FALSE;
    if (listLength_surface(*head)==0){
      if (deepest_AET_accepting_TO_WF_in_RZ==NULL){
        run_AET_correction = TRUE;
      }
      if (deepest_AET_accepting_TO_WF_in_RZ_depth_before_move>root_zone_depth || deepest_AET_accepting_TO_WF_in_RZ_depth_before_move==0.0 || deepest_AET_accepting_TO_WF_in_RZ->to_bottom==TRUE || deepest_AET_accepting_TO_WF_in_RZ->is_WF_GW==FALSE){
        run_AET_correction = TRUE;
      }
    }
  }


  // if (listLength_TO_WFs_in_rz_nonzero_depth(root_zone_depth, *head)>0){
    // if ((surf_frac_rz<0.99999) && (PET_timestep_cm>0.0) && (cumulative_ET_from_TO_WFs_cm<1e-6)) {
  if (run_AET_correction){
    if ((surf_frac_rz<0.99999) && (PET_timestep_cm>0.0)) {
      struct wetting_front *shallowest_eligible_TO_WF = listFindFront(listLength_TO_WFs_above_surface_WFs(*head) + listLength_surface(*head) + 1, *head, NULL);
      struct wetting_front *next = shallowest_eligible_TO_WF->next;


      if (listLength_surface(*head)==0){
        while ((shallowest_eligible_TO_WF->next!=NULL) && (shallowest_eligible_TO_WF->depth_cm<root_zone_depth)){
            shallowest_eligible_TO_WF = shallowest_eligible_TO_WF->next;
        }
      }
      while(shallowest_eligible_TO_WF->to_bottom==TRUE){
        if (shallowest_eligible_TO_WF->next==NULL){
          break;
        }
        shallowest_eligible_TO_WF = shallowest_eligible_TO_WF->next;
      }
      if (shallowest_eligible_TO_WF->next==NULL){
        next = shallowest_eligible_TO_WF;
      }
      else{
        next = shallowest_eligible_TO_WF->next;
      }
      int WF_num = shallowest_eligible_TO_WF->front_num;
      double WF_thickness_cm = fmax(0,root_zone_depth*(1-deepest_surf_depth_at_start/root_zone_depth-already_computed_TO_RZ_frac));

      double temp_AET_value = calc_aet_for_individual_TO_WFs(WF_num, WF_thickness_cm, root_zone_depth, PET_timestep_cm, timestep_h, wilting_point_psi_cm, field_capacity_psi_cm,
                  soil_type, soil_properties, head);

        if (verbosity.compare("high") == 0) {
          printf("shallowest_eligible_TO_WF->front_num: %d \n", shallowest_eligible_TO_WF->front_num);
          printf("shallowest_eligible_TO_WF->depth_cm: %.17lf \n", shallowest_eligible_TO_WF->depth_cm);
          printf("mass at start of TO AET extraction below rz: %lf \n", lgar_calc_mass_bal(cum_layer_thickness_cm, *head));
          printf("temp_AET_value: %lf \n", temp_AET_value);
          listPrint(*head);
        }

        double mass_before_AET_extraction = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        if (shallowest_eligible_TO_WF->layer_num==num_layers+1000){//...the idea here is to temporarily turn off this code, because when working with TO WFs, the "simple case" is the lowest layer not the highest layer. Anyway the code in the else clause should work generally. 
          shallowest_eligible_TO_WF->theta -= temp_AET_value / WF_thickness_cm;
          int soil_num = soil_type[num_layers];
          double theta_r = soil_properties[soil_num].theta_r;
          if (shallowest_eligible_TO_WF->theta<theta_r){
            double temp_mass_before = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
            shallowest_eligible_TO_WF->theta = theta_r + 1.E-12;
            double temp_mass_after = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
            temp_AET_value = temp_AET_value - (temp_mass_after - temp_mass_before);
          }
        }
        else{
          int layer_num = shallowest_eligible_TO_WF->layer_num;
          int soil_num = soil_type[layer_num];
          double vg_a_k, vg_m_k, vg_n_k;
          double theta_e_k, theta_r_k;

          double *delta_thetas    = (double *)malloc(sizeof(double)*(layer_num+1));
          double *delta_thickness = (double *)malloc(sizeof(double)*(layer_num+1));

          double psi_cm = shallowest_eligible_TO_WF->psi_cm;
          double psi_cm_below = next->psi_cm;

          double new_mass = (shallowest_eligible_TO_WF->depth_cm - cum_layer_thickness_cm[layer_num-1]) * (shallowest_eligible_TO_WF->theta - 0*next->theta);
          double prior_mass = new_mass; 

          // compute mass in the layers above the current wetting front
          // use the psi of the current wetting front and van Genuchten parameters of
          // the respective layers to get the total mass above the current wetting front

          double AET_demand_cm;
            AET_demand_cm = temp_AET_value;

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

              new_mass += layer_thickness * (theta - 0*theta_below);
            prior_mass += layer_thickness * (theta - 0*theta_below);

            delta_thetas[k] = 0*theta_below; 
            delta_thickness[k] = layer_thickness;

          }

          delta_thetas[layer_num] = 0*next->theta;
          delta_thickness[layer_num] = shallowest_eligible_TO_WF->depth_cm - cum_layer_thickness_cm[layer_num-1];
          
          new_mass = prior_mass;
          prior_mass = prior_mass - AET_demand_cm;
          
          double theta_new = 0.0;
          if (listLength_surface(*head)==0){
            theta_new = lgar_theta_mass_balance(layer_num, soil_num, shallowest_eligible_TO_WF->psi_cm, new_mass, prior_mass, &AET_demand_cm, 
                      delta_thetas, delta_thickness, soil_type, soil_properties);
            shallowest_eligible_TO_WF->theta = theta_new;
          }
          else {

            double psi_cm_loc = psi_cm; // location psi
            double tolerance = 1.E-12;
            new_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
            prior_mass = new_mass - AET_demand_cm;
            double psi_cm_loc_prev = psi_cm_loc;
            double delta_mass = fabs(new_mass - prior_mass);
            double theta = 0.0;
            int iter = 0;
            double factor = 1.0;
            bool switched = false; // flag that determines capillary head to be incremented or decremented
            bool iter_aug_flag=FALSE;
            double delta_mass_prev   = delta_mass;
            int count_no_mass_change = 0;
            int break_no_mass_change = 5;


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
                
              }

              theta = calc_theta_from_h(psi_cm_loc, soil_properties[soil_num].vg_alpha_per_cm, soil_properties[soil_num].vg_m,
                      soil_properties[soil_num].vg_n,soil_properties[soil_num].theta_e,
                      soil_properties[soil_num].theta_r);
              

              shallowest_eligible_TO_WF->psi_cm = psi_cm_loc;
              shallowest_eligible_TO_WF->theta = theta;
              theta_new = theta;

              struct wetting_front *current = listFindFront(listLength(*head), *head, NULL);
              struct wetting_front *above = listFindFront(listLength(*head)-1, *head, NULL);

              for (int k=listLength(*head); k!=0; k--){
                if (above!=NULL){
                  if (above->to_bottom==TRUE){
                    above->psi_cm = current->psi_cm;

                    int soil_num_k  = soil_type[above->layer_num];
                    
                    double theta_e_k = soil_properties[soil_num_k].theta_e;
                    double theta_r_k = soil_properties[soil_num_k].theta_r;
                    double vg_a_k    = soil_properties[soil_num_k].vg_alpha_per_cm;
                    double vg_m_k    = soil_properties[soil_num_k].vg_m;
                    double vg_n_k    = soil_properties[soil_num_k].vg_n;
                    above->theta = calc_theta_from_h(above->psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
                  }
                }
                current = above;
                above = listFindFront(current->front_num - 1, *head, NULL);
                if (above==NULL){
                  break;
                }
              }

              new_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
              delta_mass = fabs(new_mass - prior_mass);



              if (fabs(psi_cm_loc - psi_cm_loc_prev) < 1E-15 && factor < 1E-13) break;

              // another condition to avoid infinite loop when the error does not improve
              if (fabs(delta_mass - delta_mass_prev) < 1E-15)
                count_no_mass_change++;
              else
                count_no_mass_change = 0;

              // break the loop if the mass does not change in the five consecutive iterations.
              if (count_no_mass_change == break_no_mass_change)
                break;
              
              if (psi_cm_loc > 1e7){//there are rare cases where theta is very close to theta_r, and delta_mass - delta_mass_prev will change extremely slowly. Convergence might be possible but the model will take hours to converge rather than seconds. 
              //an alternative solution was to change the threshold in if (fabs(delta_mass - delta_mass_prev) < 1e-15) to 1e-11, but that solution is somewhat slow. 
                break;
              }

              // -ve pressure will return NAN, so terminate the loop if previous psi is way small and current psi is zero
              // the wetting front is almost saturated
              if (psi_cm_loc <= 0 && psi_cm_loc_prev < 0) break;

              delta_mass_prev = delta_mass;


            }
          }

          free(delta_thetas);
          free(delta_thickness);
        }


        int soil_num  = soil_type[shallowest_eligible_TO_WF->layer_num];
        
        double theta_e = soil_properties[soil_num].theta_e;
        double theta_r = soil_properties[soil_num].theta_r;
        double vg_a    = soil_properties[soil_num].vg_alpha_per_cm;
        double vg_m    = soil_properties[soil_num].vg_m;
        double vg_n    = soil_properties[soil_num].vg_n;
        double Se = calc_Se_from_theta(shallowest_eligible_TO_WF->theta,theta_e,theta_r);
        shallowest_eligible_TO_WF->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
        
        struct wetting_front *current = listFindFront(listLength(*head), *head, NULL);
        struct wetting_front *above = listFindFront(listLength(*head)-1, *head, NULL);

        for (int k=listLength(*head); k!=0; k--){
          if (above!=NULL){
            if (above->to_bottom==TRUE){
              above->psi_cm = current->psi_cm;

              int soil_num_k  = soil_type[above->layer_num];
              
              double theta_e_k = soil_properties[soil_num_k].theta_e;
              double theta_r_k = soil_properties[soil_num_k].theta_r;
              double vg_a_k    = soil_properties[soil_num_k].vg_alpha_per_cm;
              double vg_m_k    = soil_properties[soil_num_k].vg_m;
              double vg_n_k    = soil_properties[soil_num_k].vg_n;
              above->theta = calc_theta_from_h(above->psi_cm, vg_a_k, vg_m_k, vg_n_k, theta_e_k, theta_r_k);
            }
          }
          current = above;
          if (current!=NULL){
            above = listFindFront(current->front_num - 1, *head, NULL);
          }
          if (above==NULL){
            break;
          }
        }

        double mass_after_AET_extraction = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
        cumulative_ET_from_TO_WFs_cm += (mass_before_AET_extraction - mass_after_AET_extraction);

    }//makes it so that water is extracted from below the RZ under certain circumstances 
  }

  return(cumulative_ET_from_TO_WFs_cm);

}


#endif
