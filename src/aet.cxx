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
{

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
    double head_at_which_PET_equals_AET_cm = field_capacity_psi_cm; //340.9 is 0.33 atm, expressed in water depth, which is a good field capacity for most soils.
    //Coarser soils like sand will have a field capacity of 0.1 atm or so, which would be 103.3 cm.
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
    
    // return actual_ET_demand;
    int wetting_front_free_drainage_temp = wetting_front_free_drainage(head);
    // surf_AET_vec[0] = actual_ET_demand;
    surf_AET_vec[wetting_front_free_drainage_temp] = actual_ET_demand;
    *surf_frac_rz = 0.0;
  }

  else {//TO mode is on
    double actual_ET_demand = 0.0;

    if (listLength_surface(head)==0){
      actual_ET_demand = 0.0;
      // printf("listLength_surface(head)==0 \n");
    }
    else{

      struct wetting_front *current;

      double prior_depth = 0.0;
      double *surf_thicknesses_vec = (double *)malloc(sizeof(double)*(listLength_surface(head)+1));

      current = listFindFront(listLength_TO_WFs_above_surface_WFs(head) + 1, head, NULL);
      double vector_sum = 0.0;
      double frac_of_rz_occupied_by_surf_WFs = 0.0;
      double frac_of_rz_occupied_by_top_surf_WF = 0.0;
      // printf("loop to calc AET for surf WFs starts at: listLength_TO_WFs_above_surface_WFs(head) + 1: %d \n", listLength_TO_WFs_above_surface_WFs(head) + 1);
      // printf("loop ends at: %d \n", (listLength(head)));
      for (int wf = listLength_TO_WFs_above_surface_WFs(head) + 1; wf != (listLength(head)); wf++){
        // printf("in loop that does surf_thicknesses_vec and surf AET stuff. On WF number: %d \n", wf);
        int listLength_TO_WFs_above_surface_WFs_temp = listLength_TO_WFs_above_surface_WFs(head);
        if (current->is_WF_GW==TRUE){
          // printf("breaking because current->is_WF_GW==TRUE \n");
          break;
        }

        surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] += current->depth_cm - prior_depth; 
        prior_depth = current->depth_cm;

        if (current->to_bottom==TRUE){
          surf_thicknesses_vec[wf + 1 - listLength_TO_WFs_above_surface_WFs_temp] += surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp];
          surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] = 0.0;
        }

        // if (soil_properties[soil_type[1]].Ksat_cm_per_h>0.5){
        struct wetting_front *lowest_eligible_TO_wf;

        lowest_eligible_TO_wf = listFindFront(listLength_surface(head), head, NULL);

        while (lowest_eligible_TO_wf->next!=NULL){

          lowest_eligible_TO_wf = lowest_eligible_TO_wf->next;

          if (lowest_eligible_TO_wf->depth_cm>root_zone_depth_cm){
            break;
          }
        }






        if (correct_AET_out_of_root_zone_via_surf_WFs && listLength_surface(head)>0){
          // printf("listLength_surface(head): %d \n", listLength_surface(head));
          // if ((listLength_TO_WFs_in_rz_nonzero_depth(root_zone_depth_cm, head)==0) && (wf==(listLength_TO_WFs_above_surface_WFs(head) + listLength_surface(head)))){
          if ((listLength_TO_WFs_in_rz_nonzero_depth(root_zone_depth_cm, head)==0) && (wf==(listLength_TO_WFs_above_surface_WFs(head) + listLength_surface(head)))){
              // maybe if AET from surf WFs>0? or perhaps this condition in the other code that extracts from below RZ? just checking if surf_frac > 0 in the other code might not work when a new surf WF was created that has no AET component. And in fact, that already might make that calc a bit off ... 
          // if (wf==(listLength_TO_WFs_above_surface_WFs(head) + listLength_surface(head))){
            surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] = 0.0;
            double sum = 0.0;
            // for (int ii=0; ii<=sizeof(surf_thicknesses_vec)/sizeof(surf_thicknesses_vec[0]); ii++) { //This code, which was in LGARTO, does not correctly determine the length of surf_thicknesses_vec. However, I think it should just be listlengthsurface or listlengthsurface+1
            for (int ii=0; ii<listLength_surface(head); ii++) {
              sum = sum + surf_thicknesses_vec[ii];
            }
            surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] = root_zone_depth_cm - sum;
            // surf_thicknesses_vec[wetting_front_free_drainage(head)] = root_zone_depth_cm - sum;
            // printf("surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp]: %lf \n", surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp]);
          }//makes it so that if there are surface WFs that do not span the entire RZ but there are no eligible TO WFs in the RZ, the deepest surface WF takes the remainder of the rz
        // } 

          // if ((listLength_TO_WFs_in_rz_nonzero_depth(root_zone_depth_cm, head)>0) && (wf==(listLength_TO_WFs_above_surface_WFs(head) + listLength_surface(head)))){
          //     // maybe if AET from surf WFs>0? or perhaps this condition in the other code that extracts from below RZ? just checking if surf_frac > 0 in the other code might not work when a new surf WF was created that has no AET component. And in fact, that already might make that calc a bit off ... 
          // // if (wf==(listLength_TO_WFs_above_surface_WFs(head) + listLength_surface(head))){
          //   surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] = 0.0;
          //   double sum = 0.0;
          //   double sum_next = 0.0;
          //   int count_index = 0;
          //   // for (int ii=0; ii<=sizeof(surf_thicknesses_vec)/sizeof(surf_thicknesses_vec[0]); ii++) { //This code, which was in LGARTO, does not correctly determine the length of surf_thicknesses_vec. However, I think it should just be listlengthsurface or listlengthsurface+1
          //   for (int ii=0; ii<=listLength_surface(head); ii++) {
          //     sum = sum + surf_thicknesses_vec[ii];
          //     count_index = ii;
          //   }

          //   sum_next = sum + surf_thicknesses_vec[count_index+1];
              
          //   printf("sum, which should be equal to the depth of surf WFs above: %lf \n", sum);
          //   printf("frac_rz_TO_below_lowest_eligible_TO_wf, which should be equal to the frac of the RZ that is TO: %lf \n", frac_rz_TO_below_lowest_eligible_TO_wf);
          //   // surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] = root_zone_depth_cm*frac_rz_TO_below_lowest_eligible_TO_wf - sum;
          //   surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] = sum_next - sum;
          //   printf("surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp]: %lf \n", surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp]);
          // }//makes it so that

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
              // printf("breaking because current->next==NULL \n");
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

        layer_num = 1;//current->layer_num; //5 Septmber 2023: trying to suss out issue with too large recharge (negative recharge issue)
        // layer_num = current->layer_num;
        // //issue was that HYDRUS uses a single P50 value for whole soil profile, so for now we do that too
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
        if (verbosity.compare("high") == 0) {
          printf("this is for running with HYDRUS: this is P50 in HYDRUS \n");
          printf("%lf \n",-1*psi_wp_cm); //this is for running with HYDRUS: this is P50 in HYDRUS 
          // abort();
        }

        double h_ratio = 1.0 + pow(current->psi_cm/psi_wp_cm, 3.0);

        frac_of_rz_occupied_by_surf_WFs += fmin(1, surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] / root_zone_depth_cm);
        // printf("frac_of_rz_occupied_by_surf_WFs: %lf \n", frac_of_rz_occupied_by_surf_WFs);
        // printf("surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp]: %lf \n", surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp]);
        
        double frac_of_rz_occupied_by_this_WF = fmin(1, surf_thicknesses_vec[wf - listLength_TO_WFs_above_surface_WFs_temp] / root_zone_depth_cm);

        // actual_ET_demand = frac_of_rz_occupied_by_surf_WFs * PET_timestep_cm * (1/h_ratio) * time_step_h; //names are slightly off but are correct: PET_timestep_cm in this case has units of cm/h, so actual_ET_demand will have units of cm
        actual_ET_demand = frac_of_rz_occupied_by_this_WF * PET_timestep_cm * (1/h_ratio) * time_step_h;
        if (actual_ET_demand>1e10){//remove before commit
          printf("actual_ET_demand: %lf \n", actual_ET_demand);
          printf("frac_of_rz_occupied_by_this_WF: %lf \n", frac_of_rz_occupied_by_this_WF);
          printf("PET_timestep_cm: %lf \n", PET_timestep_cm);
          printf("1/h_ratio: %lf \n", 1/h_ratio);
          printf("time_step_h: %lf \n", time_step_h);
          abort();
        }
        if (frac_of_rz_occupied_by_top_surf_WF == 0.0){
          frac_of_rz_occupied_by_top_surf_WF = frac_of_rz_occupied_by_this_WF;
        }

        surf_AET_vec[wf] = actual_ET_demand;
        // printf("\n");
        // printf("actual_ET_demand as computed in calc_aet: %.17lf \actual_ET_demand", actual_ET_demand);
        // printf("\n");

        // printf("wf: %d \n", wf);
        // printf("associated actual_ET_demand: %.17lf \n", actual_ET_demand);

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
          // printf("breaking because front_num_deepest_surface_WF==current->front_num \n");
          // printf("should probably break here for front_num_deepest_surface_WF>=current->front_num instead \n");
          break;
        }

        if ((current->depth_cm > root_zone_depth_cm) && (current->to_bottom==FALSE)){
          // printf("breaking because current->depth_cm > root_zone_depth_cm \n");
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
  // layer_num = current->layer_num; 
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

  if (verbosity.compare("high") == 0) {
    double p_hydrus;
    p_hydrus = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
    printf("this is P50 for HYDRUS, but in a TO WF, so possibly in a layer deeper than top: %lf \n",p_hydrus);
    // abort();
    printf("theta_fc: %lf \n",theta_fc);
    printf("wp_head_theta: %lf \n",wp_head_theta);
    printf("theta_wp: %lf \n",theta_wp);
    printf("Se: %lf \n",Se);
    printf("theta: %lf \n",current->theta);
  }

  
  double psi_wp_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n); 
  double h_ratio = 1.0 + pow(current->psi_cm/psi_wp_cm, 3.0); //1 July 2024 experiment
  // double psi_to_use = 0.0;
  // if (current->next==NULL){
  //   psi_to_use = current->psi_cm;
  // }
  // else{
  //   psi_to_use = current->next->psi_cm;
  // }
  // double h_ratio = 1.0 + pow(psi_to_use/psi_wp_cm, 3.0);

  // printf("psi_wp_cm: %lf \n", psi_wp_cm);

  actual_ET_demand = PET_timestep_cm * (1/h_ratio) * time_step_h;

  // printf("PET_timestep_cm =  %14.17f \n", PET_timestep_cm);
  // printf("1/h_ratio =  %14.17f \n", 1/h_ratio);
  // printf("time_step_h =  %14.17f \n", time_step_h);

  if (actual_ET_demand<0.0){
    actual_ET_demand=0.0;
  }
  else if (actual_ET_demand>(PET_timestep_cm*time_step_h))
    actual_ET_demand = PET_timestep_cm*time_step_h;

  // printf("frac =  %14.17f \n", frac);
  
  actual_ET_demand = actual_ET_demand * frac;

  if (current->depth_cm==0.0){//this is just to fix an apparent rounding error 
    actual_ET_demand = 0.0; ///ok, but now we allow AET from 1 specific depth = 0 TO WF ...
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
      listDeleteFront(current->front_num, head, soil_type, soil_properties);
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

  // printf("surf_frac_rz before aug: %lf \n", surf_frac_rz);
  // printf("deepest_surf_depth_at_start: %lf \n", deepest_surf_depth_at_start);

  if (surf_frac_rz >= 1.0){
    return 0;
  }

  double cumulative_ET_from_TO_WFs_cm = 0.0;

  int front_num_deepest_surface_WF = listLength_surface(*head) + listLength_TO_WFs_above_surface_WFs(*head);
  if (front_num_deepest_surface_WF>0){
    surf_frac_rz = listFindFront(front_num_deepest_surface_WF, *head, NULL)->depth_cm / root_zone_depth;
  }
  else{
    surf_frac_rz = 0.0;//17 June 2024 add
  }

  if (verbosity.compare("high") == 0) {
    printf("front_num_deepest_surface_WF: %d \n", front_num_deepest_surface_WF);
    printf("surf_frac_rz: %lf \n", surf_frac_rz);
  }

  // printf("surf_frac_rz after aug: %lf \n", surf_frac_rz);

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
  }//2 July AET fix
  double deepest_AET_accepting_TO_WF_in_RZ_depth_before_move = deepest_AET_accepting_TO_WF_in_RZ->depth_cm;

  // printf("deepest_AET_accepting_TO_WF_in_RZ->front_num: %d \n", deepest_AET_accepting_TO_WF_in_RZ->front_num);

  if ( (root_zone_depth>0) && (surf_frac_rz<(1-1e-5) ) ){

    int WF_num = 1; //17 June 2024: added = 0.0
    double WF_thickness_cm = 0.0; //17 June 2024: added = 0.0
    double temp_AET_value = 0.0;
    double thickness_of_WF_crossing_rooting_zone = 0.0;
    double *thicknesses_vec = (double *)malloc(sizeof(double)*(listLength(*head)+1)); //the +1 is to avoid the error "malloc: Heap corruption detected, free list is damaged"
    double frac_vec_sum = 0.0; 

    thicknesses_vec[0] = 0.0;

    // printf("WF_thickness_cm: %lf \n", WF_thickness_cm);

    /////for AET in TO stuff, first we make a vector thicknesses_vec that will be used to determine the fraction of AET each WF in the rooting zone contributes. 
    for (int wf = 1; wf != 1+(listLength(*head)); wf++) { 

      WF_num = wf;
      struct wetting_front *current = listFindFront(WF_num, *head, NULL);

      if (current->depth_cm >= root_zone_depth){
        break;//2 July 2024 addition 
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
        // printf("WF_thickness_cm calc 1: %lf \n", WF_thickness_cm);
      }
      else{
        WF_thickness_cm = current->depth_cm - listFindFront(WF_num-1, *head, NULL)->depth_cm;
        // printf("WF_thickness_cm calc 2: %lf \n", WF_thickness_cm);
      }

      if ( (current->to_bottom==TRUE) && (current->depth_cm <= root_zone_depth)){
        WF_thickness_cm = 0.0;
        // printf("WF_thickness_cm calc 3: %lf \n", WF_thickness_cm);
      }

      struct wetting_front *previous = listFindFront(WF_num - 1, *head, NULL);

      if (previous!=NULL){
        while(previous->to_bottom == TRUE){
          struct wetting_front *previous_previous = listFindFront(previous->front_num - 1, *head, NULL);
          if (previous_previous!=NULL){
            WF_thickness_cm += previous->depth_cm - previous_previous->depth_cm;
            // printf("WF_thickness_cm calc 4: %lf \n", WF_thickness_cm);
          }
          else{
            WF_thickness_cm += previous->depth_cm;
            // printf("WF_thickness_cm calc 5: %lf \n", WF_thickness_cm);
          }
          previous = previous_previous;
          if (previous==NULL){
            break;
          }
        }
      }







      //find the remainder of the RZ thickness below that WF 
      //If there are non to bottom TO WFs in the RZ, then add that thickness to the deepest. If not, do the theta correction code. 
      if (deepest_AET_accepting_TO_WF_in_RZ!=NULL){
        if (deepest_AET_accepting_TO_WF_in_RZ->depth_cm<=root_zone_depth && deepest_AET_accepting_TO_WF_in_RZ->depth_cm!=0.0 && current->front_num==deepest_AET_accepting_TO_WF_in_RZ->front_num){
          WF_thickness_cm += root_zone_depth - deepest_AET_accepting_TO_WF_in_RZ->depth_cm;
        }
      }//2 July AET fix





      if (wf==listLength(*head)){
        WF_thickness_cm = (1 - frac_vec_sum - surf_frac_rz) * root_zone_depth;
        // printf("WF_thickness_cm calc 6: %lf \n", WF_thickness_cm);
      }
      else if ( current->depth_cm >= root_zone_depth ){//should never happen now because we break if depth>root zone depth earlier 
        WF_thickness_cm = (1 - frac_vec_sum - surf_frac_rz) * root_zone_depth;
        // printf("WF_thickness_cm calc 7: %lf \n", WF_thickness_cm);
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

      // if (WF_thickness_cm<0.0){//this was commented out in LGARTO code
      //   printf("WF_thickness_cm<0.0, which is not allowed. Issue might be that surf_frac_rz is set before surface WF movement, and then surf_frac_rz is used in lgarto_calc_aet_from_TO_WFs after surface WF movement. Solution might be to just redefine surf_frac_rz in lgarto_calc_aet_from_TO_WFs \n ");
      //   abort();
      // }

      thicknesses_vec[WF_num] = WF_thickness_cm;

      if (thicknesses_vec[0]>1e10){
        printf("thicknesses_vec[0]>1e10 \n");
        abort();
      }

      // printf("setting thicknesses_vec[WF_num] as WF_thickness_cm, printing WF_thickness_cm: %lf \n", WF_thickness_cm);
      // printf("and that was for WF: %d \n", WF_num);

      if (verbosity.compare("high") == 0) {
        printf("thicknesses_vec[WF_num]: %.17lf \n",thicknesses_vec[WF_num]);
      }

      if (current->depth_cm > root_zone_depth){
        break;
      }
      
    }

/////commenting out because a value of 0 in TO frac is now allowed is now allowed 
    // if ( ((frac_vec_sum+surf_frac_rz)<0.999) || ((frac_vec_sum+surf_frac_rz)>1.001) ){// just because frac_vec_sum should be 1, but it can for example be something like 1 + 1e-9 as well
    //   if (){
    //     printf("fractions of thicknesses of WFs contributing to AET do not sum to 1 \n");
    //     printf("frac_vec_sum: %10.20lf \n",frac_vec_sum);
    //     printf("surf_frac_rz: %10.20lf \n", surf_frac_rz);
    //     listPrint(*head);
    //     abort();
    //   }
    // } //changing to check if the sum of thicknesses in the root zone divided by the root zone depth is 1, because now we take ET from surface and GW WFs at the same time 




    /////second for the AET in TO WFs process, we calculate the contribution each WF in the rooting zone yields for the AET, and then the depths of these WFs are adjusted, and the AET flux is adjusted too
    int second_to_last_thicknesses_vec_index_yielding_nonzero_result = 0;
    int last_thicknesses_vec_index_yielding_nonzero_result = 0;

    int n = listLength(*head);

    for (int ii = 0; ii<n; ii++){
      if (thicknesses_vec[ii] != 0.0){
        second_to_last_thicknesses_vec_index_yielding_nonzero_result = last_thicknesses_vec_index_yielding_nonzero_result; //not correctly named, 
        last_thicknesses_vec_index_yielding_nonzero_result = ii;
        // second_to_last_thicknesses_vec_index_yielding_nonzero_result = ii-1;
      }
    }

    for (int ii = 0; ii<last_thicknesses_vec_index_yielding_nonzero_result; ii++){
      already_computed_TO_RZ_frac = already_computed_TO_RZ_frac + thicknesses_vec[ii];
    }
    already_computed_TO_RZ_frac = already_computed_TO_RZ_frac / root_zone_depth;

    // printf("thicknesses_vec[0] before aug: %.27lf \n", thicknesses_vec[0]);
    // printf("thicknesses_vec[1] before aug: %.27lf \n", thicknesses_vec[1]);
    // printf("thicknesses_vec[2] before aug: %.27lf \n", thicknesses_vec[2]);
    // printf("thicknesses_vec[3] before aug: %.27lf \n", thicknesses_vec[3]);
    // printf("thicknesses_vec[4] before aug: %.27lf \n", thicknesses_vec[4]);


    // // if (){
    //   thicknesses_vec[second_to_last_thicknesses_vec_index_yielding_nonzero_result] += thicknesses_vec[last_thicknesses_vec_index_yielding_nonzero_result];
    //   thicknesses_vec[last_thicknesses_vec_index_yielding_nonzero_result] = 0.0;
    // // }

    // int second_to_last = -1;

    // for (int ii = 1; ii<n; ii++){
    //   if ((thicknesses_vec[ii] != 0.0) && (ii!=last_thicknesses_vec_index_yielding_nonzero_result)){
    //     second_to_last = ii;
    //   }
    // }

    // if (thicknesses_vec[1]==0.0 && thicknesses_vec[2]==0.0 && listFindFront(1, *head, NULL)->depth_cm==0.0 && listFindFront(2, *head, NULL)->depth_cm==0.0){
    //   second_to_last = 2;
    // }

    // printf("second_to_last: %d \n", second_to_last);

    // if (second_to_last<n && second_to_last>-1){
    //   while(listFindFront(second_to_last, *head, NULL)->next->depth_cm==0.0){
    //     printf("second_to_last: %d \n", second_to_last);
    //     if (listFindFront(second_to_last + 1, *head, NULL)->next==NULL){
    //       break;
    //     }
    //     second_to_last = second_to_last + 1;
    //   }
    // }

    // printf("thicknesses_vec[second_to_last] before aug: %lf \n", thicknesses_vec[second_to_last]);

    // if (second_to_last>-1){
    //   thicknesses_vec[second_to_last] += thicknesses_vec[last_thicknesses_vec_index_yielding_nonzero_result];
    //   thicknesses_vec[last_thicknesses_vec_index_yielding_nonzero_result] = 0.0;
    // } //26 June 2024 removal 

    if (thicknesses_vec[second_to_last_thicknesses_vec_index_yielding_nonzero_result]>1e10){//TODO: fix error where second_to_last_thicknesses_vec_index_yielding_nonzero_result is an index greater than the length of thicknesses_vec
      thicknesses_vec[second_to_last_thicknesses_vec_index_yielding_nonzero_result] = 0.0;//I think this code can be deleted
    }

    // printf("last_thicknesses_vec_index_yielding_nonzero_result: %d \n", last_thicknesses_vec_index_yielding_nonzero_result);

    // printf("thicknesses_vec[0] after aug: %.17lf \n", thicknesses_vec[0]);
    // printf("thicknesses_vec[1] after aug: %.17lf \n", thicknesses_vec[1]);
    // printf("thicknesses_vec[2] after aug: %.17lf \n", thicknesses_vec[2]);
    // printf("thicknesses_vec[3] after aug: %.17lf \n", thicknesses_vec[3]);
    // printf("thicknesses_vec[4] after aug: %.17lf \n", thicknesses_vec[4]);

    
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

        // printf("WF_thickness_cm: %.17lf \n", WF_thickness_cm);

        // printf("\n ###################################################################################################################################################################################### \n");
        // int size = sizeof(thicknesses_vec) / sizeof(thicknesses_vec[0]);
        // printf("size: %d \n", size);
        // for(int i = 0; i < size+10; i++) {
        //   printf("Element of thicknesses_vec %d: %.17lf\n", i, thicknesses_vec[i]);
        // }
        // printf("\n ###################################################################################################################################################################################### \n");

        if ( (WF_thickness_cm>0.0) && (current->to_bottom==FALSE) ){
          // printf("tag 1 \n");
          temp_AET_value = calc_aet_for_individual_TO_WFs(WF_num, WF_thickness_cm, root_zone_depth, PET_timestep_cm, timestep_h, wilting_point_psi_cm, field_capacity_psi_cm,
                soil_type, soil_properties, head);

          if (verbosity.compare("high") == 0) {
            printf("temp_AET_value: %lf \n",temp_AET_value); 
            printf("WF_thickness_cm: %lf \n",WF_thickness_cm);
          }

          // if (WF_thickness_cm>1e10){
          //   printf("WF_thickness_cm>1e10 \n");
          //   // abort();
          // }

          // printf("temp_AET_value: %.17lf \n", temp_AET_value);
          // printf("and that is for WF_num: %d \n", WF_num);
          // cumulative_ET_from_TO_WFs_cm += temp_AET_value; //8 July 2024: moving to after height is subtracted from TO WFs, to address marginal case where AET has to be augmented
          // printf("updated cumulative_ET_from_TO_WFs_cm: %.17lf \n", cumulative_ET_from_TO_WFs_cm);


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
                // current->theta -= temp_AET_value / WF_thickness_cm; // was previously extracting AET by changing theta, now doing it based on height.
                current->depth_cm += temp_AET_value / (current->next->theta - current->theta); 
                // printf("new depth: %lf \n", current->depth_cm);
              }
              else{
                // current->theta -= temp_AET_value / thickness_of_WF_crossing_rooting_zone;
                current->depth_cm += temp_AET_value / (current->next->theta - current->theta); 
                // printf("new depth: %lf \n", current->depth_cm);
              }
            }
          }
          else{
            if (thickness_of_WF_crossing_rooting_zone==0.0){
              // current->theta -= temp_AET_value / WF_thickness_cm; // was previously extracting AET by changing theta, now doing it based on height. 
              current->depth_cm += temp_AET_value / (current->next->theta - current->theta); 
              // printf("new depth: %lf \n", current->depth_cm);
            }
            else{
              // current->theta -= temp_AET_value / thickness_of_WF_crossing_rooting_zone;
              current->depth_cm += temp_AET_value / (current->next->theta - current->theta); 
              // printf("new depth: %lf \n", current->depth_cm);
            }
          }

          // if (current->depth_cm > cum_layer_thickness_cm[num_layers]){
          //   current->depth_cm = cum_layer_thickness_cm[num_layers];
          // }

          // if (current->depth_cm > cum_layer_thickness_cm[num_layers] && current->next!=NULL){//trying to eliminate huge movements when two adjacent TO WFs have very close theta vals, 3 July 2024
          //   double AET_adj = (current->next->theta - current->theta)*(current->depth_cm - cum_layer_thickness_cm[num_layers]);
          //   current->depth_cm = cum_layer_thickness_cm[num_layers];
          //   temp_AET_value -= AET_adj;
          // }

          cumulative_ET_from_TO_WFs_cm += temp_AET_value;
        }

      // if ( (current->depth_cm > root_zone_depth) && (cumulative_ET_from_TO_WFs_cm!=0.0) ){
        // break; //25 June 2024: not necessary I think, and also causes an error when AET causes a WF to get quite deep in 1 time step. Maybe try breaking in another way ...
      // }

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

    free(thicknesses_vec); //should this type of thing be more common for efficiency? 
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




  // if (listLength_surface(*head)==0){
  //   struct wetting_front *current = *head;
  //   while (current->next->depth_cm==0.0){
  //     if (current->next==NULL){
  //       break;
  //     }
  //     current = current->next;
  //   }
  //   if (current->depth_cm==0.0 && current->is_WF_GW==TRUE){
  //     current->depth_cm = current->depth_cm + 1.E-9;
  //   }
  // }




  // if (soil_properties[soil_type[1]].Ksat_cm_per_h<=0.5){
  // if (surf_frac_rz==0.0){


  bool run_AET_correction = FALSE; //2 July AET fix
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
      // printf("shallowest_eligible_TO_WF->front_num second time: %d \n", shallowest_eligible_TO_WF->front_num);
      // printf("aet check 2 \n");
      while(shallowest_eligible_TO_WF->to_bottom==TRUE){
        // printf("aet check 3 \n");
        if (shallowest_eligible_TO_WF->next==NULL){
          break;
        }
        shallowest_eligible_TO_WF = shallowest_eligible_TO_WF->next;
        // printf("aet check 3.1 \n");
      }
      // printf("shallowest_eligible_TO_WF->front_num third time: %d \n", shallowest_eligible_TO_WF->front_num);
      if (shallowest_eligible_TO_WF->next==NULL){
        next = shallowest_eligible_TO_WF;
      }
      else{
        next = shallowest_eligible_TO_WF->next;
      }
      // printf("aet check 4 \n");
      // printf("aet check 5 \n");
      int WF_num = shallowest_eligible_TO_WF->front_num;
      // printf("aet check 6 \n");

      // double already_computed_TO_RZ_frac = 
      // printf("already_computed_TO_RZ_frac: %lf \n", already_computed_TO_RZ_frac);
      // printf("surf_frac_rz: %lf \n", surf_frac_rz);
      // double WF_thickness_cm = fmax(0,root_zone_depth*(1-surf_frac_rz-already_computed_TO_RZ_frac));
      double WF_thickness_cm = fmax(0,root_zone_depth*(1-deepest_surf_depth_at_start/root_zone_depth-already_computed_TO_RZ_frac));
      
      // printf("WF_thickness_cm: %.17lf \n", WF_thickness_cm);
      // printf("aet check 7 \n");
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
          // printf("theta_new in extracting below RZ, top layer: %lf \n", shallowest_eligible_TO_WF->theta);
        }
        else{
          // printf("theta_new in extracting below RZ, not top layer: %lf \n", shallowest_eligible_TO_WF->theta);
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

          // double new_mass = mass_before_AET_extraction;
          // double prior_mass = mass_before_AET_extraction;
          double AET_demand_cm;
          // if (surf_frac_rz>0.0){
            AET_demand_cm = temp_AET_value;//*(1/surf_frac_rz);//26 June 2024: added *(1/surf_frac_rz)
          // }
          // else{
          //   AET_demand_cm = temp_AET_value;//26 June 2024: added *(1/surf_frac_rz)
          // }

          // printf("checck 1 \n");


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

          // new_mass = fabs(new_mass);
          // prior_mass = fabs(new_mass);
          
          new_mass = prior_mass;
          prior_mass = prior_mass - AET_demand_cm;
          // printf("right before lgar_theta_mass_balance, here is new_mass: %lf \n", new_mass);
          // printf("right before lgar_theta_mass_balance, here is prior_mass: %lf \n", prior_mass);
          
          double theta_new = 0.0;
          if (listLength_surface(*head)==0){
            theta_new = lgar_theta_mass_balance(layer_num, soil_num, shallowest_eligible_TO_WF->psi_cm, new_mass, prior_mass, &AET_demand_cm, 
                      delta_thetas, delta_thickness, soil_type, soil_properties);
            shallowest_eligible_TO_WF->theta = theta_new;
          }
          else {







            
            // double psi_cm = shallowest_eligible_TO_WF->psi_cm;
            double psi_cm_loc = psi_cm; // location psi
            double tolerance = 1.E-12;
            new_mass = lgar_calc_mass_bal(cum_layer_thickness_cm, *head);
            prior_mass = new_mass - AET_demand_cm;
            // printf("right before lgar_theta_mass_balance like code when surf_WFs, here is new_mass: %lf \n", new_mass);
            // printf("right before lgar_theta_mass_balance like code when surf_WFs, here is prior_mass: %lf \n", prior_mass);
            double psi_cm_loc_prev = psi_cm_loc;
            double delta_mass = fabs(new_mass - prior_mass);
            double theta = 0.0;
            int iter = 0;
            double factor = 1.0;//was 1.0 previously. This code is far faster and seems to avoid loops with >10000 iterations. Low van Genuchten n values can cause this
            bool switched = false; // flag that determines capillary head to be incremented or decremented
            // printf("delta_mass: %lf \n", delta_mass);
            bool iter_aug_flag=FALSE;
            double delta_mass_prev   = delta_mass;
            int count_no_mass_change = 0;
            int break_no_mass_change = 5;


            while (delta_mass > tolerance) {
              iter++;
              // if (iter>10000){
              //   printf("iter: %d \n", iter);
              //   listPrint(*head);
              //   abort();
              // }

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

            //     if (psi_cm_loc<0.0){
            //       wanted_to_saturate_flag = TRUE;
            //     }
                
            //     if (psi_cm_loc < 0 && psi_cm_loc_prev != 0) {
            // /* this is for the extremely rare case when iterative psi_cm_loc calculation temporarily
            //   yields a negative value and the actual answer for psi_cm_loc is nonzero. For example
            //   when a completely saturated wetting front with a tiny amount of ET should yield a resulting
            //   theta that is slightly below saturation. */
            //       psi_cm_loc = psi_cm_loc_prev * 0.1;
            //     }
                
              }

              // double theta_layer;
              // double mass_layers= 0.0;

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








              // mass_layers += delta_thickness[layer_num] * (theta - delta_theta[layer_num]);

              // for (int k=1; k<layer_num; k++) {
              //   int soil_num_loc =  soil_type[k]; // _loc denotes the variable is local to the loop

              //   theta_layer = calc_theta_from_h(psi_cm_loc, soil_properties[soil_num_loc].vg_alpha_per_cm,
              //           soil_properties[soil_num_loc].vg_m, soil_properties[soil_num_loc].vg_n,
              //           soil_properties[soil_num_loc].theta_e, soil_properties[soil_num_loc].theta_r);

              //   mass_layers += delta_thickness[k] * (theta_layer - delta_theta[k]);
              // }

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




          // printf("theta_new in extracting below RZ, multilayer: %lf \n", theta_new);
          // printf("was AET_demand_cm changed? here it is: %lf \n", AET_demand_cm);
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
        // printf("cumulative_ET_from_TO_WFs_cm before below rz aug: %.17lf \n", cumulative_ET_from_TO_WFs_cm);
        cumulative_ET_from_TO_WFs_cm += (mass_before_AET_extraction - mass_after_AET_extraction);//5 June 2024 takeout
        // printf("mass_before_AET_extraction: %.17lf \n", mass_before_AET_extraction);
        // printf("mass_after_AET_extraction: %.17lf \n", mass_after_AET_extraction);
        // printf("cumulative_ET_from_TO_WFs_cm after below rz aug: %.17lf \n", cumulative_ET_from_TO_WFs_cm);

    }//makes it so that water is extracted from below the RZ under certain circumstances 
  }
  

  return(cumulative_ET_from_TO_WFs_cm);

}


#endif
