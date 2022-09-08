#include "../include/all.h"

#define DEBUG 0
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

// calculate precipitation mass that will be added to the layers
extern double potential_infiltration()
{
  // in the main -- this should only be called for the condition (self.precip_data[i]>0)&(self.precip_data[i-1]>0) 

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

#if DEBUG > 1
  printf("wettign_front_free_drainage, layer_num = %d %d \n", wf_that_supplies_free_drainage_demand, current->layer_num);
#endif
}

extern void lgar_move_fronts(double *ponded_depth_cm, double time_step_s, int wf_free_drainage_demand, double old_mass, int num_layers , double *AET_demand_cm, double *cum_layer_thickness_cm, int *soil_type, struct soil_properties_ *soil_properties)
{

#if DEBUG > 1
  printf("wetting fronts before moving \n");
  listPrint();
#endif
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
  
  for (int l= number_of_wetting_fronts; l != 0; l--) {

#if DEBUG > 1
    printf("Looping over wetting fronst = %d %d \n", l);
#endif
    
    if (l == 1 && number_of_wetting_fronts >0) {
      current = listFindFront(l,NULL);
      next = listFindFront(l+1,NULL);
      previous = NULL;
      
      current_old = listFindFront(l,state_previous);
      next_old = listFindFront(l+1,state_previous);
      previous_old = NULL;
    }
    else if (l <number_of_wetting_fronts) {
      current = listFindFront(l,NULL);
      next = listFindFront(l+1,NULL);
      previous = listFindFront(l-1,NULL);
      
      current_old = listFindFront(l,state_previous);
      next_old = listFindFront(l+1,state_previous);
      previous_old = listFindFront(l-1,state_previous);
    }
    else if (l==number_of_wetting_fronts) {
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

#if DEBUG > 1
    printf ("****************** Cases ***************** \n");
    printf ("Layers (wf, layer, above, below) == %d %d %d %d \n", l ,layer_num, layer_number_above, layer_number_below);
#endif
    
    double free_drainage_demand = 0.0;
    double actual_ET_demand = *AET_demand_cm;
    
    // case to check if the wetting front is at the interface, i.e. deepest wetting front within a layer
    // todo. this condition can be replace by current->to_depth = FALSE && l<last_wetting_front_index
    /*************************************************************************************/
    if ( (l<last_wetting_front_index) && (layer_number_below!=layer_num) ) {
#if DEBUG > 1
      printf("case (deepest wetting front) : layer_num_below (%d) != layer_num (%d) \n", layer_num, layer_number_below);
#endif
      
      current->theta = calc_theta_from_h(next->psi_cm, vg_a,vg_m, vg_n, theta_e, theta_r);
      current->psi_cm = next->psi_cm;
    }

    // case to check if the number of wetting fronst are equal to the number of layers, i.e., one wetting front per layer
    /*************************************************************************************/
    
    if (l == number_of_wetting_fronts && layer_number_below != layer_num && number_of_wetting_fronts == num_layers) {
      
#if DEBUG > 1
      printf("case (number_of_wetting_fronts equal to num_layers) : l (%d) == num_layers (%d) == num_wetting_fronts(%d) \n", l, num_layers,number_of_wetting_fronts);
#endif
      
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
      
#if DEBUG > 1
      printf ("prior_mass_layer = %6.12f \n", prior_mass);
#endif
      
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
      
      //double xx = lgar_calc_mass_bal(0,cum_layer_thickness_cm);
      delta_thetas[layer_num] = 0.0; 
      delta_thickness[layer_num] = current->depth_cm - cum_layer_thickness_cm[layer_num-1];
      
      double free_drainage_demand = 0;
      
      if (wf_free_drainage_demand == l)
	prior_mass += precip_mass_to_add - (free_drainage_demand + actual_ET_demand);
      
      
      double theta_new = lgar_theta_mass_balance_2(layer_num, soil_num, psi_cm, new_mass, prior_mass, current->depth_cm, delta_thetas, delta_thickness, soil_type, soil_properties);
      
      current->theta = fmin(theta_new, theta_e);
      
      double Se = calc_Se_from_theta(current->theta,theta_e,theta_r);
      current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
      
    }


    // case to check if the wetting fronst is within the layer
    /*************************************************************************************/
    
    if ( (l < last_wetting_front_index) && (layer_number_below == layer_num) ) {
      
#if DEBUG > 1
      printf("case (wetting front within layer) : layer_num (%d) == layer_num_below (%d) \n", layer_num,layer_number_below);
#endif
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
	
	double theta_new = lgar_theta_mass_balance_2(layer_num, soil_num, psi_cm, new_mass, prior_mass, current->depth_cm, delta_thetas, delta_thickness, soil_type, soil_properties);
	
	current->theta = fmin(theta_new, theta_e);
	
      }
      
      
      double Se = calc_Se_from_theta(current->theta,theta_e,theta_r);
      current->psi_cm = calc_h_from_Se(Se, vg_a, vg_m, vg_n);
    }
    
    
#if DEBUG > 1
    printf("*********** Cases for mass balance of wetting fronts is done!! ************** \n");
#endif
    
    // if f_p (predicted infiltration) causes theta > theta_e, mass correction is needed. depth of the wetting front with theta > theta_e is increased to close the mass balance
    // this should be moved out of here to a subroutine; add a call to that subroutine
    if (l == 1) {
      int soil_num_k1  = soil_type[wf_free_drainage_demand];
      double theta_e_k1 = soil_properties[soil_num_k1].theta_e;
      
      struct wetting_front *wf_free_drainage = listFindFront(wf_free_drainage_demand,NULL);
      
      double mass_timestep = (old_mass + precip_mass_to_add) - (actual_ET_demand+free_drainage_demand);
      
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
    //    lgar_merge_wetting_fronts()
#if DEBUG > 1
    if (next != NULL) {
      printf("********** Merging wetting fronts ********** \n");
      listPrint();
    }
    else
      printf("********** No merging is needed ********** \n");
#endif

    // THIS WILL NEED SOME MAJOR CLEANUP or a separate module if needed
    
    double column_depth = cum_layer_thickness_cm[num_layers];
    
    if (next != NULL) {

      struct wetting_front *next_to_next = listFindFront(l+2,NULL); // next->next;

      
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

#if DEBUG > 1
	printf ("Deleting wetting front (before)... \n");
	listPrint();
#endif
	listDeleteFront(next->front_num);

#if DEBUG > 1
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
	//printf("A3: merging time ............ %lf \n",cum_layer_thickness_cm[layer_num]);
	// wetting front moving within a layer; update values
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
	*ponded_depth_cm = (current->theta - next->theta) *  (current->depth_cm - next->depth_cm);
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


    }
   
  }
  /*******************************************************************/
  // end of the for loop
  /*******************************************************************/

  
  // reset current to head to fix any mass balance issues and dry-over-wet wetting fronts conditions
  current = head;
  next = current->next;

  double mass_change = 0.0;
  
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
	mass_change += fabs(mass_after - mass_before);

	// note: mass_before is less when we have wetter front over drier front condition, however, lgar_calc_mass_bal returns mass_before > mass_after due to fabs(theta_current - theta_next); for mass_before the functions compuates more than the actual mass; removing fabs in that function might be one option, but for now we are adding fabs to mass_change to make sure we added extra water back to AET after deleting the drier front
	//printf("Mass change = %6.15f %6.15f %6.15f \n", mass_change, mass_after, mass_before);
	
      }
      current = current->next;
      next = current->next;
    }
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
  //printf("Mass Balance = %6.12f %6.12f %6.12f %6.12f\n" , mass_before_move, mass_after_move, *AET_demand_cm, precip_mass_to_add);

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
  
}


extern double lgar_insert_water(double *ponded_depth_cm, double *volin_this_timestep, double precip_timestep_cm, double dry_depth, int nint, double time_step_s, int wf_free_drainage_demand,
                              int *soil_type, struct soil_properties_ *soil_properties, double *cum_layer_thickness_cm)
{
  //printf("inserting water.... %lf \n",*ponded_depth_cm);

  int wf_that_supplies_free_drainage_demand = wf_free_drainage_demand;
  // This subroutine inserts precipitation into the soil
  // note ponded_depth_cm is a pointer.   Access it's value as (*ponded_depth_cm).

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

  double hp_cm_max=0.0; //h_p_max = 0.0;

  double h_p = fmax(*ponded_depth_cm - precip_timestep_cm * time_step_s,0.0); // water ponded on the surface
  // printf("P1 = %lf %lf %lf \n", h_p, *ponded_depth_cm, precip_timestep_cm * time_step_s);
  current = head;
  current_free_drainage = listFindFront(wf_that_supplies_free_drainage_demand,NULL);
  current_free_drainage_next = listFindFront(wf_that_supplies_free_drainage_demand+1,NULL);

  /*  
  if(current->to_bottom == TRUE )  
    {
      // first layer does not have an infiltration front, add one
      printf("aborting in insert water \n"); abort();
      if(debug_flag) printf("dry_depth= %lf\n",dry_depth);
      
      if(dry_depth * delta_theta > (*ponded_depth_cm))  // all the rain will fit in
	
	listInsertFirst((*ponded_depth_cm)/delta_theta, theta, 1, 1, FALSE);
      
      current= head;
      layer_num          = current->layer_num;
      soil_num           = soil_type[layer_num];
      theta              = soil_properties[soil_num].theta_e;
      theta_r           = soil_properties[soil_num].theta_r;
      vg_a               = soil_properties[soil_num].vg_alpha_per_cm;
      vg_m               = soil_properties[soil_num].vg_m;
      vg_n               = soil_properties[soil_num].vg_n;
      Se                 = calc_Se_from_theta(theta, theta_e, theta_r);
      Ksat_cm_per_s  = soil_properties[soil_num].Ksat_cm_per_s;
      current->psi_cm    = calc_h_from_Se(Se, vg_a, vg_m, vg_n);   
      current->K_cm_per_s = calc_K_from_Se(Se, Ksat_cm_per_s, vg_m);
    }
  
*/
  int number_of_wetting_fronts = listLength();
  double time_step = time_step_s;
  
  //printf ("Number of WF: %d \n", number_of_wetting_fronts);

  int l = number_of_wetting_fronts;
  int last_wetting_front_index = number_of_wetting_fronts;
  int num_layers = 3; // hacked
  int layer_num_fp = current_free_drainage->layer_num;

  
  double Geff;
  
  if (number_of_wetting_fronts == num_layers)
    Geff = 0.0;
  else {
   
    //double wf_below_free_drainage_one = 0.0;
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
  
    //Geff = calc_Geff(theta1, theta2, theta_e, theta_r, vg_alpha_per_cm, vg_n, vg_m, h_min_cm, Ksat_cm_per_s, nint);
    //printf("values1 = %d %lf  %lf %lf %lf \n", soil_num, theta, theta_below, theta_e, theta_r);
    //printf("values2 = %lf  %lf %lf %lf %lf %d \n", vg_a, vg_n, vg_m, h_min_cm*10, Ksat_cm_per_s*10, nint);
    Geff = calc_Geff(theta_below, theta_e, theta_e, theta_r, vg_a, vg_n, vg_m, h_min_cm, Ksat_cm_per_s, nint);
    //printf("Geff = %lf \n", Geff*10);
  }
   
  if (layer_num_fp == 1) {
    //    f_p = Ksat_cm_per_s * (1 + (Geff+ 1*(*ponded_depth_cm))/current_free_drainage->depth_cm);

    f_p = Ksat_cm_per_s * (1 + (Geff+ h_p)/current_free_drainage->depth_cm);
    
    //printf("A1 = %lf %lf %lf %lf %lf \n", f_p*10, Ksat_cm_per_s*10, Geff*10, current_free_drainage->depth_cm, (*ponded_depth_cm));
  }
  else {
    
    double bottom_sum = (current_free_drainage->depth_cm-cum_layer_thickness_cm[layer_num_fp-1])/Ksat_cm_per_s;
    double xx = bottom_sum;
    //double denominatorA = bottom_sum;
    //printf("layerXX = %d %d \n", layer_num, layer_num_fp);
    for (int k = 1; k < layer_num_fp; k++) {
      int soil_numA = soil_type[layer_num_fp-k];
      double Ksat_cm_per_s1 = soil_properties[soil_numA].Ksat_cm_per_s;
      
      bottom_sum += (cum_layer_thickness_cm[layer_num_fp - k] - cum_layer_thickness_cm[layer_num_fp - (k+1)])/ Ksat_cm_per_s1;
      //printf("check = %lf \n", (cum_layer_thickness_cm[layer_num_fp - k] - cum_layer_thickness_cm[layer_num_fp - (k-1)]));
      //printf("aborting here \n"); abort();
    }
    
    //f_p = (current_free_drainage->depth_cm + Geff + 1*(*ponded_depth_cm)) / bottom_sum;

    f_p = (current_free_drainage->depth_cm + Geff + h_p) / bottom_sum;
    //printf ("X0 = %lf \n",(current_free_drainage->depth_cm + Geff + 0*(*ponded_depth_cm)));
    //printf("X00 = %lf %lf \n", xx, bottom_sum-xx);
    //printf("f_p in loop = %6.8f %lf %6.8f %lf \n", f_p, current_free_drainage->depth_cm , Geff, bottom_sum);
  }
  
  // checkpoint # AJ
  int soil_numB = soil_type[head->layer_num];
  double theta_e1 = soil_properties[soil_numB].theta_e; // saturated theta of top layer
  
  // what is 2 in the python code??
  if ((layer_num_fp==num_layers) && (current_free_drainage->theta == theta_e1) && (num_layers == number_of_wetting_fronts)) 
    f_p = 0.0;

  double ponded_depth_temp = *ponded_depth_cm;

  //double theta_value_for_free_drainage = self.current_states[-1][1];
  //double free_drainage_demand = temp_derivs[-1][0]*config.time_step*theta_value_for_free_drainage
    
  double free_drainage_demand = 0;
  
  // checkpoint # AJ
  if ((layer_num_fp==num_layers) && (num_layers == number_of_wetting_fronts))
    ponded_depth_temp = *ponded_depth_cm - f_p*time_step - free_drainage_demand*0;
  else
    ponded_depth_temp = *ponded_depth_cm - f_p*time_step - free_drainage_demand*0;
    
  ponded_depth_temp = fmax(ponded_depth_temp, 0.0);
  
  //printf("f_p value = %lf %lf \n", f_p*10, *ponded_depth_cm*10);
  
  //printf("f_p valueA = %lf %lf \n", f_p, *ponded_depth_cm);
  // checkpoint # AJ
  
  double fp_cm = f_p * time_step + free_drainage_demand/time_step; // infiltration in cm
  
  double precip_mass_to_add;

  if (hp_cm_max>0.0 ) {
    //    printf("got here = %lf %lf %lf \n", ponded_depth_temp, hp_cm_max, h_p);
    if (ponded_depth_temp < hp_cm_max) {
      runoff = 0.0;
      *volin_this_timestep = fmin(*ponded_depth_cm, fp_cm);
      *ponded_depth_cm = *ponded_depth_cm - *volin_this_timestep;
      return runoff;
    }
    else if (ponded_depth_temp > hp_cm_max ) {
      runoff = ponded_depth_temp - hp_cm_max;
      //precip_mass_to_add = fp_cm;
      *ponded_depth_cm = hp_cm_max;
      *volin_this_timestep=fp_cm;
      //printf("HEREss %lf %lf %lf \n", ponded_depth_temp, fp_cm, h_p);
      //abort();
      return runoff;
    }


  }

  else {

    // if it got to this point, no ponding is allowed, either infiltrate or runoff
    // order is important here; assign zero to ponded depth once done using it computing volume in and runoff
    *volin_this_timestep = fmin(*ponded_depth_cm, fp_cm); //
    runoff = *ponded_depth_cm < fp_cm ? 0.0 : (*ponded_depth_cm - *volin_this_timestep);
    *ponded_depth_cm = 0.0;
    
    
    /*
  if ( (*ponded_depth_cm) < fp_cm) {
    //precip_mass_to_add = *ponded_depth_cm;
    //*volin_this_timestep= *ponded_depth_cm;
    //*ponded_depth_cm = 0.0;
    runoff = 0.0;
  }
  else {
    //precip_mass_to_add = fp_cm;
    //*volin_this_timestep = fp_cm;
    //runoff = *ponded_depth_cm - precip_mass_to_add;
    runoff = *ponded_depth_cm - *volin_this_timestep;
    
    //*ponded_depth_cm = precip_mass_to_add;
    //*ponded_depth_cm = 0.0; //*ponded_depth_cm - *volin_this_timestep;
  }
*/
  }
  //printf("Runoff = %lf %lf %lf %lf \n", runoff*10, *ponded_depth_cm*10, precip_mass_to_add*10, f_p*10);

  return runoff;
}

extern void lgar_create_surfacial_front(double *ponded_depth_cm, double *volin, double dry_depth, double theta1,
                                        int *soil_type, struct soil_properties_ *soil_properties, 
                                        double *cum_layer_thickness_cm, int nint,double time_step_s)
{
  // This subroutine is called iff there is no surfacial front, it creates a new front and inserts ponded depth
  // into the soil.  Note ponded_depth_cm is a pointer.   Access it's value as (*ponded_depth_cm).
  //printf("Calculating supercial front now... \n");
  // local vars
  double theta_e,Se,theta_r;
  double delta_theta,hp_cm;
  double vg_alpha_per_cm, vg_m, vg_n, Ksat_cm_per_s;
  double h_min_cm;
  double Geff;
  bool debug_flag=FALSE;
  bool to_bottom = FALSE;
  struct wetting_front *current;
  int layer_num,soil_num,front_num;
  char word[100];
  double guess,amount,diff,Dd,tau,Ktry;
  
  current = head;
  
  amount= (*ponded_depth_cm);
  
  layer_num = 1;   // we only create new surfacial fronts in the first layer
  soil_num = soil_type[layer_num];
  front_num = 1;   // we are creating a new surfacial front, which by definition must be front #1
  
  theta_e = soil_properties[soil_num].theta_e;  // rhs of the new front, assumes theta_e as per Peter
  theta_r = soil_properties[soil_num].theta_r;
  delta_theta =  theta_e - theta1;
  h_min_cm = soil_properties[soil_num].h_min_cm;
  
  
  if(debug_flag) printf("dry_depth= %lf\n",dry_depth);
  //printf("New theta: %lf %lf %lf %lf \n", delta_theta, theta_e, theta1, dry_depth);
  double theta_new;

  if(dry_depth * delta_theta > (*ponded_depth_cm))  // all the ponded depth enters the soil
    {
      
      *volin = *ponded_depth_cm;
      //double time_step = 1.0;
      //theta_new = fmin(theta1 + (*ponded_depth_cm) * (time_step_s/3600.) /dry_depth, theta_e);
      //theta_new = fmin(theta1 + (*ponded_depth_cm) * time_step_s /dry_depth, theta_e);
      theta_new = fmin(theta1 + (*ponded_depth_cm) /dry_depth, theta_e);
      //printf("InFilA: %lf %lf %lf %lf, %lf \n", theta_new, theta1, (*ponded_depth_cm)/dry_depth, dry_depth, (*ponded_depth_cm));
      //listInsertFirst((*ponded_depth_cm)/delta_theta, theta_e, front_num, layer_num, to_bottom);
      
      listInsertFirst(dry_depth, theta_new, front_num, layer_num, to_bottom); // AJ
      //printf ("printing surficial list \n");
      //listPrint();
      *ponded_depth_cm = 0.0;
      hp_cm =0.0;
    }
  else  // not all ponded depth fits in
    {
      //printf("test--- not all water fits in \n");
      //abort();
      *volin = dry_depth * delta_theta;
      *ponded_depth_cm -= dry_depth * delta_theta;
      theta_new = theta_e; //fmin(theta1 + (*ponded_depth_cm) /dry_depth, theta_e);
      // printf("InFilB: %lf %lf %lf %lf %lf %d \n", theta_new, dry_depth, (*ponded_depth_cm), (*ponded_depth_cm)/delta_theta, time_step_s, to_bottom);
      if (dry_depth < cum_layer_thickness_cm[1])
	listInsertFirst(dry_depth, theta_e, front_num, layer_num, to_bottom);
      else
	listInsertFirst(dry_depth, theta_e, front_num, layer_num, 1);
      hp_cm = *ponded_depth_cm;
    } 
  
  current = head;  // must ddo this again because listInsertFirst() created a new *head
  /* // AJ
    Se = 1.0;  // assumed saturated
    current->psi_cm    = 0.0; // assumed saturated
    
  */
  vg_alpha_per_cm    = soil_properties[soil_num].vg_alpha_per_cm;
  vg_m               = soil_properties[soil_num].vg_m;
  vg_n               = soil_properties[soil_num].vg_n;
  Ksat_cm_per_s      = soil_properties[soil_num].Ksat_cm_per_s;

  Se = calc_Se_from_theta(theta_new,theta_e,theta_r);
  current->psi_cm = calc_h_from_Se(Se, vg_alpha_per_cm , vg_m, vg_n);

  //current->K_cm_per_s = soil_properties[soil_num].Ksat_cm_per_s;
  current->K_cm_per_s = calc_K_from_Se(Se, Ksat_cm_per_s, vg_m); // AJ - K_temp in python version for 1st layer

  // AJ
  //Geff = calc_Geff(theta1, theta_e, theta_e, theta_r, vg_alpha_per_cm, vg_n, vg_m, h_min_cm, Ksat_cm_per_s, nint);
  Geff = calc_Geff(theta_new, theta1, theta_e, theta_r, vg_alpha_per_cm, vg_n, vg_m, h_min_cm, Ksat_cm_per_s, nint);
  h_min_cm = 0.0;
  //printf("surficial G+: %lf %lf %lf \n",Geff, theta_new, theta1);
  //current->dzdt_cm_per_s = 1.0/delta_theta*(Ksat_cm_per_s*(Geff+hp_cm)/current->depth_cm+current->K_cm_per_s); AJ

  current->dzdt_cm_per_s = 0.0; //1.0/fabs(theta_new -theta1)*(Ksat_cm_per_s*(Geff+hp_cm)/current->depth_cm+current->K_cm_per_s);
  
  //printf("Surficial: %lf %lf \n",current->dzdt_cm_per_s*3600*10,Ksat_cm_per_s);
  // here I tried a manual scheme to estimate how much theta and depth increase jointly, and conserve mass instead
  // of assuming that the soil saturates.  Didn't work....
  // while (TRUE)
  //    {
  //    printf("Enter guess between %8.6lf and %8.6lf: ", theta1, theta_e);
  //    fgets(word, sizeof(word), stdin);
  //    guess=atof(word);
  //    Se=calc_Se_from_theta(guess,theta_e,theta_r);
  //    Ktry=calc_K_from_Se(Se,Ksat_cm_per_s,vg_m);
  //    delta_theta = guess - theta1;
  //    tau = time_step_s * Ktry/delta_theta;
  //    Geff=calc_Geff(theta1, guess, theta_e, theta_r, vg_alpha_per_cm, vg_n, vg_m, h_min_cm, Ksat_cm_per_s, nint);
  //    Dd = 0.5*(tau+sqrt(tau*tau+4.0*tau*Geff));
  //    diff=amount-Dd*delta_theta;
  //    printf(" Geff=%8.5lf Ktry=%7.5e Dd=%8.6lf diff= %lf\n",Geff,Ktry,Dd,diff);
  //    }
  
  
  return;
  
}


extern double lgar_calc_dry_depth(int nint, double time_step_s, int *soil_type, 
                                  struct soil_properties_ *soil_properties, double *cum_layer_thickness_cm,
                                  double *delta_theta)
{
  /*******************************************/
  /* This routine calculates the "dry depth" */
  /* of a newly created wetting front in     */
  /* the top soil layer after a non-rainy    */
  /* period or a big increase in rainrate    */
  /* on an unsaturated first layer.          */
  
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
  int    fred;                       
  
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
  theta1      = current->theta;                 // water content of the first (most surficial) existing wetting front
  theta_e     = soil_properties[soil_num].theta_e;
  theta2 = theta_e;
  
  *delta_theta = theta_e - current->theta;  // return the delta_theta value to the calling function
  
  tau       = time_step_s * Ksat_cm_per_s/(theta_e-current->theta); //3600
  //printf("Tau = %lf %lf %lf %lf \n",time_step_s,Ksat_cm_per_s,theta_e, current->theta);
  //abort();
  if(theta1>theta_e) {fred=1;}
  
  // printf("theta_e = %lf %lf %lf \n", theta1, theta2, theta_e);
  //  Geff      = calc_Geff(theta1, theta2, theta_e, theta_r, vg_alpha_per_cm, vg_n, vg_m, h_min_cm, Ksat_cm_per_s, nint);

  Geff      = calc_Geff(theta1, theta2, theta_e, theta_r, vg_alpha_per_cm, vg_n, vg_m, h_min_cm, Ksat_cm_per_s, nint);
 
 //###!!! note that dry depth originally has a factor of 0.5 in front
 dry_depth = 0.5*(tau + sqrt( tau*tau + 4.0*tau*Geff) ); 

 //printf("DD, Tau, G : %lf %lf %lf \n", dry_depth, tau, Geff);

 dry_depth = fmin(cum_layer_thickness_cm[layer_num], dry_depth);
 
 return(dry_depth);

 // iff dry depth greater than layer 1 thickness, set dry depth to layer 1 thickness.  Unlikely to occur.
// if(dry_depth > cum_layer_thickness_cm[1])  dry_depth=cum_layer_thickness_cm[1];
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
  //printf("MB layer = %d \n", layer);
  if(current->next != NULL)                  // this is not the last entry in the list
    {
    next=current->next;
    if(next->layer_num == current->layer_num)
      {
	//sum+=(current->depth_cm - base_depth) * fabs(current->theta - next->theta);
      sum+=(current->depth_cm - base_depth) * (current->theta - next->theta);
      //if ((current->theta - next->theta) < 0) abort();
      //printf("in mass balanceA = %lf %lf %lf %lf\n",current->depth_cm, base_depth,current->theta , next->theta);
      }
    else
      {
      sum+=(current->depth_cm - base_depth) * current->theta;
      //printf("in mass balanceB = %lf %lf %lf \n",current->depth_cm, base_depth,current->theta);
      }
    }
  else // this is the last entry in the list.  This must be the deepest front in the final layer
    {
    layer=current->layer_num;
    sum+=current->theta * (current->depth_cm - base_depth);
    //printf("in mass balanceC = %lf %lf %lf\n",current->depth_cm, base_depth,current->theta);
    }
  
  //if(current->next != NULL)
  current=current->next;
  
  } while(current != NULL );   // putting conditional at end of do looop makes sure it executes once.
return sum;
}

extern void lgar_read_vG_param_file(char *vG_param_file_name, int num_soil_types, double wilting_point_psi_cm,
                                    struct soil_properties_ *soil_properties)
{
//###################################################################
// READ THE SOIL PARAMS ******************************************
// OPEN FILE TO READ IN THE vG parameters ffor standard soil types
//###################################################################

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
  //  soil_properties[soil].Ksat_cm_per_s   = Ksat_cm_per_h/3600.0;  // convert from cm/h to cm/s
  soil_properties[soil].Ksat_cm_per_s   = Ksat_cm_per_h;  // convert from cm/h to cm/s
  
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
  //struct wetting_front* current_new;
  //struct wetting_front* next_new;
  //struct wetting_front* previous_new = NULL;

  
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
  bool   fred=FALSE;
  
  
  if(head == NULL)
    {
      return nfronts_analyzed;  // No wetting fronts!!
    }
  
  //current_new = head;  // head is global, and ALWAYS points to the first (most superficial) structure in the linked list
  
  // make sure to use previous state values as current state is updated during the timestep (that's how it is done is Peter's python version)
  
  //current = state_previous;
  current = head;
  // printf("checking.... %lf %lf \n", current->depth_cm, current->theta);
  do  // loop through the wetting fronts
    {
      dzdt = 0.0;
      nfronts_analyzed++;
      
      // copy structure elements into shorter variables names to increase readability
      
      // WETTING FRONT PROPERTIES
      front_num    = current->front_num;    // the front number
      theta        = current->theta;        // water content of this front
      layer_num    = current->layer_num;    // what layer the front is in
      K_cm_per_s   = current->K_cm_per_s;   // K(theta)
      
      //printf("V1: %d %lf %lf \n", current->layer_num, current->K_cm_per_s*3600*10,K_cm_per_s*3600*10);
      if (K_cm_per_s <=0){
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
      //printf("A_dz: %lf %lf \n", theta2, theta1);
      //        max_depth_vec = params.max_depth_vec   //FLO NOT SURE WHAT THIS IS.
      // comments from Peter's python code:
      //#currently, LGAR works with exactly 3 layers (which can be set to have the same properties). Derivatives have
      //#different forms for wetting fronts in different layers, as moisture in above layers must be considered for wetting
      //#fronts in deeper layers. Therefore dZdt calculations look similar between layers, albeit there are more terms in 
	  //#parts of dZdt for deeper layers. Eventually this will probably be replaced with code that generally calculates dZdt
	  //#based on what layer the wetting front is in.
	  
	   bottom_sum = 0.0;  // needed ffor multi-layered dz/dt equation.  Equal to sum from n=1 to N-1 of (L_n/K_n(theta_n))
      //printf("calculating DzDt - to_bottom: %d %d %lf \n",current->front_num, current->to_bottom, current->depth_cm);
      if(current->to_bottom == TRUE) // checkpoint # AJ
	//if(current->to_bottom == FALSE) // checkpoint # AJ
	{
	  //printf("AA: to_bottom: %d %d %lf \n",current->front_num, current->to_bottom,current->dzdt_cm_per_s);
	  //listPrint();
	  if(layer_num > 1) {
	    //bottom_sum += (cum_layer_thickness_cm[layer_num]-cum_layer_thickness_cm[layer_num-1])/K_cm_per_s;
	    current->dzdt_cm_per_s = 0.0;
	  }
	  else
	    current->dzdt_cm_per_s = 0.0;
	  //listPrint();
	  //printf("to_bottom: %d %d %lf \n",current->front_num, current->to_bottom,current->dzdt_cm_per_s);
	  previous = current;
	  current = current->next;  // point to the next link
	  //printf("to_bottomXC: %d %d %lf \n",current->front_num, current->to_bottom,current->dzdt_cm_per_s);
	  continue;                 // go to next front, this one fully penetrates the layer
	}
      else {
	
	if(layer_num > 1) {
	  double K_cm_per_s_prev = previous->K_cm_per_s; 
	  //bottom_sum += (cum_layer_thickness_cm[layer_num]-cum_layer_thickness_cm[layer_num-1])/K_cm_per_s;
	  bottom_sum += (current->depth_cm-cum_layer_thickness_cm[layer_num-1])/K_cm_per_s;
	  //printf("bottom sum = %lf %lf %lf \n", (current->depth_cm-cum_layer_thickness_cm[layer_num-1])*10, K_cm_per_s*3600*10, bottom_sum);
	}
	//printf("D1 (bottom sum) = %d %lf %lf %lf %lf %lf \n",current->front_num, cum_layer_thickness_cm[layer_num], cum_layer_thickness_cm[layer_num-1], (cum_layer_thickness_cm[layer_num]-cum_layer_thickness_cm[layer_num-1]), K_cm_per_s*3600*10);
	// printf("D2 = %lf %lf %lf %lf \n", current->depth_cm, cum_layer_thickness_cm[layer_num-1], K_cm_per_s*3600*10);
	//printf("D3 = %lf \n",bottom_sum/3600); 
      }
      
      if(theta1 > theta2)
	{
	  fred=TRUE;  // this should never happen
	}
      //printf("A: %lf %lf %lf %lf \n", theta1, theta2, K_cm_per_s, K_cm_per_s * 3600*10);
      Geff = calc_Geff(theta1, theta2, theta_e, theta_r, vg_alpha_per_cm, vg_n, vg_m, h_min_cm, Ksat_cm_per_s, nint);
      //printf("G_temp = %lf \n", Geff);
      
      delta_theta = current->theta - next->theta;
      //printf("delta_theta = %lf \n", delta_theta);
      if(current->layer_num == 1)  // this front is in the upper layer
	{
	  
	  dzdt = 1.0/delta_theta*(Ksat_cm_per_s*(Geff+h_p)/current->depth_cm+current->K_cm_per_s);
	  //printf("DzDT1: %lf \n", delta_theta);
	  //printf("DzDT2: %lf %lf %lf %lf \n", Ksat_cm_per_s*3600*10, Geff*10, h_p,current->depth_cm*10);
	  //printf("DzDT3: %lf %lf \n", dzdt*10, current->K_cm_per_s*10);
	}
      else  // we are in the second or greater layer
	{
	  //printf("DzDt layer number = %d \n", layer_num);
	  //double theta_prevA = 0.0;
	  double denominatorA = bottom_sum;
	  //printf("theta loop = %lf %lf \n", theta, current->psi_cm*10);
	  for (int k = 1; k < layer_num; k++) {
	    int soil_numA = soil_type[layer_num-k];
	    double theta_prevA = calc_theta_from_h(current->psi_cm, soil_properties[soil_numA].vg_alpha_per_cm, soil_properties[soil_numA].vg_m, soil_properties[soil_numA].vg_n,soil_properties[soil_numA].theta_e,soil_properties[soil_numA].theta_r);

	    
	    double Se_prevA = calc_Se_from_theta(theta_prevA,soil_properties[soil_numA].theta_e,soil_properties[soil_numA].theta_r);

	    double K_cm_per_s_prevA = calc_K_from_Se(Se_prevA,soil_properties[soil_numA].Ksat_cm_per_s, soil_properties[soil_numA].vg_m);
	    //printf("m = %lf %lf \n", soil_properties[soil_numA].vg_m, soil_properties[soil_numA].Ksat_cm_per_s*3600*10);
	    denominatorA += (cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k-1])/ K_cm_per_s_prevA;
	    //printf("X1 = %lf %lf \n",  (cum_layer_thickness_cm[k] - cum_layer_thickness_cm[k-1]), K_cm_per_s_prevA*10);
	  }

	  
	  double theta_prev = calc_theta_from_h(current->psi_cm, soil_properties[soil_num-1].vg_alpha_per_cm, soil_properties[soil_num-1].vg_m, soil_properties[soil_num-1].vg_n,soil_properties[soil_num-1].theta_e,soil_properties[soil_num-1].theta_r);
	  //printf("psi current = %lf %lf \n", current->psi_cm, theta_prev);
	  double Se_prev = calc_Se_from_theta(theta_prev,soil_properties[soil_num-1].theta_e,soil_properties[soil_num-1].theta_r);


	  //current->psi_cm = calc_h_from_Se(Se, vg_alpha_per_cm , vg_m, vg_n);
	  //current->K_cm_per_s = soil_properties[soil_num].Ksat_cm_per_s;
	   double K_cm_per_s_prev = calc_K_from_Se(Se_prev,soil_properties[soil_num-1].Ksat_cm_per_s, soil_properties[soil_num-1].vg_m);
  
	  //double K_cm_per_s_prev =     K_cm_per_s = calc_K_from_Se(Se, Ks_cm_per_s, vg_m);
	  //previous->K_cm_per_s;   // K(theta_previous)
	  //printf("K_temp_above = %lf %lf %lf \n", K_cm_per_s_prev * 3600*10, current->K_cm_per_s*3600*10, previous->depth_cm);
	  //abort();
	  //## f_temp = [(1/abs(theta-theta_below))*(Z+Z_above+(G_temp+h_p)*K_s/K_temp)*
	  //##                                                      (1/(max_depth_vec[0]/K_temp_above+Z/K_temp))#,
	  //## the hydraulic conductivity in the layer above corresponding to the current wetting front is necessary 
	  //## to calculate dZdt.
	      // checkpoint #2 (in python there is G_temp + h_p
	      Z_cm = depth_cm - cum_layer_thickness_cm[layer_num -1];
	  //dzdt = 1.0/delta_theta *(Z_cm + Geff * Ksat_cm_per_s / K_cm_per_s)/
	  //                        ((Z_cm-cum_layer_thickness_cm[layer_num-1])/K_cm_per_s + bottom_sum);
	  
	   double numerator = depth_cm + (Geff +h_p)* Ksat_cm_per_s / K_cm_per_s;
	  double denominator = cum_layer_thickness_cm[layer_num -1] / K_cm_per_s_prev + bottom_sum;
	  //dzdt = 1.0/delta_theta * numerator / denominator;

	 dzdt = 1.0/delta_theta * numerator / denominatorA;
	 //printf ("layer num: %d %lf %lf \n ", current->layer_num, depth_cm, cum_layer_thickness_cm[layer_num -1]);
	 //printf ("K_s and K_temp = %lf %lf \n", Ksat_cm_per_s  * 3600*10, K_cm_per_s * 3600*10);
	 
	 //printf("N/D: %lf %lf \n", numerator, 1/(denominator/3600.)*10);
	 //printf("D1 = %lf %lf %lf \n", cum_layer_thickness_cm[layer_num -1], K_cm_per_s_prev*3600*10,  cum_layer_thickness_cm[layer_num -1] / K_cm_per_s_prev /3600.);
	 //printf("D2 = %lf \n", bottom_sum);
	 //printf("Full D/N = %lf %lf \n", delta_theta,  numerator/denominatorA*10);
	 //printf("DzDt2 : %lf \n", dzdt*10);
	 //abort();
       }

     previous = current;
     current->dzdt_cm_per_s = dzdt;
     // bottom_sum = 0.0;

     // copy to the current state as well

     //current_new->dzdt_cm_per_s = dzdt;
     
     //printf("DICH: %lf %lf %lf \n", dzdt*3600*10, current->dzdt_cm_per_s*3600*10, current_new->dzdt_cm_per_s*3600*10);
     current = current->next;  // point to the next link
     //current_new = current_new->next;
     //printf("Aborting in DzDt............ \n");
     //abort();
     
   } while(current != NULL );   // putting conditional at end of do looop makes sure it executes at least once

      
 return nfronts_analyzed;
}


extern double lgar_theta_mass_balance(int layer_num, int soil_num, double psi_cm, double new_mass, double prior_mass, double theta_below, double theta_below_above, double wf_depth_cm, double layer_thickness_cm, struct soil_properties_ *soil_properties)
{

  double x = psi_cm;
  double delta_mass = fabs(new_mass - prior_mass);
  double tolerance = 1e-12;
  //printf ("in mass balance: %lf %lf \n", new_mass, prior_mass);

  double factor = 1.0;
  bool switched = false;

  double theta_of_same_wf_in_layer_above,  new_mass_in_layer_above, new_mass_in_same_layer;
  double theta = 0; // this will be returned
  
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

    theta_of_same_wf_in_layer_above = calc_theta_from_h(x, soil_properties[soil_num-1].vg_alpha_per_cm, soil_properties[soil_num-1].vg_m, soil_properties[soil_num-1].vg_n,soil_properties[soil_num-1].theta_e,soil_properties[soil_num-1].theta_r);
    
    theta = calc_theta_from_h(x, soil_properties[soil_num].vg_alpha_per_cm, soil_properties[soil_num].vg_m, soil_properties[soil_num].vg_n,soil_properties[soil_num].theta_e,soil_properties[soil_num].theta_r);
    
    new_mass_in_layer_above = layer_thickness_cm *(theta_of_same_wf_in_layer_above - theta_below_above);
    new_mass_in_same_layer = wf_depth_cm * (theta-theta_below);
    new_mass = new_mass_in_layer_above + new_mass_in_same_layer;

    delta_mass = fabs(new_mass -prior_mass);
  
  }

  //printf ("******************in mass balance: %lf %lf %lf \n", new_mass, prior_mass, theta);
  return theta;
  
}

/*
extern double lgar_theta_mass_balance_3(double current_mass, double old_mass, double *cum_layer_thickness_cm, wetting_front *wf_free_drainage)
{
  double mass_balance_error = fabs(current_mass - old_mass); // mass error

  double factor = 1.0;
  bool switched = false;
  double tolerance = 1e-10;
 
  // check if the difference is less than the tolerance
  if (mass_balance_error <= tolerance) {
    return current_mass;
  }

  while (mass_balance_error > tolerance) {
   
    if (current_mass>old_mass) {
      x = x + 0.01 * factor;
      switched = false;
    }
    else {
      if (!switched) {
	switched = true;
	factor = factor * 0.01;
      }
      x = x - 0.01 * factor;

    }

    wf_free_drainage->depth_cm += 0.01;
    
    current_mass = lgar_calc_mass_bal(0,cum_layer_thickness_cm);
    mass_balance_error =current_mass - (old_mass + precip_mass_to_add) + (actual_ET_demand+free_drainage_demand);
    
   current_mass = current_mass - old_mass


    new_mass = mass_layers;
    
    delta_mass = fabs(new_mass - prior_mass);
    
  }
  
}
*/
extern double lgar_theta_mass_balance_2(int layer_num, int soil_num, double psi_cm, double new_mass, double prior_mass, double depth_cm_old, double *delta_theta, double *delta_thickness, int *soil_type, struct soil_properties_ *soil_properties)
{

  double x = psi_cm;
  double delta_mass = fabs(new_mass - prior_mass);
  double tolerance = 1e-12;
  //printf ("In mass balance (before): %lf %lf \n", new_mass, prior_mass);


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

  //printf("_2 = %d %lf %lf %0.6e\n", layer_num,delta_thickness[layer_num], delta_thickness[1], delta_mass);
  //abort();
  //printf("in mass balance = %d \n", soil_num);
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
    //printf("factor = %.6e \n", factor);
    //x = fmax(x,0.0);
    double theta_layer;
    double mass_layers= 0.0;
    
    theta = calc_theta_from_h(x, soil_properties[soil_num].vg_alpha_per_cm, soil_properties[soil_num].vg_m, soil_properties[soil_num].vg_n,soil_properties[soil_num].theta_e,soil_properties[soil_num].theta_r);
    mass_layers += delta_thickness[layer_num] * (theta - delta_theta[layer_num]);
    //printf("theta1 = %lf %lf %lf %lf %lf \n", theta, delta_theta[layer_num], delta_thickness[layer_num], (theta - delta_theta[layer_num]), mass_layers );
    
    for (int k=1; k<layer_num; k++) {
      int soil_numA =  soil_type[k];
      //theta_layer = calc_theta_from_h(x, soil_properties[soil_num-k].vg_alpha_per_cm, soil_properties[soil_num-k].vg_m, soil_properties[soil_num-k].vg_n,soil_properties[soil_num-k].theta_e,soil_properties[soil_num-k].theta_r);

      theta_layer = calc_theta_from_h(x, soil_properties[soil_numA].vg_alpha_per_cm, soil_properties[soil_numA].vg_m, soil_properties[soil_numA].vg_n,soil_properties[soil_numA].theta_e,soil_properties[soil_numA].theta_r);

      mass_layers += delta_thickness[k] * (theta_layer - delta_theta[k]);
      //printf("thetaA = %d %lf \n", k, delta_thickness[k] * (theta - delta_theta[k]));
      //printf("V0 = %lf %lf %lf \n", theta,x,soil_properties[soil_num].vg_m);
    }


    new_mass = mass_layers;
    
    delta_mass = fabs(new_mass - prior_mass);
    
  }
 
  //printf ("**** In mass balance (after) : %6.15f %6.15f %0.6e \n", new_mass, prior_mass, delta_mass);
  return theta;
  
}
 

  /*
extern double lgar_theta_mass_balance_2(int layer_num, int soil_num, double psi_cm, double new_mass, double prior_mass, double depth_cm_old, double *delta_theta, double *cum_layer_thickness_cm, struct soil_properties_ *soil_properties)
{

  double x = psi_cm;
  double delta_mass = fabs(new_mass - prior_mass);
  double tolerance = 1e-10;
  printf ("In mass balance (before): %lf %lf \n", new_mass, prior_mass);


  double factor = 1.0;
  bool switched = false;

  double theta_of_wf_in_layer1_above, theta_of_wf_in_layer2_above;
  double new_mass_in_layer1_above, new_mass_in_layer2_above, new_mass_in_same_layer;
  double theta = 0; // this will be updated and returned


  // check if the difference is less than the tolerance
  if (delta_mass <= tolerance) {
    theta = calc_theta_from_h(x, soil_properties[soil_num].vg_alpha_per_cm, soil_properties[soil_num].vg_m, soil_properties[soil_num].vg_n,soil_properties[soil_num].theta_e,soil_properties[soil_num].theta_r);
    return theta;
  }

    
  double x1 = (cum_layer_thickness_cm[layer_num-1] - cum_layer_thickness_cm[layer_num-2]);
  double x2 = (cum_layer_thickness_cm[layer_num-2] - cum_layer_thickness_cm[layer_num-3]);
  double x3 = (depth_cm_old - cum_layer_thickness_cm[layer_num-1]);
  
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

    theta_of_wf_in_layer1_above = calc_theta_from_h(x, soil_properties[soil_num-1].vg_alpha_per_cm, soil_properties[soil_num-1].vg_m, soil_properties[soil_num-1].vg_n,soil_properties[soil_num-1].theta_e,soil_properties[soil_num-1].theta_r);

    theta_of_wf_in_layer2_above = calc_theta_from_h(x, soil_properties[soil_num-2].vg_alpha_per_cm, soil_properties[soil_num-2].vg_m, soil_properties[soil_num-2].vg_n,soil_properties[soil_num-2].theta_e,soil_properties[soil_num-2].theta_r);
    
    theta = calc_theta_from_h(x, soil_properties[soil_num].vg_alpha_per_cm, soil_properties[soil_num].vg_m, soil_properties[soil_num].vg_n,soil_properties[soil_num].theta_e,soil_properties[soil_num].theta_r);
    
    new_mass_in_layer1_above = x1 * (theta_of_wf_in_layer1_above - delta_theta[0]);
    new_mass_in_layer2_above = x2 * (theta_of_wf_in_layer2_above - delta_theta[1]);
    new_mass_in_same_layer = x3 * (theta-delta_theta[2]);
    new_mass = new_mass_in_layer1_above + new_mass_in_layer2_above + new_mass_in_same_layer;
    //printf("theta = %lf %lf \n", theta,x);
    delta_mass = fabs(new_mass - prior_mass);
  
  }

  printf ("****************** In mass balance (after) : %lf %lf %lf \n", new_mass, prior_mass, theta);
  return theta;
  
}
*/




      /*************************************************************************************/
/*


      if ( ((l<last_wetting_front_index) && (layer_number_above!=layer_num) && layer_number_below == layer_num) && false) {
	printf("case 4: layer_num_above != layer_num (interface above) \n");
	
	abort();
	//printf("moving layer number above = %d \n", layer_num);
	//listPrint();
	printf("** Depth update (before) = %d %lf %lf \n", layer_num, current->depth_cm, current->dzdt_cm_per_s*3600*10);
	current->depth_cm +=  current->dzdt_cm_per_s * time_step_s;
	printf("** Depth update (after) = %lf %lf \n", current->depth_cm, current->dzdt_cm_per_s*3600*10);
	if (layer_num == 2) {
	  
	  double psi_cm_old = current_old->psi_cm;
	  double psi_cm_below_old = next_old->psi_cm;
	  //printf ("checking1 = %lf %lf \n", theta_old,theta_below_old);

	  //printf ("checking2 = %lf %lf %lf %lf \n", psi_cm_old*10, psi_cm_below_old*10, current_old->psi_cm*10, next_old->psi_cm*10);
	  
	  double theta_of_same_wf_in_layer_above_previous = calc_theta_from_h(psi_cm_old, soil_properties[soil_num-1].vg_alpha_per_cm, soil_properties[soil_num-1].vg_m, soil_properties[soil_num-1].vg_n,soil_properties[soil_num-1].theta_e,soil_properties[soil_num-1].theta_r); 
	  
	  double theta_of_wf_below_in_layer_above_previous = calc_theta_from_h(psi_cm_below_old, soil_properties[soil_num-1].vg_alpha_per_cm, soil_properties[soil_num-1].vg_m, soil_properties[soil_num-1].vg_n,soil_properties[soil_num-1].theta_e,soil_properties[soil_num-1].theta_r); 
	  
	  //printf ("checking3 = %lf %lf \n", theta_of_same_wf_in_layer_above_previous,theta_of_wf_below_in_layer_above_previous);
	  
	  
	  double prior_mass_in_layer_above = cum_layer_thickness_cm[layer_num-1] * (theta_of_same_wf_in_layer_above_previous-theta_of_wf_below_in_layer_above_previous);
	    
	  double prior_mass_in_layer = (current_old->depth_cm - cum_layer_thickness_cm[layer_num-1]) * (current_old->theta-next_old->theta);
	    
	  //printf("checking4 = %lf %lf %lf \n", cum_layer_thickness_cm[layer_num-1], current_old->depth_cm, (current_old->depth_cm - cum_layer_thickness_cm[layer_num-1]));
	    
	  //printf ("checking5 = %lf %lf \n", prior_mass_in_layer_above*10, prior_mass_in_layer*10);
	  
	  // checkpoint # : AJ
	  double free_drainage_demand = 0;
	  double actual_ET_demand = 0;
	  
	  double prior_mass = prior_mass_in_layer_above + prior_mass_in_layer - (free_drainage_demand+actual_ET_demand);
	  //*(1 if l==wf_that_supplies_free_drainage_demand else 0) + precip_mass_to_add*(1 if l==wf_that_supplies_free_drainage_demand else 0)
	  //printf("prior mass = %lf \n", prior_mass*10);
	  
	  
	  double psi_cm = current->psi_cm;
	  double psi_cm_below = next->psi_cm;
	  
	  double theta_of_wf_in_layer_above = calc_theta_from_h(psi_cm, soil_properties[soil_num-1].vg_alpha_per_cm, soil_properties[soil_num-1].vg_m, soil_properties[soil_num-1].vg_n,soil_properties[soil_num-1].theta_e,soil_properties[soil_num-1].theta_r);
	  
	  double theta_of_wf_below_in_layer_above = calc_theta_from_h(psi_cm_below, soil_properties[soil_num-1].vg_alpha_per_cm, soil_properties[soil_num-1].vg_m, soil_properties[soil_num-1].vg_n,soil_properties[soil_num-1].theta_e,soil_properties[soil_num-1].theta_r); 
	  
	  //printf("checking4: = %lf %lf \n", psi_cm*10, psi_cm_below*10);
	  
	  //printf("checking5: = %lf %lf \n", theta_of_wf_in_layer_above, theta_of_wf_below_in_layer_above);
	    
	  double wf_depth_cm = (current->depth_cm - cum_layer_thickness_cm[layer_num-1]);
	  double new_mass_in_layer_above = cum_layer_thickness_cm[layer_num-1] * (theta_of_wf_in_layer_above - theta_of_wf_below_in_layer_above);
	  double new_mass_in_same_layer = wf_depth_cm * (current->theta-next->theta);
	  double new_mass = new_mass_in_layer_above + new_mass_in_same_layer;
	  
	  //printf("checking6: = %lf %lf %lf %lf %lf \n", wf_depth_cm,cum_layer_thickness_cm[layer_num], 10*current->depth_cm, current->theta, next->theta);
	  //abort();
	  //printf("checking7: = %lf  %lf \n", new_mass_in_layer_above*10, new_mass_in_same_layer*10);
	  
	  printf("new/prior mass = %lf %lf \n", new_mass*10, prior_mass*10);
	  
	  new_mass = lgar_theta_mass_balance(layer_num, soil_num, psi_cm, new_mass, prior_mass, next->theta, theta_of_wf_below_in_layer_above, wf_depth_cm, cum_layer_thickness_cm[1], soil_properties);
	  
	  current->theta = fmin(new_mass, theta_e);
	  printf("theta after balance = %lf \n", new_mass);
	}
	
	if (layer_num == 3) {

	    double psi_cm_old = current_old->psi_cm;
	    double psi_cm_below_old = next_old->psi_cm;
	    
	  
	    double theta_of_same_wf_in_layer1_above_previous = calc_theta_from_h(psi_cm_old, soil_properties[soil_num-1].vg_alpha_per_cm, soil_properties[soil_num-1].vg_m, soil_properties[soil_num-1].vg_n,soil_properties[soil_num-1].theta_e,soil_properties[soil_num-1].theta_r); 
	  
	    double theta_of_wf_below_in_layer1_above_previous = calc_theta_from_h(psi_cm_below_old, soil_properties[soil_num-1].vg_alpha_per_cm, soil_properties[soil_num-1].vg_m, soil_properties[soil_num-1].vg_n,soil_properties[soil_num-1].theta_e,soil_properties[soil_num-1].theta_r); 

	    double theta_of_same_wf_in_layer2_above_previous = calc_theta_from_h(psi_cm_old, soil_properties[soil_num-2].vg_alpha_per_cm, soil_properties[soil_num-2].vg_m, soil_properties[soil_num-2].vg_n,soil_properties[soil_num-2].theta_e,soil_properties[soil_num-2].theta_r); 
	    
	    double theta_of_wf_below_in_layer2_above_previous = calc_theta_from_h(psi_cm_below_old, soil_properties[soil_num-2].vg_alpha_per_cm, soil_properties[soil_num-2].vg_m, soil_properties[soil_num-2].vg_n,soil_properties[soil_num-2].theta_e,soil_properties[soil_num-2].theta_r);
	    
	    double x1 = (cum_layer_thickness_cm[layer_num-1] - cum_layer_thickness_cm[layer_num-2]);
	    double x2 = (cum_layer_thickness_cm[layer_num-2] - cum_layer_thickness_cm[layer_num-3]);
	    double x3 = (current->depth_cm - cum_layer_thickness_cm[layer_num-1]);
	    
	    //printf ("checking3 = %lf %lf \n", theta_of_same_wf_in_layer_above_previous,theta_of_wf_below_in_layer_above_previous);
	    double prior_mass_in_layer1_above = (cum_layer_thickness_cm[layer_num-1] - cum_layer_thickness_cm[layer_num-2]) * (theta_of_same_wf_in_layer1_above_previous-theta_of_wf_below_in_layer1_above_previous);
	    
	    double prior_mass_in_layer2_above = (cum_layer_thickness_cm[layer_num-2] - cum_layer_thickness_cm[layer_num-3]) * (theta_of_same_wf_in_layer2_above_previous-theta_of_wf_below_in_layer2_above_previous);
	    //printf("checking4 = %lf %lf %lf \n", cum_layer_thickness_cm[layer_num], (cum_layer_thickness_cm[layer_num-2] - cum_layer_thickness_cm[layer_num-3]),cum_layer_thickness_cm[layer_num-3]);
	    
	    
	    double prior_mass_in_layer = (current_old->depth_cm - cum_layer_thickness_cm[layer_num-1]) * (current_old->theta-next_old->theta);
	    
	    
	    //printf("checking4A = %lf %lf %lf %lf \n", x2 , prior_mass_in_layer2_above*10, theta_of_same_wf_in_layer2_above_previous, theta_of_wf_below_in_layer2_above_previous);
	    //abort();
	    //printf("checking4B = %lf %lf  \n", x1 , prior_mass_in_layer1_above*10);
	    //printf("checking4C = %lf %lf  %lf %lf \n", x3 , prior_mass_in_layer*10,cum_layer_thickness_cm[layer_num-1], current_old->depth_cm);
	    // checkpoint # : AJ
	    double free_drainage_demand = 0;
	    double actual_ET_demand = 0;
	    
	    double prior_mass = prior_mass_in_layer1_above +  prior_mass_in_layer2_above + prior_mass_in_layer - (free_drainage_demand+actual_ET_demand);
	    //*(1 if l==wf_that_supplies_free_drainage_demand else 0) + precip_mass_to_add*(1 if l==wf_that_supplies_free_drainage_demand else 0)
	    //printf("checking3 = %lf %lf %lf %lf \n", prior_mass_in_layer2_above*10, prior_mass_in_layer1_above*10, prior_mass_in_layer*10, prior_mass*10);
	    //abort();
	    //printf("prior mass = %lf \n", prior_mass*10);
	    

	    double psi_cm = current->psi_cm;
	    double psi_cm_below = next->psi_cm;
	    printf("checking4: = %lf %lf \n", psi_cm*10, psi_cm_below*10);
	    

	    double theta_of_wf_in_layer1_above = calc_theta_from_h(psi_cm, soil_properties[soil_num-1].vg_alpha_per_cm, soil_properties[soil_num-1].vg_m, soil_properties[soil_num-1].vg_n,soil_properties[soil_num-1].theta_e,soil_properties[soil_num-1].theta_r);

	    double theta_of_wf_in_layer2_above = calc_theta_from_h(psi_cm, soil_properties[soil_num-2].vg_alpha_per_cm, soil_properties[soil_num-2].vg_m, soil_properties[soil_num-2].vg_n,soil_properties[soil_num-2].theta_e,soil_properties[soil_num-2].theta_r);
	    
	    double theta_of_wf_below_in_layer1_above = calc_theta_from_h(psi_cm_below, soil_properties[soil_num-1].vg_alpha_per_cm, soil_properties[soil_num-1].vg_m, soil_properties[soil_num-1].vg_n,soil_properties[soil_num-1].theta_e,soil_properties[soil_num-1].theta_r);
	    
	    double theta_of_wf_below_in_layer2_above = calc_theta_from_h(psi_cm_below, soil_properties[soil_num-2].vg_alpha_per_cm, soil_properties[soil_num-2].vg_m, soil_properties[soil_num-2].vg_n,soil_properties[soil_num-2].theta_e,soil_properties[soil_num-2].theta_r); 

	    //printf("checking4: = %lf %lf \n", psi_cm*10, psi_cm_below*10);
	    
	    //printf("checking5: = %lf %lf \n", theta_of_wf_in_layer_above, theta_of_wf_below_in_layer_above);
	    
	    // these are set like this because there is no wetting front below the current one
	    double theta_below = next->theta;
	    //	double theta_of_wf_below_in_layer1_above = 0;
	    //double theta_of_wf_below_in_layer2_above = 0;
	    
	    double new_mass_in_layer2_above = x1 * (theta_of_wf_in_layer2_above - theta_of_wf_below_in_layer2_above);
	    double new_mass_in_layer1_above = x2 * (theta_of_wf_in_layer1_above - theta_of_wf_below_in_layer1_above);
	    double new_mass_in_same_layer = x3 * (theta-theta_below);
	    double new_mass = new_mass_in_layer2_above + new_mass_in_layer1_above + new_mass_in_same_layer;
	    printf("checking5 = %lf %lf %lf \n",new_mass_in_layer2_above*10,new_mass_in_layer1_above*10, new_mass_in_same_layer*10);
	    printf("checking6+ = %lf %lf %lf \n", x3, theta, theta_below);
	    printf("proir/new mass = %lf %lf \n", prior_mass*10, new_mass*10);
	    
	    double delta_theta[] = {theta_of_wf_below_in_layer1_above, theta_of_wf_below_in_layer2_above, theta_below};
	    
	    new_mass = lgar_theta_mass_balance_2(layer_num, soil_num, psi_cm, new_mass, prior_mass, current->depth_cm, delta_theta, cum_layer_thickness_cm, soil_properties);
	    
	    current->theta = fmin(new_mass, theta_e);
	    printf("new mass after balance = %lf \n", new_mass);
	    
	    
	}
	//printf("Quiting....................B \n");
	// abort();
      }



 */
/*
create plane
edit mode --> ctrl + R (to add line) or use loop cut

edge loop : alt + select edge

alt + ctrl to select all faces

*/
