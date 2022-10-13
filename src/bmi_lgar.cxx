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


string verbosity="none";

void BmiLGAR::
Initialize (std::string config_file)
{
  if (config_file.compare("") != 0 ) {
    this->model = new lgar_model_;
    lgar_initialize(config_file, model);
  }
  
  num_giuh_ordinates = model->lgar_bmi_params.num_giuh_ordinates;

  giuh_ordinates = new double[num_giuh_ordinates];
  giuh_runoff_queue = new double[num_giuh_ordinates+1];
  
  for (int i=0; i<num_giuh_ordinates;i++)
    giuh_ordinates[i] = model->lgar_bmi_params.giuh_ordinates[i+1]; // note lgar use 1 indexing

  for (int i=0; i<=num_giuh_ordinates;i++)
    giuh_runoff_queue[i] = 0.0;

}

void BmiLGAR::
Update()
{
  if (verbosity.compare("none") != 0) {
    std::cout<<"*** LASAM BMI Update... ***  \n";
  }
  //   lgar_update(this->model); 
  //listPrint();
  double mm_to_cm = 0.1;

  // local variables for readibility
  int subcycles = model->lgar_bmi_params.forcing_interval;
  int num_layers = model->lgar_bmi_params.num_layers;
  
  // full timestep (timestep of the forcings)
  double precip_timestep_cm = 0.0;
  double PET_timestep_cm = 0.0;
  double ponded_depth_cm = 0.0;
  double AET_timestep_cm = 0.0;
  double volstart_timestep_cm = 0.0;
  double volend_timestep_cm = lgar_calc_mass_bal(num_layers,model->lgar_bmi_params.cum_layer_thickness_cm); //0.0; // this should not be reset to 0.0 in the for loop
  double volin_timestep_cm = 0.0;
  double volon_timestep_cm = 0.0;
  double volrunoff_timestep_cm = 0.0;
  double volrech_timestep_cm = 0.0;
  double surface_runoff_timestep_cm = 0.0; // direct surface runoff
  double volrunoff_giuh_timestep_cm = 0.0;
  double volQ_timestep_cm = 0.0;
  
  // subtimestep (timestep of the model)
  double precip_subtimestep_cm;
  double precip_subtimestep_cm_per_h;
  double PET_subtimestep_cm;
  double ponded_depth_subtimestep_cm;
  double AET_subtimestep_cm;
  double volstart_subtimestep_cm;
  double volend_subtimestep_cm = volend_timestep_cm; //0.0; // this should not be reset to 0.0 in the for loop
  double volin_subtimestep_cm;
  double volon_subtimestep_cm;
  double volrunoff_subtimestep_cm;
  double volrech_subtimestep_cm;
  double surface_runoff_subtimestep_cm; // direct surface runoff
  double precip_previous_subtimestep_cm;
  double volrunoff_giuh_subtimestep_cm;
  double volQ_subtimestep_cm;
  
  double subtimestep_h = model->lgar_bmi_params.timestep_h;
  int nint = model->lgar_bmi_params.nint;
  double wilting_point_psi_cm = model->lgar_bmi_params.wilting_point_psi_cm;
  double AET_thresh_Theta = 0.85;    // scaled soil moisture (0-1) above which AET=PET (fix later!)
  double AET_expon = 1.0;  // exponent that allows curvature of the rising portion of the Budyko curve (fix later!)

  
  for (int cycle=1; cycle <= subcycles; cycle++) {

    if (verbosity.compare("high") == 0 || verbosity.compare("medium") == 0) {
      std::cout<<"*** ----------------- Subcycle ------------------: "<< cycle <<" of "<<subcycles<<"\n";
    }
    
    state_previous = NULL;
    state_previous = listCopy(head);
    
    precip_subtimestep_cm_per_h = model->lgar_bmi_params.precipitation_cm_per_h * mm_to_cm / double(subcycles); // rate; cm/hour
    PET_subtimestep_cm = model->lgar_bmi_params.PET_cm_per_h * mm_to_cm / double(subcycles);
    ponded_depth_subtimestep_cm = precip_subtimestep_cm_per_h * subtimestep_h;

    if (verbosity.compare("high") == 0 || verbosity.compare("medium") == 0) {
      std::cout<<"Pr [cm/h], Pr [cm], subtimestep [h] = "<<model->lgar_bmi_params.precipitation_cm_per_h<<", "<< precip_subtimestep_cm<<", "<< subtimestep_h<<" ("<<subtimestep_h*3600<<" sec)"<<"\n";
    }


    AET_subtimestep_cm = 0.0;
    volstart_subtimestep_cm = 0.0;
    volin_subtimestep_cm = 0.0;
    volon_subtimestep_cm= 0.0;
    volrunoff_subtimestep_cm = 0.0;
    volrech_subtimestep_cm = 0.0;
    surface_runoff_subtimestep_cm = 0.0;

    precip_previous_subtimestep_cm = model->lgar_bmi_params.precip_previous_timestep_cm;
    
    num_layers = model->lgar_bmi_params.num_layers;
    double delta_theta;   // the width of a front, such that its volume=depth*delta_theta
    double dry_depth;
    
    
    if (PET_subtimestep_cm > 0.0) {
      // Calculate AET from PET and root zone soil moisture.  Note PET was reduced if raining
      
      AET_subtimestep_cm = calc_aet(PET_subtimestep_cm, subtimestep_h, wilting_point_psi_cm, model->soil_properties, model->lgar_bmi_params.layer_soil_type, AET_thresh_Theta, AET_expon);
    }
    

    precip_subtimestep_cm = precip_subtimestep_cm_per_h * subtimestep_h;
    precip_timestep_cm += precip_subtimestep_cm;
    PET_timestep_cm += fmax(PET_subtimestep_cm,0.0); // ensures non-negative PET
    volstart_subtimestep_cm = lgar_calc_mass_bal(num_layers,model->lgar_bmi_params.cum_layer_thickness_cm);

    
    int soil_num = model->lgar_bmi_params.layer_soil_type[head->layer_num];
    double theta_e = model->soil_properties[soil_num].theta_e;
    bool is_top_wf_saturated = head->theta >= theta_e ? true : false;
    bool create_surficial_front = (precip_previous_subtimestep_cm == 0.0 && precip_subtimestep_cm >0.0);
    
    double mass_source_to_soil_timestep = 0.0;
    
    int wf_free_drainage_demand = wetting_front_free_drainage();

    if (verbosity.compare("high") == 0 || verbosity.compare("medium") == 0) {
      std::string flag = (create_surficial_front && !is_top_wf_saturated) == true ? "Yes" : "No";
      std::cout<<"Superficial wetting front created? "<< flag << "\n";
    }
    
    //if the follow is true, that would mean there is no wetting front in the top layer to accept the water, must create one.
    if(create_surficial_front && !is_top_wf_saturated)  {
      
      double temp_pd = 0.0; // necessary to assign zero precip due to the creation of new wetting front; AET will still be taken out of the layers
      
      lgar_move_wetting_fronts(&temp_pd, subtimestep_h, wf_free_drainage_demand, volend_subtimestep_cm, num_layers, &AET_subtimestep_cm, model->lgar_bmi_params.cum_layer_thickness_cm, model->lgar_bmi_params.layer_soil_type, model->soil_properties);
      
      dry_depth = lgar_calc_dry_depth(nint, subtimestep_h, model->lgar_bmi_params.layer_soil_type, model->soil_properties, model->lgar_bmi_params.cum_layer_thickness_cm,&delta_theta);
      
      double theta1 = head->theta;
      lgar_create_surfacial_front(&ponded_depth_subtimestep_cm, &volin_subtimestep_cm, dry_depth, theta1, model->lgar_bmi_params.layer_soil_type, model->soil_properties, model->lgar_bmi_params.cum_layer_thickness_cm, nint, subtimestep_h);
      
      state_previous = NULL;
      state_previous = listCopy(head);
      
      volin_timestep_cm += volin_subtimestep_cm;

      if (verbosity.compare("high") == 0) {
	std::cout<<"New wetting front created...\n";
	std::cout<<" "<<"\n";
	listPrint();
      }
    }


    if (ponded_depth_subtimestep_cm > 0 && !create_surficial_front) {
      
      volrunoff_subtimestep_cm = lgar_insert_water(&ponded_depth_subtimestep_cm, &volin_subtimestep_cm, precip_subtimestep_cm_per_h, dry_depth, nint, subtimestep_h, wf_free_drainage_demand, model->lgar_bmi_params.layer_soil_type, model->soil_properties, model->lgar_bmi_params.cum_layer_thickness_cm);

      volin_timestep_cm += volin_subtimestep_cm;
      volrunoff_timestep_cm += volrunoff_subtimestep_cm;
      volrech_subtimestep_cm = volin_subtimestep_cm; // this gets updated later, probably not needed here
      volon_timestep_cm += ponded_depth_subtimestep_cm;
      
      //printf("Mass in = %lf %lf %lf \n", volin_subtimestep_cm, volrech_subtimestep_cm, volrunoff_subtimestep_cm);
      if (volrunoff_subtimestep_cm < 0) abort();  
    }
    else {
      //printf("wetting front created = %lf %d \n", ponded_depth_cm ,!create_surficial_front );
      double hp_cm_max = 0.0;
      
      if (ponded_depth_subtimestep_cm < hp_cm_max) {
	volrunoff_timestep_cm += 0.0;
	volon_timestep_cm = ponded_depth_subtimestep_cm;
	ponded_depth_subtimestep_cm = 0.0;
	volrunoff_subtimestep_cm = 0.0;
      }
      else {
	volrunoff_subtimestep_cm = (ponded_depth_subtimestep_cm - hp_cm_max);
	volrunoff_timestep_cm += (ponded_depth_subtimestep_cm - hp_cm_max);
	volon_timestep_cm = hp_cm_max;
	ponded_depth_subtimestep_cm = hp_cm_max;
      }
    }
    

    if (!create_surficial_front) {
      lgar_move_wetting_fronts(&volin_subtimestep_cm, subtimestep_h, wf_free_drainage_demand, volend_subtimestep_cm, num_layers, &AET_subtimestep_cm, model->lgar_bmi_params.cum_layer_thickness_cm, model->lgar_bmi_params.layer_soil_type, model->soil_properties);
      
      // this is the volume of water leaving through the bottom
      volrech_subtimestep_cm = volin_subtimestep_cm;
      volrech_timestep_cm += volrech_subtimestep_cm;
    }
    
    
    int num_dzdt_calculated = lgar_dzdt_calc(nint, model->lgar_bmi_params.layer_soil_type, model->soil_properties, model->lgar_bmi_params.cum_layer_thickness_cm, ponded_depth_subtimestep_cm);

    AET_timestep_cm += AET_subtimestep_cm;
    volrech_timestep_cm += volrech_subtimestep_cm;
    
    volend_subtimestep_cm = lgar_calc_mass_bal(num_layers,model->lgar_bmi_params.cum_layer_thickness_cm);
    volend_timestep_cm = volend_subtimestep_cm;
    model->lgar_bmi_params.precip_previous_timestep_cm = precip_subtimestep_cm;
    
    double local_mb = volstart_subtimestep_cm + precip_subtimestep_cm - volrunoff_subtimestep_cm - AET_subtimestep_cm - volon_subtimestep_cm - volrech_subtimestep_cm - volend_subtimestep_cm;

    
    // compute giuh runoff for the sub-timestep
    surface_runoff_subtimestep_cm = volrunoff_subtimestep_cm;
    volrunoff_giuh_subtimestep_cm = giuh_convolution_integral(volrunoff_subtimestep_cm, num_giuh_ordinates, giuh_ordinates, giuh_runoff_queue);

    surface_runoff_timestep_cm += surface_runoff_subtimestep_cm ;
    volrunoff_giuh_timestep_cm += volrunoff_giuh_subtimestep_cm;

    // total mass of water leaving the system, at this time it is the giuh-only, but later will add groundwater component as well.

    volQ_timestep_cm += volrunoff_giuh_subtimestep_cm;

    if (verbosity.compare("high") == 0 || verbosity.compare("medium") == 0) {
      std::cout<<"Printing wetting fronts at this subtimestep... \n";
      listPrint();
    }

    bool unexpected_local_error = fabs(local_mb) > 1.0e-7 ? true : false;
    
    if (verbosity.compare("high") == 0 || verbosity.compare("medium") == 0 || unexpected_local_error) {
      printf("\nLocal mass balance at this timestep... \n\
      Error         = %14.10f \n\
      Initial water = %14.10f \n\
      Water added   = %14.10f \n\
      Infiltration  = %14.10f \n\
      Runoff        = %14.10f \n\
      AET           = %14.10f \n\
      Percolation   = %14.10f \n\
      Final water   = %14.10f \n", local_mb, volstart_subtimestep_cm, precip_subtimestep_cm, volin_timestep_cm, volrunoff_subtimestep_cm, AET_subtimestep_cm, volrech_subtimestep_cm, volend_subtimestep_cm);

      if (unexpected_local_error) {
	printf("Local mass balance (in this timestep) is %14.10f, larger than expected, needs some debugging...\n ",local_mb);
	abort();
      }
	
    }
    
    if (fabs(local_mb) >1e-7) {
      printf("\nLocal mass balance at this timestep... \n\
      Error         = %14.10f \n\
      Initial water = %14.10f \n\
      Water added   = %14.10f \n\
      Infiltration  = %14.10f \n\
      Runoff        = %14.10f \n\
      AET           = %14.10f \n\
      Percolation   = %14.10f \n\
      Final water   = %14.10f \n", local_mb, volstart_subtimestep_cm, precip_subtimestep_cm, volin_timestep_cm, volrunoff_subtimestep_cm, AET_subtimestep_cm, volrech_subtimestep_cm, volend_subtimestep_cm);
      printf("Local mass balance (in this timestep) is %14.10f, larger than expected, needs some debugging...\n ",local_mb);
      abort();
    }
    
    assert (head->depth_cm > 0.0); // check on negative layer depth


  } // end of cycle loop

  //update number of wetting fronts
  model->lgar_bmi_params.num_wetting_fronts = listLength();
  model->lgar_bmi_params.soil_thickness_wetting_fronts = new double[model->lgar_bmi_params.num_wetting_fronts];
  model->lgar_bmi_params.soil_moisture_wetting_fronts = new double[model->lgar_bmi_params.num_wetting_fronts];

  // update thickness/depth and soil moisture of wetting fronts (used for model coupling)
  struct wetting_front *current = head;
  for (int i=0; i<model->lgar_bmi_params.num_wetting_fronts; i++) {
    assert (current != NULL);
    model->lgar_bmi_params.soil_moisture_wetting_fronts[i] = current->theta;
    model->lgar_bmi_params.soil_thickness_wetting_fronts[i] = current->depth_cm * model->units.cm_to_m;
    current = current->next;
  }
  
  // add to mass balance timestep variables
  model->lgar_mass_balance.volprecip_timestep_cm = precip_timestep_cm;
  model->lgar_mass_balance.volin_timestep_cm = volin_timestep_cm;
  model->lgar_mass_balance.volend_timestep_cm = volend_timestep_cm;
  model->lgar_mass_balance.volAET_timestep_cm = AET_timestep_cm;
  model->lgar_mass_balance.volrech_timestep_cm = volrech_timestep_cm;
  model->lgar_mass_balance.volrunoff_timestep_cm = volrunoff_timestep_cm;
  model->lgar_mass_balance.volrunoff_giuh_timestep_cm = volrunoff_giuh_timestep_cm;
  model->lgar_mass_balance.volQ_timestep_cm = volQ_timestep_cm;
  model->lgar_mass_balance.volPET_timestep_cm = PET_timestep_cm;
  
  // add to mass balance accumulated variables
  model->lgar_mass_balance.volprecip_cm += precip_timestep_cm;
  model->lgar_mass_balance.volin_cm += volin_timestep_cm;
  model->lgar_mass_balance.volend_cm = volend_timestep_cm;
  model->lgar_mass_balance.volAET_cm += AET_timestep_cm;
  model->lgar_mass_balance.volrech_cm += volrech_timestep_cm;
  model->lgar_mass_balance.volrunoff_cm += volrunoff_timestep_cm;
  model->lgar_mass_balance.volrunoff_giuh_cm += volrunoff_giuh_timestep_cm;
  model->lgar_mass_balance.volQ_cm += volQ_timestep_cm;//model->lgar_mass_balance.volQ_timestep_cm;
  model->lgar_mass_balance.volPET_cm += PET_timestep_cm;


  // unit conversion
   // add to mass balance timestep variables
  bmi_unit_conv.volprecip_timestep_m = precip_timestep_cm * model->units.cm_to_m;
  bmi_unit_conv.volin_timestep_m = volin_timestep_cm * model->units.cm_to_m;
  bmi_unit_conv.volend_timestep_m = volend_timestep_cm * model->units.cm_to_m;
  bmi_unit_conv.volAET_timestep_m = AET_timestep_cm * model->units.cm_to_m;
  bmi_unit_conv.volrech_timestep_m = volrech_timestep_cm * model->units.cm_to_m;
  bmi_unit_conv.volrunoff_timestep_m = volrunoff_timestep_cm * model->units.cm_to_m;
  bmi_unit_conv.volrunoff_giuh_timestep_m = volrunoff_giuh_timestep_cm * model->units.cm_to_m;
  bmi_unit_conv.volQ_timestep_m = volQ_timestep_cm * model->units.cm_to_m;
  bmi_unit_conv.volPET_timestep_m = PET_timestep_cm * model->units.cm_to_m;
    
}

/*
void BmiLGAR::
Update()
{
  lgar_update(this->model);
  
  //double surface_runoff_m = this->model->lgar_mass_balance.volrunoff_timestep_cm;

  //double giuh_runoff_cm = 10 * giuh_convolution_integral(surface_runoff_m,num_giuh_ordinates,
  //						 giuh_ordinates,giuh_runoff_queue);
  //std::cout<<"giuh runoff "<<giuh_runoff_cm<<" "<<surface_runoff_m<<"\n";

  //this->model->lgar_mass_balance = this->model->lgar_mass_balance.volrunoff_giuh_timestep_cm;
  //massbal_struct->vol_out_giuh+=giuh_runoff_m;
  //abort();
}

*/

struct lgar_model_* BmiLGAR::get_model()
{
  return model;
}

void BmiLGAR::
global_mass_balance()
{
  lgar_global_mass_balance(this->model, giuh_runoff_queue);
}

void BmiLGAR::
UpdateUntil(double t)
{
  //model->LGARUpdate();
}


void BmiLGAR::
Finalize()
{
  // if (this->model)
  //  this->model->~LGAR();
}


int BmiLGAR::
GetVarGrid(std::string name)
{
  if (name.compare("soil_storage_model") == 0)   // int
    return 0;
  else if (name.compare("precipitation_rate") == 0 || name.compare("precipitation") == 0)
    return 1;
  else if (name.compare("potential_evapotranspiration_rate") == 0 || name.compare("potential_evapotranspiration") == 0 || name.compare("actual_evapotranspiration") == 0) // double
    return 1;
  else if (name.compare("surface_runoff") == 0 || name.compare("giuh_runoff") == 0 || name.compare("soil_storage") == 0) // double
    return 1;
  else if (name.compare("total_discharge") == 0 || name.compare("infiltration") == 0 || name.compare("percolation") == 0) // double
    return 1;
  else if (name.compare("soil_moisture_layer") == 0 || name.compare("soil_thickness_layer") == 0) // array of doubles (fixed length)
    return 2;
  else if (name.compare("soil_moisture_wetting_fronts") == 0 || name.compare("soil_thickness_wetting_fronts") == 0) // array of doubles (dynamic length)
    return 3; 
  else 
    return -1;
}


std::string BmiLGAR::
GetVarType(std::string name)
{
  int var_grid = GetVarGrid(name);
  
  if (var_grid == 0)
    return "int";
  else if (var_grid == 1 || var_grid == 2 || var_grid == 3) 
    return "double";
  else
    return "";
}


int BmiLGAR::
GetVarItemsize(std::string name)
{
  int var_grid = GetVarGrid(name);

   if (var_grid == 0)
    return sizeof(int);
  else if (var_grid == 1 || var_grid == 2 || var_grid == 3)
    return sizeof(double);
  else
    return 0;
}


std::string BmiLGAR::
GetVarUnits(std::string name)
{
  if (name.compare("precipitation_rate") == 0 || name.compare("potential_evapotranspiration_rate") == 0)
    return "mm h^-1";
  else if (name.compare("precipitation") == 0 || name.compare("potential_evapotranspiration") == 0 || name.compare("actual_evapotranspiration") == 0) // double
    return "m";
  else if (name.compare("surface_runoff") == 0 || name.compare("giuh_runoff") == 0 || name.compare("soil_storage") == 0) // double
    return "m";
   else if (name.compare("total_discharge") == 0 || name.compare("infiltration") == 0 || name.compare("percolation") == 0) // double
    return "m";
  else if (name.compare("soil_moisture_layer") == 0 || name.compare("soil_moisture_wetting_fronts") == 0) // array of doubles 
    return "none";
  else if (name.compare("soil_thickness_layer") == 0 || name.compare("soil_thickness_wetting_fronts") == 0) // array of doubles 
    return "m"; 
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
  else if (name.compare("surface_runoff") == 0 || name.compare("giuh_runoff") == 0 || name.compare("soil_storage") == 0) // double
    return "node";
   else if (name.compare("total_discharge") == 0 || name.compare("infiltration") == 0 || name.compare("percolation") == 0) // double
    return "node";
  else if (name.compare("soil_moisture_layer") == 0 || name.compare("soil_moisture_wetting_fronts") == 0) // array of doubles 
    return "node";
  else if (name.compare("soil_thickness_layer") == 0 || name.compare("soil_thickness_wetting_fronts") == 0) // array of doubles 
    return "node"; 
  else 
    return "none";
}


void BmiLGAR::
GetGridShape(const int grid, int *shape)
{
  if (grid == 2) {
    shape[0] = this->model->lgar_bmi_params.shape[0];
  }
}


void BmiLGAR::
GetGridSpacing (const int grid, double * spacing)
{
  if (grid == 0) {
    spacing[0] = this->model->lgar_bmi_params.spacing[0];
  }
}


void BmiLGAR::
GetGridOrigin (const int grid, double *origin)
{
  if (grid == 0) {
    origin[0] = this->model->lgar_bmi_params.origin[0];
  }
}


int BmiLGAR::
GetGridRank(const int grid)
{
  if (grid == 0 || grid == 1 || grid == 2 || grid == 3)
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
    return this->model->lgar_bmi_params.num_layers;
  else if (grid == 3) // number of wetting fronts (dynamic)
    return this->model->lgar_bmi_params.num_wetting_fronts;
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
 
  if (name.compare("precipitation_rate") == 0) {
    return (void*)(&this->model->lgar_bmi_params.precipitation_cm_per_h);
  }
  else if (name.compare("precipitation") == 0) {
    //return (void*)(&this->model->lgar_mass_balance.volprecip_timestep_cm);
    return (void*)(&bmi_unit_conv.volprecip_timestep_m);
  }
  else  if (name.compare("potential_evapotranspiration_rate") == 0) {
    return (void*)(&this->model->lgar_bmi_params.PET_cm_per_h);
  }
  else  if (name.compare("potential_evapotranspiration") == 0) {
    // return (void*)(&this->model->lgar_bmi_params.PET_cm);
    return (void*)(&bmi_unit_conv.volPET_timestep_m);
  }
  else  if (name.compare("actual_evapotranspiration") == 0) {
    // return (void*)(&this->model->lgar_bmi_params.AET_cm);
    return (void*)(&bmi_unit_conv.volAET_timestep_m);
  }
  else  if (name.compare("surface_runoff") == 0) {
    //return (void*)(&this->model->lgar_mass_balance.volrunoff_timestep_cm);
    return (void*)(&bmi_unit_conv.volrunoff_timestep_m);
  }
  else  if (name.compare("giuh_runoff") == 0) {
    // return (void*)(&this->model->lgar_mass_balance.volrunoff_giuh_timestep_cm);
    return (void*)(&bmi_unit_conv.volrunoff_giuh_timestep_m);
  }
  else  if (name.compare("soil_storage") == 0) {
    // return (void*)(&this->model->lgar_mass_balance.volrunoff_giuh_timestep_cm);
    return (void*)(&bmi_unit_conv.volend_timestep_m);
  }
  else  if (name.compare("total_discharge") == 0) {
    // return (void*)(&this->model->lgar_mass_balance.volQ_timestep_cm);
    return (void*)(&bmi_unit_conv.volQ_timestep_m);
  }
  else  if (name.compare("infiltration") == 0) {
    // return (void*)(&this->model->lgar_mass_balance.volin_timestep_cm);
    return (void*)(&bmi_unit_conv.volin_timestep_m);
  }
  else  if (name.compare("percolation") == 0) {
    // return (void*)(&this->model->lgar_mass_balance.volrech_timestep_cm);
    return (void*)(&bmi_unit_conv.volrech_timestep_m);
  }
  else if (name.compare("soil_moisture_layer") == 0)
    return (void*)this->model->lgar_bmi_params.soil_moisture_layer;
  else if (name.compare("soil_thickness_layer") == 0)
    return (void*)this->model->lgar_bmi_params.soil_moisture_layer;
  else if (name.compare("soil_moisture_wetting_fronts") == 0)
    return (void*)this->model->lgar_bmi_params.soil_moisture_wetting_fronts;
  else if (name.compare("soil_thickness_wetting_fronts") == 0)
    return (void*)this->model->lgar_bmi_params.soil_thickness_wetting_fronts;
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
  return "LASAM (lumped arid/semi-arid model";
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
  return 0.0;
}


double BmiLGAR::
GetCurrentTime () {
  return 0.0;
}


std::string BmiLGAR::
GetTimeUnits() {
  return "s";
}


double BmiLGAR::
GetTimeStep () {
  return 0;
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
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridY(const int grid, double *y)
{
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridZ(const int grid, double *z)
{
  throw bmi_lgar::NotImplemented();
}


int BmiLGAR::
GetGridNodeCount(const int grid)
{
  throw bmi_lgar::NotImplemented();
  /*
  if (grid == 0)
    return this->model->shape[0];
  else
    return -1;
  */
}


int BmiLGAR::
GetGridEdgeCount(const int grid)
{
  throw bmi_lgar::NotImplemented();
}


int BmiLGAR::
GetGridFaceCount(const int grid)
{
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridEdgeNodes(const int grid, int *edge_nodes)
{
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridFaceEdges(const int grid, int *face_edges)
{
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridFaceNodes(const int grid, int *face_nodes)
{
  throw bmi_lgar::NotImplemented();
}


void BmiLGAR::
GetGridNodesPerFace(const int grid, int *nodes_per_face)
{
  throw bmi_lgar::NotImplemented();
}

#endif
