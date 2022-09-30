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
  
  for (int i=0; i<num_giuh_ordinates;i++) {
    giuh_ordinates[i] = model->lgar_bmi_params.giuh_ordinates[i+1]; // note lgar use 1 indexing
    std::cout<<"giuh = "<<num_giuh_ordinates<<" "<<giuh_ordinates[i]<<"\n";
  }

  for (int i=0; i<=num_giuh_ordinates;i++) {
    giuh_runoff_queue[i] = 0.0;
  }

}

void BmiLGAR::
Update()
{
  std::cout<<"LGAT BMI Update... \n";
  //listPrint();
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
  double surface_runoff_m;
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
    surface_runoff_m = 0.0;

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

    
    // compute giuh runoff for the sub-timestep
    surface_runoff_m = this->model->lgar_mass_balance.volrunoff_timestep_cm;
    model->lgar_mass_balance.volrunoff_giuh_timestep_cm = 10 * giuh_convolution_integral(surface_runoff_m,num_giuh_ordinates,
											 giuh_ordinates,giuh_runoff_queue);

    model->lgar_mass_balance.volrunoff_giuh_cm +=  model->lgar_mass_balance.volrunoff_giuh_timestep_cm;

    if(VERBOSE > 1) {
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


  } // end of cycle loop 
  
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

void BmiLGAR::
global_mass_balance()
{
  lgar_global_mass_balance(this->model, giuh_runoff_queue);
}

void BmiLGAR::
UpdateUntil(double t)
{
  //model->LGARUpdate();
  /*
  if (model->soil_storage_model == "conceptual" || model->soil_storage_model == "Conceptual") {
    model->LGARFromConceptualReservoir();
  }
  else if (model->soil_storage_model == "layered" || model->soil_storage_model == "Layered") {
    model->LGARFromLayeredReservoir();
  }
  else {
    std::stringstream errMsg;
    errMsg << "Soil moisture profile OPTION provided in the config file is " << model->soil_storagemodel<< ", which should be either \'concepttual\' or \'layered\' " <<"\n";
    throw std::runtime_error(errMsg.str());
    
    }*/
  //  this->model->LGARVertical();
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
  if (name.compare("soil_storage_model") == 0)   // int (ke
    return 0;
  else if (name.compare("precipitation") == 0 || name.compare("potential_evapotranspiration") == 0 || name.compare("actual_evapotranspiration") == 0) // double
    return 1; 
  else if (name.compare("soil_moisture_layer") == 0 || name.compare("soil_moisture_wetting_front") == 0) // array of doubles 
    return 2; 
  else 
    return -1;
}


std::string BmiLGAR::
GetVarType(std::string name)
{
  if (name.compare("soil_storage_model") == 0)
    return "int";
  else if (name.compare("precipitation") == 0 || name.compare("potential_evapotranspiration") == 0 || name.compare("actual_evapotranspiration") == 0)
    return "double";
  else if (name.compare("soil_moisture_layer") == 0 || name.compare("soil_moisture_wetting_front") == 0)
    return "double";
  else
    return "";
}


int BmiLGAR::
GetVarItemsize(std::string name)
{
  if (name.compare("soil_storage_model") == 0)
    return sizeof(int);
  else if (name.compare("precipitation") == 0 || name.compare("potential_evapotranspiration") == 0 || name.compare("actual_evapotranspiration") == 0)
    return sizeof(double);
  else if (name.compare("soil_moisture_layer") == 0 || name.compare("soil_moisture_wetting_front") == 0)
    return sizeof(double);
  else
    return 0;
}


std::string BmiLGAR::
GetVarUnits(std::string name)
{
  if (name.compare("precipitation") == 0 || name.compare("potential_evapotranspiration") == 0 || name.compare("actual_evapotranspiration") == 0)
    return "m";
  else if (name.compare("soil_moisture_layer") == 0 || name.compare("soil_moisture_wetting_front") == 0)
    return "none";
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
  if (name.compare("precipitation") == 0 || name.compare("potential_evapotranspiration") == 0 || name.compare("actual_evapotranspiration") == 0)
    return "node";
  else if (name.compare("soil_moisture_layer") == 0 || name.compare("soil_moisture_wetting_front") == 0)
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
  if (grid == 0 || grid == 1 || grid == 2)
    return 1;
  else
    return -1;
}


int BmiLGAR::
GetGridSize(const int grid)
{
  if (grid == 0 || grid == 1)
    return 1;
  else if (grid == 2)
    return this->model->lgar_bmi_params.shape[0];
  else
    return -1;
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
 
  if (name.compare("precipitation") == 0)
    return (void*)(&this->model->lgar_bmi_params.precipitation_cm);
  else  if (name.compare("potential_evapotranspiration") == 0)
    return (void*)(&this->model->lgar_bmi_params.PET_cm);
  else  if (name.compare("actual_evapotranspiration") == 0)
    return (void*)(&this->model->lgar_bmi_params.AET_cm);
  else if (name.compare("soil_moisture_layer") == 0)
    return (void*)this->model->lgar_bmi_params.soil_moisture_layer;
  else if (name.compare("soil_moisture_wetting_front") == 0)
    return (void*)this->model->lgar_bmi_params.soil_moisture_wetting_front;
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
  return "LGAR model";
}


int BmiLGAR::
GetInputItemCount()
{
  return this->input_var_name_count;

  /* // this is for dynamically setting input vars
  std::vector<std::string>* names_m = model->InputVarNamesModel();
  int input_var_name_count_m = names_m->size();
  
  assert (this->input_var_name_count >= input_var_name_count_m);
  return input_var_name_count_m;
 */
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

  /* // this is for dynamically setting input vars
  std::vector<std::string>* names_m = model->InputVarNamesModel();
  
  for (int i=0; i<this->input_var_name_count; i++) {
    if (std::find(names_m->begin(), names_m->end(), this->input_var_names[i]) != names_m->end()) {
      names.push_back(this->input_var_names[i]);
    }
  }
  */
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

#endif
