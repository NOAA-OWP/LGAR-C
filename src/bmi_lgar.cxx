#ifndef BMI_LGAR_CXX_INCLUDED
#define BMI_LGAR_CXX_INCLUDED


#include <stdio.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include "../bmi/bmi.hxx"
#include "../include/bmi_lgar.hxx"
#include "../include/all.h"
//#include "../include/soil_moisture_profile.hxx"


void BmiLGAR::
Initialize (std::string config_file)
{
  if (config_file.compare("") != 0 ) {
    this->model = new lgar_model;
    //soil_moisture_profile::LGAR(config_file, model);
    lgar_initialize(config_file, model);
  }
}


void BmiLGAR::
Update()
{
  //  model->LGARUpdate();
  LGARUpdate(model);
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
  if (name.compare("soil_storage_model") == 0)   // int
    return 0;
  else if (name.compare("soil_storage") == 0 || name.compare("soil_storage_change") == 0 || name.compare("soil_water_table") == 0) // double
    return 1; 
  else if (name.compare("soil_moisture_profile") == 0 || name.compare("soil_moisture_layered") == 0) // array of doubles 
    return 2; 
  else 
    return -1;
}


std::string BmiLGAR::
GetVarType(std::string name)
{
  if (name.compare("soil_storage_model") == 0)
    return "int";
  else if (name.compare("soil_storage") == 0 || name.compare("soil_storage_change") == 0 || name.compare("soil_water_table") == 0)
    return "double";
  else if (name.compare("soil_moisture_profile") == 0 || name.compare("soil_moisture_layered") == 0)
    return "double";
  else
    return "";
}


int BmiLGAR::
GetVarItemsize(std::string name)
{
  if (name.compare("soil_storage_model") == 0)
    return sizeof(int);
  else if (name.compare("soil_storage") == 0 || name.compare("soil_storage_change") == 0 || name.compare("soil_water_table") == 0)
    return sizeof(double);
  else if (name.compare("soil_moisture_profile") == 0 || name.compare("soil_moisture_layered") == 0)
    return sizeof(double);
  else
    return 0;
}


std::string BmiLGAR::
GetVarUnits(std::string name)
{
  if (name.compare("soil_storage") == 0 || name.compare("soil_storage_change") == 0  || name.compare("soil_water_table") == 0)
    return "m";
  else if (name.compare("soil_moisture_profile") == 0 || name.compare("soil_moisture_layered") == 0)
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
  if (name.compare("soil_storage") == 0 || name.compare("soil_storage_change") == 0 || name.compare("soil_water_table") == 0)
    return "node";
  else if (name.compare("soil_moisture_profile") == 0 || name.compare("soil_moisture_layered") == 0)
    return "node";
  else
    return "none";
}


void BmiLGAR::
GetGridShape(const int grid, int *shape)
{
  if (grid == 2) {
    shape[0] = this->model->shape[0];
  }
}


void BmiLGAR::
GetGridSpacing (const int grid, double * spacing)
{
  if (grid == 0) {
    spacing[0] = this->model->spacing[0];
  }
}


void BmiLGAR::
GetGridOrigin (const int grid, double *origin)
{
  if (grid == 0) {
    origin[0] = this->model->origin[0];
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
    return this->model->shape[0];
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
  throw coupler::NotImplemented();
}


void BmiLGAR::
GetGridY(const int grid, double *y)
{
  throw coupler::NotImplemented();
}


void BmiLGAR::
GetGridZ(const int grid, double *z)
{
  throw coupler::NotImplemented();
}


int BmiLGAR::
GetGridNodeCount(const int grid)
{
  throw coupler::NotImplemented();
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
  throw coupler::NotImplemented();
}


int BmiLGAR::
GetGridFaceCount(const int grid)
{
  throw coupler::NotImplemented();
}


void BmiLGAR::
GetGridEdgeNodes(const int grid, int *edge_nodes)
{
  throw coupler::NotImplemented();
}


void BmiLGAR::
GetGridFaceEdges(const int grid, int *face_edges)
{
  throw coupler::NotImplemented();
}


void BmiLGAR::
GetGridFaceNodes(const int grid, int *face_nodes)
{
  throw coupler::NotImplemented();
}


void BmiLGAR::
GetGridNodesPerFace(const int grid, int *nodes_per_face)
{
  throw coupler::NotImplemented();
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
  if (name.compare("soil_storage") == 0)
    return (void*)(&this->model->soil_storage);
  else if (name.compare("soil_storage_change") == 0)
    return (void*)(&this->model->soil_storage_change_per_timestep);
  else  if (name.compare("soil_water_table") == 0)
    return (void*)(&this->model->water_table_thickness);
  else if (name.compare("soil_moisture_profile") == 0)
    return (void*)this->model->soil_moisture_profile;
  else if (name.compare("soil_moisture_layered") == 0)
    return (void*)this->model->soil_moisture_layered;
  else if (name.compare("soil_storage_model") == 0)
    return (void*)(&this->model->soil_storage_model);
  else {
    std::stringstream errMsg;
    errMsg << "variable "<< name << " does not exist";
    throw std::runtime_error(errMsg.str());
    return NULL;
  }
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
  return "LGAR BMI";
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
