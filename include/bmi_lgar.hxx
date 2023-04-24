#ifndef BMI_LGAR_HXX_INCLUDED
#define BMI_LGAR_HXX_INCLUDED

/*
  authors : Ahmad Jan and Fred Ogden
  year    : 2022
  email   : ahmad.jan@noaa.gov
  Description: BMI for the Lumped Aric/semi-arid model (LASAM)
*/

using namespace std;

#include <string.h>
#include "../bmi/bmi.hxx"
#include "all.hxx"
#include <stdexcept>

extern "C" {
#include "../giuh/giuh.h"
}

namespace bmi_lgar {
class NotImplemented : public std::logic_error {
  public:
  NotImplemented() : std::logic_error("Not Implemented Function in LGAR") { };
};

}

class BmiLGAR : public bmixx::Bmi {
public:
  BmiLGAR() {
    this->input_var_names[0] = "precipitation_rate";
    this->input_var_names[1] = "potential_evapotranspiration_rate";
    this->input_var_names[2] = "soil_temperature_profile";
    
    this->output_var_names[0] = "soil_moisture_layers";
    this->output_var_names[1] = "soil_moisture_wetting_fronts";
    this->output_var_names[2] = "soil_thickness_layers";
    this->output_var_names[3] = "soil_thickness_wetting_fronts";
    this->output_var_names[4] = "soil_num_wetting_fronts";
    
    // vis outputs
    this->output_var_names[5]  = "precipitation";
    this->output_var_names[6]  = "potential_evapotranspiration";
    this->output_var_names[7]  = "actual_evapotranspiration";
    this->output_var_names[8]  = "surface_runoff"; // direct surface runoff
    this->output_var_names[9]  = "giuh_runoff";
    this->output_var_names[10] = "soil_storage";
    this->output_var_names[11] = "total_discharge";
    this->output_var_names[12] = "infiltration";
    this->output_var_names[13] = "percolation";
    this->output_var_names[14] = "groundwater_to_stream_recharge";
    this->output_var_names[15] = "mass_balance";
    
    /*
    this->output_var_names[13] = "cum_precipitation";
    this->output_var_names[14] = "cum_potential_evapotranspiration";
    this->output_var_names[15] = "cum_actual_evapotranspiration";
    this->output_var_names[16] = "cum_surface_runoff"; // direct surface runoff
    this->output_var_names[17] = "cum_giuh_runoff";
    this->output_var_names[18] = "cum_soil_storage";
    this->output_var_names[19] = "cum_total_discharge";
    this->output_var_names[20] = "cum_infiltration";
    this->output_var_names[21] = "cum_percolation";
    */
  };
  
  void Initialize(std::string config_file);
  
  void Update();
  void UpdateUntil(double time);
  void Finalize();

  std::string GetComponentName();
  int GetInputItemCount();
  int GetOutputItemCount();
  std::vector<std::string> GetInputVarNames();
  std::vector<std::string> GetOutputVarNames();
  
  int GetVarGrid(std::string name);
  std::string GetVarType(std::string name);
  int GetVarItemsize(std::string name);
  std::string GetVarUnits(std::string name);
  int GetVarNbytes(std::string name);
  std::string GetVarLocation(std::string name);
  
  double GetCurrentTime();
  double GetStartTime();
  double GetEndTime();
  std::string GetTimeUnits();
  double GetTimeStep();
  
  void GetValue(std::string name, void *dest);
  void *GetValuePtr(std::string name);
  void GetValueAtIndices(std::string name, void *dest, int *inds, int count);
  
  void SetValue(std::string name, void *src);
  void SetValueAtIndices(std::string name, int *inds, int len, void *src);
  
  int GetGridRank(const int grid);
  int GetGridSize(const int grid);
  std::string GetGridType(const int grid);
  
  void GetGridShape(const int grid, int *shape);
  void GetGridSpacing(const int grid, double *spacing);
  void GetGridOrigin(const int grid, double *origin);
  
  void GetGridX(const int grid, double *x);
  void GetGridY(const int grid, double *y);
  void GetGridZ(const int grid, double *z);
  
  int GetGridNodeCount(const int grid);
  int GetGridEdgeCount(const int grid);
  int GetGridFaceCount(const int grid);
  
  void GetGridEdgeNodes(const int grid, int *edge_nodes);
  void GetGridFaceEdges(const int grid, int *face_edges);
  void GetGridFaceNodes(const int grid, int *face_nodes);
  void GetGridNodesPerFace(const int grid, int *nodes_per_face);
  void global_mass_balance();
  struct model_state* get_model();
  
private:
  struct model_state* state;
  static const int input_var_name_count = 3;
  static const int output_var_name_count = 16;
  
  std::string input_var_names[input_var_name_count];
  std::string output_var_names[output_var_name_count];
  int num_giuh_ordinates;
  double *giuh_ordinates;
  double *giuh_runoff_queue;

  // unit conversion
  //struct unit_conversion units;
  struct bmi_unit_conversion {
    double volprecip_timestep_m;
    double volin_timestep_m;
    double volend_timestep_m;
    double volAET_timestep_m;
    double volrech_timestep_m;
    double volrunoff_timestep_m;
    double volrunoff_giuh_timestep_m;
    double volQ_timestep_m;
    double volQ_gw_timestep_m;
    double volPET_timestep_m;
    double mass_balance_m;
  };

  struct bmi_unit_conversion bmi_unit_conv;
};


#ifdef NGEN
extern "C"
{
  
  /**
   * Construct this BMI instance as a normal C++ object, to be returned to the framework.
   *
   * @return A pointer to the newly allocated instance.
   */
  
  BmiLGAR *bmi_model_create()
  {
    return new BmiLGAR();
  }
  
  /**
   * @brief Destroy/free an instance created with @see bmi_model_create
   * 
   * @param ptr 
   */
  
  void bmi_model_destroy(BmiLGAR *ptr)
  {
    delete ptr;
  }
  
}

#endif

#endif
