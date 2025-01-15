#ifndef _ALL_HXX

#define _ALL_HXX

/*
  authors : Ahmad Jan and Fred Ogden and Peter La Follette
  year    : 2022
  email   : ahmad.jan@noaa.gov
  - This header file constains functions' definitions used in the lgar.cxx, bmi_lasam.cxx and in other files.
  - Parts of the code were originally written by Fred Ogden and modified/extended by Ahmad Jan for LGAR/LASAM
    implementation and to make it bmi compliant.
  - The struct model_state is needed by bmi_lasam, which encloses many structs needed for advancing a timestep via bmi.
  - LASAM : Lumped Arid/semi-arid Model
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <stddef.h>
#include <assert.h>
#include <stdbool.h>
#include <vector>
#include <time.h>
#include <sstream>

using namespace std;

#define TRUE 1
#define FALSE 0
#define ONE 1

extern string verbosity;

#define use_bmi_flag FALSE       // TODO set to TRUE to run in BMI environment

#define MAX_NUM_SOIL_LAYERS 4
#define MAX_NUM_SOIL_TYPES 25 //changed back to 25 from 15, because the file that loads soil types for Bushland relies on entries 16, 17, and 18 in the .dat file.
//and generally, although we might want to restrict soil types to the 12 recognized ones and a few invalid soil types, in theory what if a user wanted to use a custom .dat file with a large number of soils?
#define MAX_SOIL_NAME_CHARS 25
#define MAX_NUM_WETTING_FRONTS 300


// Define a data structure to hold everything that describes a wetting front
struct wetting_front
{
  double depth_cm;         // depth down from the land surface (absolute depth)
  double theta;            // water content of the soil moisture block
  double psi_cm;           // psi calculated at rhs of the current wetting front
  double K_cm_per_h;       // the value of K(theta) associated with the wetting front
  int    layer_num;        // the layer containing this wetting front.
  int    front_num;        // the wetting front number (might be irrelevant), but useful to debug
  bool   to_bottom;        // TRUE iff this wetting front is in contact with the layer bottom
  double dzdt_cm_per_h;    // use to store the calculated wetting front speed
  bool   is_WF_GW;         // true if WF is in contact with GW, false if wetting front is only in contact with surface 
  struct wetting_front *next;  // pointer to the next wetting front.
};

/* head is a GLOBALLY defined pointer to the first link in the wetting front list.
   Making it a local variable in main() makes all linked list operations
   in subroutines a pain of referencing.  Since it is just one thing,
   making it global just makes everything easier.
*/


// Define a data structure to hold properties and parameters for each soil type
struct soil_properties_  /* note the trailing underscore on the name.  It is just part of the name */
{
  char soil_name[MAX_SOIL_NAME_CHARS];  // string to hold the soil name
  double theta_r;          // residual water content
  double theta_e;          // water content at effective saturation <= porosity
  double vg_alpha_per_cm;  // van Genuchten  "alpha" cm^(-1)
  double vg_n;             // van Genuchten  "n"
  double vg_m;             // van Genuchten  "m"
  double bc_lambda;        // Brooks & Corey pore distribution index
  double bc_psib_cm;       // Brooks & Corey bubbling pressure head (cm)
  double h_min_cm;         // the minimum Geff calculated as per Morel-Seytoux and Khanji
  double Ksat_cm_per_h;    // saturated hydraulic conductivity cm/s
  double theta_wp;         // water content at wilting point [-]
};


// Define a struct for unit conversion
struct unit_conversion
{
  double cm_to_mm = 10;
  double mm_to_cm = 0.1;
  double cm_to_m = 0.01;

  double hr_to_sec = 3600; // hour to seconds
};

// Define a data structure for parameters set by the bmi and not through the config file
// the structure holds non-pointer bmi input variables only; other bmi input variable can be put in lgar_bmi_parameters struct
struct lgar_bmi_input_parameters
{
  double precipitation_mm_per_h; // rainfall precip in mm per hour (input)
  double PET_mm_per_h;           // potential evapotranspiration in mm (input)
};


// Define a data structure for parameters accessed by the bmi and use or pass to lgar modules
struct lgar_bmi_parameters
{
  int    shape[3];
  double spacing[8];
  double origin[3];
  double *layer_thickness_cm;      // 1D array of layer thicknesses in cm, read from config file and static
  int    *layer_soil_type;         // 1D (int) array of layers soil type, read from config file, each integer represent a soil type
  int    num_layers;               // number of actual soil layers
  int    num_wetting_fronts;       // number of wetting fronts
  int    num_cells_temp;           // number of cells of the discretized soil temperature profile
  double *cum_layer_thickness_cm;  // cumulative thickness of layers, allocate memory at run time
  double soil_depth_cm;            // depth of the computational domain (i.e., depth of the last/deepest soil layer from the surface)
  double initial_psi_cm;           // model initial (psi) condition
  double timestep_h;               // model timestep in hours
  double forcing_resolution_h;     // forcing resolution in hours
  double minimum_timestep_h;       // minimum time step in hours, only used if adaptive_timestep is true
  int    forcing_interval;         // = forcing_resolution_h/timestep_h
  int    num_soil_types;           // number of soil types; must be less than or equal to MAX_NUM_SOIL_TYPES
  double AET_cm;                   // actual evapotranspiration in cm

  //double *soil_moisture_layers;    // 1D array of thetas (mean soil moisture content) per layer; output option to other models if needed
  double *soil_moisture_wetting_fronts; /* 1D array of thetas (soil moisture content) per wetting front;
					   output to other models (e.g. soil freeze-thaw) */
  double *soil_depth_wetting_fronts;    /* 1D array of absolute depths of the wetting fronts [meters];
					    output to other models (e.g. soil freeze-thaw) */
  double *soil_temperature;              // 1D array of soil temperature [K]; bmi input for coupling lasam to soil freeze thaw model
  double *soil_temperature_z;            /* 1D array of soil discretization associated with temperature profile [m];
					    depth from the surface in meters */
  double *frozen_factor;                 // frozen factor added to the hydraulic conductivity due to coupling to soil freeze-thaw
  double  wilting_point_psi_cm;          // wilting point (the amount of water not available for plants or not accessible by plants)
  double  field_capacity_psi_cm;          // field capacity represented as a capillary head. Note that both wilting point and field capacity are specified for the whole model domain with single values
  double root_zone_depth_cm;             // maximum depth from which roots extract water
  bool   use_closed_form_G = false;      /* true if closed form of capillary drive calculation is desired, false if numeric integral
					    for capillary drive calculation is desired */
  bool   adaptive_timestep = false;      // if set to true, model uses adaptive timestep. In this case, the minimum timestep is the timestep specified in the config file. The maximum time step will be equal to the forcing resolution
  bool   TO_enabled = false;             // if set to true, model uses multilayer TO model for recharge 
  bool   free_drainage_enabled = false;  // When running in LGAR mode (TO_enabled is false), then free_drainage_enabled will specify whether the lower boundary condition is no flow (false), or free drainage (true). 
  double mbal_tol;                       // if a substep's mass balance error is larger than this number, the model will abort. By default it is set to a large value (10 cm).
  double ponded_depth_cm;                // amount of water on the surface unavailable for surface runoff
  double ponded_depth_max_cm;            // maximum amount of water on the surface unavailable for surface runoff
  double precip_previous_timestep_cm;    // amount of rainfall (previous time step)

  int    nint = 120;            // number of trapezoids used in integrating the Geff function
  double time_s;                // current time [s] (this is the bmi output 'time')
  double endtime_s;             // simulation endtime in seconds (bmi output endtime)
  int    timesteps;             // number of timesteps until the current time 
  int    sft_coupled = 0;       // model coupling flag. if true, lasam is coupled to soil freeze thaw model; default is uncoupled version
  
  double *giuh_ordinates;       // geomorphological instantaneous unit hydrograph
  int    num_giuh_ordinates;    // number of giuh ordinates

  int  calib_params_flag = 0;  // flag for calibratable parameters; if true, then calibratable params are updated otherwise not
  bool is_invalid_soil_type  = true;  // checks if the provided soil type is valid for the model, if not then return Q_out = precip
};

// Define a data structure for local (timestep) and global mass balance parameters
struct lgar_mass_balance_variables
{
  // for local mass balance (compute mass balance at each timestep)
  double volstart_timestep_cm;       // initial volume of water in the soil at each timestep
  double volend_timestep_cm;         // volume of water at the end of timestep
  double volprecip_timestep_cm;      // volume of rainfall at each timestep
  double volin_timestep_cm;          // volume of infiltrated water at each timestep (water that will be added to the soil)
  double volon_timestep_cm;          // volume of water on the surface (ponded water) at each timestep
  double volrunoff_timestep_cm;      // volume of water surface runoff at each timestep
  double volAET_timestep_cm;         // volume of AET at each timestep
  double volPET_timestep_cm;         // volume of PET at each timestep
  double volrech_timestep_cm;        // volume of water leaving soil to the ground water (ground water recharge)
  double volrunoff_giuh_timestep_cm; // volume of giuh runoff at each timestep
  double volQ_timestep_cm;           // total outgoing water (giuh_runoff + volrech)
  double volQ_gw_timestep_cm;        /* volume of water from groundwater to stream
					(i.e., stream water recharge from groundwater reservoir) */
  
  // for global mass balance (compute cumulative mass balance)
  double volstart_cm;         // initial volume of water in the soil (at timestep 0)
  double volend_cm;           // volume of water
  double volprecip_cm;        // volume of rainfall
  double volin_cm;            // volume of infiltrated water
  double volon_cm;            // volume of water on the surface (ponded water)
  double volrunoff_cm;        // volume of water surface runoff
  double volAET_cm;           // volume of AET
  double volPET_cm;           // volume of PET

  double volrech_cm;          // volume of water leaving soil through the bottom of the domain (ground water recharge)
  double volrunoff_giuh_cm;   // volume of giuh runoff
  double volQ_cm;             // total outgoing water
  double volQ_gw_cm;          // outgoing water from ground reservoir to stream channel
  double volchange_calib_cm;  // change in the amount of water due to calibratable parameters
  double local_mass_balance;  // local (per timestep) mass balance error
};

// Define a data structure for calibratable parameters
// the structure holds pointer bmi output variables
struct lgar_calib_parameters
{
  double *theta_e;               // theta_e = smcmax [-]
  double *theta_r;               // theta_r = smcmin [-]
  double *vg_n;                  // Van Genuchten n [-]
  double *vg_alpha;              // Van Genuchten alpha [1/cm]
  double *Ksat;                  // Hydraulic conductivity [cm/hr]
  double field_capacity_psi;    // field capacity in capillary head [cm]
  double ponded_depth_max;      // maximum ponded depth of surface water [cm]

};

// nested structure of structures; main structure for the use in bmi
struct model_state
{
  struct wetting_front*               head           = NULL; // head pointer to the current state
  struct wetting_front*               state_previous = NULL; // head pointer to the previous state,
                                                             // used in computing derivatives and mass balance
  struct soil_properties_*            soil_properties;       // dynamic allocation
  struct lgar_bmi_parameters          lgar_bmi_params;
  struct lgar_mass_balance_variables  lgar_mass_balance;
  struct unit_conversion              units;
  struct lgar_bmi_input_parameters*   lgar_bmi_input_params;
  struct lgar_calib_parameters        lgar_calib_params;
};


/* next, function prototypes. */
/* function prototypes provide the compiler with variable types and order in the calling statement */
/* any time a function is called, it must contain the same number, order, and type of variables    */

/*########################################*/
/*   Linked list code function prototypes */
/*########################################*/
// 1st  entry extern means it lives in a different source file
// 2nd entry is the type of variable it returns (void means that it returns nothing)
// third is the name of the function/subroutine
// inside parentheses are the types of require arguments, names don't matter


extern void                     listPrint(struct wetting_front* head);
extern int                      listLength(struct wetting_front* head);
extern int                      listLength_surface(struct wetting_front* head);
extern int                      listLength_TO_WFs_above_surface_WFs(struct wetting_front* head);
extern int                      listLength_TO_WFs_in_rz(double rzd, struct wetting_front* head);
extern int                      listLength_TO_WFs_in_rz_nonzero_depth(double rzd, struct wetting_front* head);
extern bool                     listIsEmpty();
extern struct wetting_front*    listDeleteFirst(struct wetting_front** head);
extern struct wetting_front*    listFindFront(int i, struct wetting_front* head, struct wetting_front* head_old);
extern struct wetting_front*    listDeleteFront(int front_num, struct wetting_front** head, int *soil_type, struct soil_properties_ *soil_properties);
extern void                     listSortFrontsByDepth(struct wetting_front *head);
extern void                     listSortFrontsByPsi(struct wetting_front *head);
extern void                     listSendToTop(struct wetting_front *head);
extern void                     listInsertFirst(double d, double t, int f, int l, bool b, struct wetting_front** head, bool is_WF_GW);
extern struct wetting_front*    listInsertFront(double d, double t, int f, int l, bool b, struct wetting_front** head, bool is_WF_GW);
extern struct wetting_front*    listInsertFrontAtDepth(int numlay, double *tvec,double d, double t, struct wetting_front* head, bool is_WF_GW);
extern void                     listReverseOrder(struct wetting_front** head_ref);
extern bool                     listFindLayer(struct wetting_front* link, int num_layers, double *cum_layer_thickness_cm,
					      int *lives_in_layer, bool *extends_to_bottom_flag);
extern struct wetting_front*    listCopy(struct wetting_front* current, struct wetting_front* state_previous=NULL);
extern void                     listDelete(struct wetting_front* head);
extern bool                     GW_fronts_among_surf_WFs(struct wetting_front *head);





/*########################################*/
/*   van Genuchten function prototypes    */
/*########################################*/
/* these are van Genuchten function prototypes.  The actual code lies below the main() function    */
extern double calc_K_from_Se(double Se,double Ks, double m);
extern double calc_h_from_Se(double Se, double alpha, double m, double n);
extern double calc_Se_from_h(double h, double alpha, double m, double n);
extern double calc_theta_from_h(double h, double alpha, double m, double n, double theta_e, double theta_r);
extern double calc_Se_from_theta(double theta,double effsat,double residual);
extern double calc_Geff(bool use_closed_form_G, double theta1, double theta2, double theta_e, double theta_r,
                        double alpha, double n, double m, double h_min, double Ks, int nint, double lambda, double bc_psib_cm);

/*########################################*/
/* LGAR calculation function prototypes   */
/*########################################*/
// computed mass balance
extern double lgar_calc_mass_bal(double *cum_layer_thickness, struct wetting_front* head);

extern double adjust_new_theta(int new_wf_num, double target_mass, double *cum_layer_thickness, int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head);

extern void lgar_clean_redundant_fronts(struct wetting_front** head, int *soil_type, struct soil_properties_ *soil_properties);

extern bool correct_close_psis(int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head);

// computes derivatives; called derivs() in Python code
extern void lgar_dzdt_calc(bool use_closed_form_G, int nint, double timestep_h, double h_p, int *soil_type, double *cum_layer_thickness_cm,
			   double *frozen_factor, struct wetting_front* head, struct soil_properties_ *soil_properties, int num_layers);

// computes dry depth
extern double lgar_calc_dry_depth(bool use_closed_form_G, bool TO_enabled, int nint, double timestep_h, double *deltheta, int *soil_type,
                                  double *cum_layer_thickness_cm, double *frozen_factor,
				  struct wetting_front* head, struct soil_properties_ *soil_properties);

// reads van Genuchten parameters from a file
extern int lgar_read_vG_param_file(char const* vG_param_file_name, int num_soil_types, double wilting_point_psi_cm,
                                    struct soil_properties_ *soil_properties);

// creates a surficial front (new top most wetting front)
extern double lgar_create_surficial_front(bool TO_enabled, int num_layers, double *ponded_depth_cm, double *volin, double dry_depth,
					double theta1, int *soil_type, double *cum_layer_thickness_cm,
					double *frozen_factor, struct wetting_front **head, struct soil_properties_ *soil_properties);

// computes the infiltration capacity, fp, of the soil
extern double lgar_insert_water(bool use_closed_form_G, int nint, double timestep_h, double *free_drainage_subtimestep_cm, double *AET_demand_cm, double *ponded_depth,
				double *volin_this_timestep, double precip_timestep_cm, int wf_free_drainge_demand,
				int num_layers, double ponded_depth_max_cm, int *soil_type, double *cum_layer_thickness_cm,
				double *frozen_factor, struct wetting_front* head, struct soil_properties_ *soil_properties);

// the subroutine moves wetting fronts, merges wetting fronts, and does the mass balance correction if needed
// extern double lgar_move_wetting_fronts(bool TO_enabled, double timestep_h, double *ponded_depth_cm, int wf_free_drainage_demand,
// 				     double old_mass, int number_of_layers, double *actual_ET_demand,
// 				     double *cum_layer_thickness_cm, int *soil_type_by_layer, double *frozen_factor,
// 				     struct wetting_front** head, struct wetting_front* state_previous, struct soil_properties_ *soil_properties);
extern double lgar_move_wetting_fronts(bool TO_enabled, double timestep_h, double *free_drainage_timestep_cm, double PET_timestep_cm, double wilting_point_psi_cm, double field_capacity_psi_cm, double root_zone_depth_cm,
             double *volin_cm, int wf_free_drainage_demand, double old_mass, int num_layers, double surf_frac_rz, double *AET_demand_cm, double *cum_layer_thickness_cm,
				     int *soil_type, double *frozen_factor, struct wetting_front** head,
				     struct wetting_front* state_previous, struct soil_properties_ *soil_properties, double *surf_AET_vec);

// the subroutine merges the wetting fronts; called from lgar_move_wetting_fronts
extern void lgar_merge_wetting_fronts(struct wetting_front** head, int *soil_type, double *frozen_factor, 
				      struct soil_properties_ *soil_properties);

// the subroutine lets wetting fronts cross soil layer boundaries; called from lgar_move_wetting_fronts
extern void lgar_wetting_fronts_cross_layer_boundary(int num_layers, double* cum_layer_thickness_cm,
						     int *soil_type, double *frozen_factor, struct wetting_front** head,
						     struct soil_properties_ *soil_properties);

/* the subroutine allows the deepest wetting front to partially leave the model through the lower boundary if necessary;
   called from lgar_move_wetting_fronts. Currently, fluxes from the lower boundary will always be 0 and this fraction of a
   wetting front will be dealt with in another way */
extern double lgar_wetting_front_cross_domain_boundary(bool TO_enabled, double* cum_layer_thickness_cm, int *soil_type, double *frozen_factor,
						       struct wetting_front** head, struct soil_properties_ *soil_properties);

// subroutine allows surface WFs to merge with TO (groundwater) WFs
extern bool lgar_merge_surface_and_TO_wetting_fronts(bool merged_in_non_top_layer, int num_layers, double* cum_layer_thickness_cm, struct wetting_front** head);

//this does some cleanup which is only necessary after the above fxn, lgar_merge_surface_and_TO_wetting_fronts, runs
extern void lgarto_cleanup_after_surface_TO_merging_in_layer_below_top(bool merged_in_non_top_layer, int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head);

//subroutine allows TO WFs to cross layer boundaries
extern void lgar_TO_wetting_fronts_cross_layer_boundary(int *front_num_with_negative_depth, int num_layers, double *cum_layer_thickness_cm, int *soil_type, double *frozen_factor, struct soil_properties_ *soil_properties, struct wetting_front** head);

//subroutine allows TO WFs to merge with other TO WFs when a higher TO WF gets too deep
extern double lgarto_TO_WFs_merge_via_depth(double initial_mass, double column_depth, double* cum_layer_thickness_cm, struct wetting_front** head, int *soil_type, struct soil_properties_ *soil_properties);

//subroutine allows TO WFs to merge with other TO WFs when a lower TO WF gets too dry. 
//in the rare case that merging based on theta is necessary but the wetting front that would be deleted has to_bottom == TRUE, then instead the mass difference is moved to the GW flux.
extern double lgarto_TO_WFs_merge_via_theta(double initial_mass, double column_depth, double* cum_layer_thickness_cm, struct wetting_front** head, int *soil_type, struct soil_properties_ *soil_properties);

// subroutine to handle wet over dry wetting fronts condtions
extern void lgar_fix_dry_over_wet_wetting_fronts(double *mass_change, double* cum_layer_thickness_cm, int *soil_type,
						 struct wetting_front** head, struct soil_properties_ *soil_properties);

//some of the TO-TO merging fxns update theta via mass balance, so after that, psi must be updated accordingly.
extern void lgar_global_psi_update(int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head);

//sometimes psi is updated, and after that theta must be updated accordingly
extern void lgar_global_theta_update(double bottom_boundary_flux_above_surface_WFs_cm, int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head);

// extern void lgarto_correct_theta_after_domain_boundary_cross(int num_layers, double* cum_layer_thickness_cm, int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head);

//if surface WFs are present, then they must provide the water necesary for the GW fluxes for the TO WFs that are in contact with the surface 
extern double lgarto_extract_TO_GW_flux_from_surface_WFs(double *bottom_boundary_flux_above_surface_WFs_cm, double bottom_boundary_flux_cm, double *AET_demand_cm,
                                                         double* cum_layer_thickness_cm, int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head);

// the subroutine addresses the rare case where, in LGARTO, a TO WF ends up between two surface WFs after surface WFs move; called from lgar_move_wetting_fronts
extern void lgarto_resolve_TO_WF_between_surf_WFs(int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head);

// sometimes, when LGARTO mode is on, a large number of wetting fronts will be generated. This keeps the total number low. 
extern void lgarto_consolidate_excessive_fronts(double* cum_layer_thickness_cm, struct wetting_front** head, int *soil_type, struct soil_properties_ *soil_properties);

// this function is used in lgarto_consolidate_excessive_fronts()
extern int lgarto_count_fronts_for_excessive_calc(struct wetting_front** head);

//corrects negative wetting front depths in LGARTO mode 
extern bool lgarto_correct_negative_depths(struct wetting_front** head);

// checks if dry over wet wetting front exists or not
extern bool lgar_check_dry_over_wet_wetting_fronts(struct wetting_front* head);

// finds free drainage wetting front (the deepest wetting front with psi value closer to zero; saturated in terms of psi)
extern int wetting_front_free_drainage(struct wetting_front* head);

// finds the largest magnitude of dzdt among TO WFs; is currently used for numerical error control
extern double lgarto_calc_largest_mag_TO_dzdt(struct wetting_front* head);

// computes updated theta (soil moisture content) after moving down a wetting front; called for each wetting front to ensure mass is conserved
extern double lgar_theta_mass_balance(int layer_num, int soil_num, double psi_cm, double new_mass,
				      double prior_mass, double *AET_demand_cm, double *delta_theta, double *layer_thickness_cm,
				      int *soil_type, struct soil_properties_ *soil_properties);

/********************************************************************/
// Bmi functions
/********************************************************************/
// functions to initialize model's state at time zero from a config file
extern void lgar_initialize(string config_file, struct model_state *state);
extern void InitFromConfigFile(string config_file, struct model_state *state);
extern vector<double> ReadVectorData(string key);
extern void InitializeWettingFronts(bool TO_enabled, int num_layers, double initial_psi_cm, int *layer_soil_type, double *cum_layer_thickness_cm,
				    double *layer_thickness_cm, double *frozen_factor, struct wetting_front** head, struct soil_properties_ *soil_properties);

/********************************************************************/
/*Other function prototypes for doing hydrology calculations, etc.  */
/********************************************************************/

extern void calc_aet(bool TO_enabled, double PET_timestep_cm, double time_step_h, double wilting_point_psi_cm, double field_capacity_psi_cm, double root_zone_depth_cm, double *surf_frac_rz, 
                       int *soil_type, double AET_thresh_Theta, double AET_expon, struct wetting_front* head,
            		       struct soil_properties_ *soil_properties, double *surf_AET_vec);

extern double calc_aet_for_individual_TO_WFs(int WF_num, double WF_thickness_cm, double rooting_zone_depth, double PET_timestep_cm, double time_step_h, double wilting_point_psi_cm,
		       double field_capacity_psi_cm, int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head);

extern double lgarto_calc_aet_from_TO_WFs(int num_layers, double deepest_surf_depth_at_start, double root_zone_depth, double PET_timestep_cm, double timestep_h, double surf_frac_rz, 
                                          double wilting_point_psi_cm, double field_capacity_psi_cm, int *soil_type, double *cum_layer_thickness_cm, struct soil_properties_ *soil_properties, struct wetting_front** head);

extern void lgarto_ensure_rooting_zone_population(double rzd, double PET_timestep_cm, int *soil_type, struct soil_properties_ *soil_properties, struct wetting_front** head);


extern int lgarto_correction_type(int num_layers, double* cum_layer_thickness_cm, struct wetting_front** head);//returns an integer that describes which type of layer boundary crossing or WF merging is necessary, for all WFs
extern int lgarto_correction_type_surf(int num_layers, double* cum_layer_thickness_cm, struct wetting_front** head);//returns an integer that describes which type of layer boundary crossing or WF merging is necessary, but only for surface WFs

/********************************************************************/
/* Input/Output functions, etc.  */
/********************************************************************/

// computes global mass balance at the end of the simulation
extern void lgar_global_mass_balance(struct model_state *state, double *giuh_runoff_queue);

// writes full state of wetting fronts (depth, theta, no. of wetting front, no. of layer, dz/dt, psi) to a file at each time step
extern void write_state(FILE *out, struct wetting_front* head);


/********************************************************************/
/* Function used in coupling with seasonally frozen soil modules  */
/********************************************************************/

// computes frozen factor for each layer (coefficient used to modify hydraulic conductivity of layers)
extern void frozen_factor_hydraulic_conductivity(struct lgar_bmi_parameters lgar_bmi_params);

/*###################################################################*/
/*   1- and 2-D int and double memory allocation function prototypes */
/*###################################################################*/
extern void itwo_alloc(int ***ptr, int x, int y);
extern void dtwo_alloc(double ***ptr, int x, int y);
extern void d_alloc(double **var,int size);
extern void i_alloc(int **var,int size);
extern void f_alloc(float **var,int size);



/*###############################*/
/*   utility function prototypes */
/*###############################*/
extern bool is_epsilon_less_than(double a, double eps);

#endif  // _ALL_HXX
