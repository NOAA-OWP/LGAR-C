#ifndef _ALL_H

#define _ALL_H


#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <assert.h>
#include <stdbool.h>
#include "X_stuff.h"
#include <time.h>

#define TRUE 1
#define FALSE 0
#define ONE 1

#define use_bmi_flag FALSE       // TODO set to TRUE to run in BMI environment

#define MAX_NUM_SOIL_LAYERS 4
#define MAX_NUM_SOIL_TYPES 16
#define MAX_SOIL_NAME_CHARS 25
#define MAX_NUM_WETTING_FRONTS 300


// DEFINE A DATA STRUCTURE TO HOLD EVERYTHING THAT DESCRIBES A WETTING FRONT
struct wetting_front 
  {
  double depth_cm;  // depth down from the land surface (absolute depth)
  double theta;         // water content of the soil moisture block
  double psi_cm;        // psi calculated at rhs of the current wetting front
  double K_cm_per_s;    // the value of K(theta) associated with the wetting front
  int layer_num;        // the layer containing this wetting front.
  int front_num;        // the wetting front number (might be irrelevant), but useful to debug
  bool to_bottom;       // TRUE iff this wetting front is in contact with the layer bottom
  double dzdt_cm_per_s; // use to store the calculated wetting front speed
  struct wetting_front *next;  // pointer to the next wetting front.
  };




// DEFINE A DATA STRUCTURE TO HOLD PROPERTIES and PARAMETERS FOR EACH SOIL TYPE 
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
  double Ksat_cm_per_s;    // saturated hydraulic conductivity cm/s
  double theta_wp;         // water content at wilting point [-]
  };





extern struct wetting_front *head;  // GLOBALLY DEFINED pointer to the first link in the wetting front list.
                                    // Making it a local variable in main() makes all linked list operations
                                    // in subroutines a pain of referencing.  Since it is just one thing, 
                                    // making it global just makes everything easier.

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


extern void                     listPrint();
extern void                     listPrintC(struct wetting_front pcurrent);
extern int                      listLength();
extern bool                     listIsEmpty();
extern struct wetting_front*    listDeleteFirst();
extern struct wetting_front*    listFindFront(int i,struct wetting_front* head_old);
extern struct wetting_front*    listDeleteFront(int i);
extern void                     listSortFrontsByDepth();
extern void                     listInsertFirst(double d, double t, int f, int l, bool b);
extern struct wetting_front*    listInsertFront(double d, double t, int f, int l, bool b);
extern struct wetting_front*    listInsertFrontAtDepth(int numlay, double *tvec,double d, double t);
extern void                     listReverseOrder(struct wetting_front** head_ref);
extern bool                     listFindLayer( struct wetting_front* link, int num_layers, 
                                               double *cum_layer_thickness_cm,
                                               int *lives_in_layer, bool *extends_to_bottom_flag);
extern struct wetting_front*    listCopy();

struct wetting_front *state_previous;
extern struct wetting_front *head_previous; //head pointer to the previous state, used in computing derivatives

/*########################################*/
/*   van Genuchten function prototypes    */
/*########################################*/
/* these are van Genuchten function prototypes.  The actual code lies below the main() function    */
extern double calc_K_from_Se(double Se,double Ks, double m);
extern double calc_h_from_Se(double Se, double alpha, double m, double n);
extern double calc_Se_from_h(double h, double alpha, double m, double n);
extern double calc_theta_from_h(double h, double alpha, double m, double n, double theta_e, double theta_r);
extern double calc_Se_from_theta(double theta,double effsat,double residual);
extern double calc_Geff(double theta1, double theta2, double theta_e, double theta_r, 
                        double alpha, double n, double m, double h_min, double Ks, int nint);

/*########################################*/
/* LGAR calculation function prototypes   */
/*########################################*/
extern double lgar_calc_mass_bal(int num_soil_layers, double *cum_layer_thickness);

extern int lgar_dzdt_calc(int nint, int *soil_type, struct soil_properties_ *soil_properties, 
                     double *cum_layer_thickness, double h_p);  // called derivs() in Python code
extern double lgar_calc_dry_depth(int nint, double time_step_s, int *soil_type, 
                                  struct soil_properties_ *soil_properties, double *cum_layer_thickness_cm,
                                  double *deltheta);
extern void lgar_read_vG_param_file(char *vG_param_file_name, int num_soil_types, double wilting_point_psi_cm,
                                    struct soil_properties_ *soil_properties);
extern void lgar_create_surfacial_front(double *ponded_depth_cm, double *volin, double dry_depth, double theta1,
                                        int *soil_type, struct soil_properties_ *soil_properties, 
                                        double *cum_layer_thickness_cm, int nint, double dt);
extern double lgar_insert_water(double *ponded_depth, double *volin_this_timestep, double precip_timestep_cm, double dry_depth, int nint, double time_step_s, int wf_free_drainge_demand,
                              int *soil_type, struct soil_properties_ *soil_properties, double *cum_layer_thickness_cm);
extern void lgar_move_fronts(double *ponded_depth_cm, double time_step_s, int wf_free_drainage_demand, double old_mass, int number_of_layers, double *actual_ET_demand, double *cum_layer_thickness_cm, int *soil_type_by_layer, struct soil_properties_ *soil_properties);

extern double lgar_theta_mass_balance(int layer_num, int soil_num, double psi_cm, double new_mass, double prior_mass,
					double depth_cm_old, double *delta_theta,
					double *layer_thickness_cm, int *soil_type, struct soil_properties_ *soil_properties);

extern int wetting_front_free_drainage();
/********************************************************************/
/*Other function prototypes for doing hydrology calculations, etc.  */
/********************************************************************/

extern double calc_aet(double PET_timestep_cm, double time_step_s, double wilting_point_psi_cm,
                       struct soil_properties_ *soil_props, int *soil_type, double AET_thresh_Theta, double AET_expon);                             
extern void write_state(FILE *out);

/*###################################################################*/
/*   1- and 2-D int and double memory allocation function prototypes */
/*###################################################################*/
extern void itwo_alloc( int ***ptr, int x, int y);
extern void dtwo_alloc( double ***ptr, int x, int y);
extern void d_alloc(double **var,int size);
extern void i_alloc(int **var,int size);
extern void f_alloc(float **var,int size);


/*###############################*/
/*   utility function prototypes */
/*###############################*/
extern bool is_epsilon_less_than(double a, double eps);

#endif  // _ALL_H
