#include "../include/all.h"  // <--- This header file contains all function prototypes and other global definitions.


// this is a C++ style comment that tells the compiler to ignore the rest of the line.

//#####################################################################################################################
//  LGAR- Layered Green and Ampt Infiltration with Redistribution.                                                   
//        Patterned after python development code by Peter La Follette, 2022.                                        
//        Adapted to C by Fred Ogden.                                                                                
//  notes:                                                                                                            
//        1. Loop indices start at 1 not 0.                                                                          
//        2. See/change #defines above for maximum number of layers and max. number of wetting fronts and soil defs. 
//        3. Set text editor to 120 columns wide. Do not type anything beyond column 120.                                                                 */
//        4. Do not use tabs in this code.  Emulate tabs with spaces.                                                
//        4. Two spaces per indent level please.  Curley braces indicate change in indentation level.  No }else{     
//
//  C is:  ( here's one of many online references: https://overiq.com/c-programming-101/ )
//  - case sensitive
//  - white space independent, so you can make your code look however you want it to look
//  - because it is white space independent, it requires punctuation.  Single statements must end in a semicolon ;
//  - instead of using indentation to infer punctuation, it uses curly braces {} to impose order.
//  - function based: func(arg1, arg2, ..., arg N)  for almost everything
//  - all variables must be declared before use.  Standard types include: char, int, double, struct, float, etc.
//  - C is dangerous in that it has direct access to memory.  A variable that contains an address in memory is called
//    a POINTER.  If you write to a random address through a pointer you can crash a machine, corrupt a disk or some
//    other nightmarish scenario.
//  - given a variable declared as "double foo", foo contains a double precision number.  That number is stored at
//    an address in memory, you can access it using the ampersand:  &foo  is the address in memory where the number 
//    stored in the variable foo is stored.
//  - if you pass the variable foo to a subroutine, you pass its value.  If you change the value of foo in that 
//    subroutine, that change will not come back to the main program.
//  - if you want to change the value of a variable in a subroutine and return that changed value to the calling
//    program, you need to pass its address &foo.  In the subroutine you access that value as *foo, because foo is
//    a pointer.  But if you change the value stored at &foo (e.g. *foo=2.0*(*foo).  <-notice I needed to use parens
//    to disambiguate the multiplication operator from the "value stored in" operator.  When control returns to the 
//    calling program, the value stored in address &foo, which is called foo in the main program is twice what it was
//    before the routine was called.
//  - all arrays start life as a pointer.   double *arr;  declares a double precision pointer that can be used to 
//    allocate memory for a 1-D array of double precision numbers.  After allocating memory for that array, the values
//    stored at those locations in memory are accessed using array indices:  arr[i].
//  - a two-dimensional array must start life as a double pointer where a double pointer is a pointer to a pointer.
//    for example: double **two_d_array; declares a double pointer (pointer to a pointer) that can be allocated as a
//    two-dimensional array.
//  - the abs() function takes the absolute value of an integer, fabs() function takes the absolute value of a real 
//    number.  Confusing them can cause heartache.
//
//#####################################################################################################################

// Global run-time option flags
bool debug_flag=TRUE;         // TODO set to TRUE to enable screen printing debugging info
bool quiet_flag=FALSE;         // TODO set to TRUE to have code print nothing to screen
bool illiterate_flag=FALSE;    // TODO set to TRUE to make code write nothing to disk


/* global variables: */

struct wetting_front *head = NULL;  //<- this pointer stores the address in memory of the first member of the linked
                                    //   list containing all the wetting fronts.  That address is called "head".  
                                    //   The contents of struct wetting_front are defined in "all.h"

struct wetting_front *head_previous = NULL;
/*################################################*/
/*################################################*/
/*################################################*/
/* main program                                   */
/* everything defined before this point is global */
/*################################################*/
/*################################################*/
/*################################################*/
 
int main(int argc, char *argv) 
{

FILE *in_forcing_fptr;    // pointer to forcing data file
FILE *outl_fptr;           // write output layers to this file pointer after we open it

FILE *outd_fptr;           // write some output data to this file pointer after we open it

// lgar local variables
struct wetting_front *current;
struct wetting_front *previous;
char jstr[256];       // a junk array of characters
char soilname[20];    // a string into which we read the name of each soil, then ignore
char vG_param_file_name[50]; // input file name from which to read in the vG parameters
int num_soil_types;   // the number of soil types with data in the input file
int layer;            // layer counter
int front;            // front counter
int num_soil_layers;  // number of actual soil layers
int load_ICs;         // set to TRUE to read in initial water content profile from a file **fixme**
int *soil_type_by_layer=NULL;  // allocate this to MAX_NUM_SOIL_LAYERS, integer equal to the soil type in each layer 
int num_fronts;       // the current number of fronts
int bottom_flag;      // TRUE iff wetting front is in contact with layer bottom
int nint;             // the number of trapezoids to use in trapezoid rule integral of Geff
int num_dzdt_calculated;    // value returned by dzdt_calc() function
int time_step_num;          // integer time step counter
int forcing_interval;       // the integer number of model timesteps per forcing timestep
int root_zone_bottom_layer; // the lowest soil layer containing roots
int yr,mo,da,hr,mi;   // integer variables to read in year,month,day,hour,minute of forcing data- not yet used

double time_step_s;
double initial_psi_cm;
double wilting_point_psi_cm; 
double Se;                           // this is the relative (scaled 0-1) water content, like Theta 
double *layer_thickness_cm=NULL;     // thickness of individual layers, allocate memory at run time
double *cum_layer_thickness_cm=NULL; // cumulative thickness of layers, allocate memory at run time
double theta_init;                   // initial water content
double fdepth,ftheta;                // fred's test variable ffor creating new wetting fronts
double theta1,theta2;                // usually the limits of integration of Geff or width of a front, theta1 <= theta2 
double dry_depth;
double ponded_depth_cm;              // more explicit variable name, commonly called h_p or H_p
double delta_theta;                  // the width of a front, such that its volume=depth*delta_theta
double precip_timestep_cm =0.0;
double precip_previous_timestep_cm;
double PET_timestep_cm; 
double AET_timestep_cm=0.0;
double AET_thresh_Theta;             // threshold scaled moisture content (0-1) above which AET=PET
double AET_expon;                    // exponent that allows curvature of the rising portion of the Budyko curve

// mass balance variables, expressed as a depth of water (cm)
double volprecip = 0.0;   //  cumulative amount of precip
double volstart  = 0.0;   //  volume of water in the soil at the beginning of the simulation
double volend    = 0.0;   //  volume of water in the soil at the end of the simulation
double volin     = 0.0;   //  volume of precip added to the ponded depth
double volrunoff = 0.0;   //  volume of water removed from the surface as runoff
double volon     = 0.0;   //  volume of water remaining on the surface at the end of the sim.
double volAET     = 0.0;   //  cumulative amount of actual ET
double volPET    = 0.0;   //  cumulative amount of potential ET
double volrech   = 0.0;   //  cumulative amount of water leaving the bottom of the soil
// note: at the end of each time step the following must be true (assuming volon=0 at t=0):
// volstart+volprecip-volAET-volin-volrunoff-volon-volrech-volend = 0.0
 
bool error;

//###############################
// Forcing variables
double precip_mm_per_15min;
double PET_mm_per_15min;

/*##########################*/
/* needed for plot routines */
/*--------------------------*/
bool PLOT=FALSE; // set to false iff you don't want to plot output to an x-window 
long nplt=1;    // the plotting frequency, set to 2 to update plot every second time step, ffor instance 
int chrome;     // the color of the water in the plot.   I like blue, which is color number 6. 
int refresh_flag=TRUE;  // erases prior profile before plotting new one 
int maxfronts=MAX_NUM_WETTING_FRONTS;      // should be enough
int num_fronts_to_plot; // should be less than maxfronts
int yes_stop=FALSE;      // change to TRUE ffor debugging purposes
double xmin;            // smallest theta_r in all layers
double xmax;            // largest theta_e in all layers
double ymax;            // depth of all soil layers, cm.  Note ymin=0.0, so we don't need to save it

// NOTE: The stuff betwen #ifdef XWINDOWS and #endif included in the source iff UNIX is defined, ellse it's ignored
//       this is called a "macro", that uses the C pre-processor to include source or ignore it.

#ifdef XWINDOWS 
//####################################                                     ###################
// X-windows stuff                   *                                     ###### FOR ########
// initialize the scaling parameters *                                     ###################
//-----------------------------------*                                     ###   #######   ###
winw = 970;         // was 770   controls width of window                  ##### ####### #####
winh = 1200;         // was 455   controls height of window                ###### ##### ######
xorg=70;            // controls x origin coord of plot                     ####### ### #######
yorg=70;            // controls y origin coord of plot                     ######## # ########
xline=480;          // controls length of x axis                           ######### #########
yline=700;          // controls length of y axis                           ######## # ########
yaxdiff=15.0;       //                                                     ####### ### #######                                    
ticlen=6;           //                                                     ###### ##### ######
ptextlen=2;         //                                                     ##### ####### #####
rescaleflag=1;      //                                                     ###   #######   ###
int *topleftx = NULL; // x-coord (pixel) of top left corner of rectangle   ###################
int *toplefty = NULL; // y-coord (pixel) of top left corner of rectangle   ##### WINDOWS #####
int *fheight  = NULL; // height of rectangle (pixels)                      ###### CODE #######
int *fwidth   = NULL; // width of rectangle (pixels)                       ###################
i_alloc(&topleftx,maxfronts);                      //                      ###################        
i_alloc(&toplefty,maxfronts);                      //                      ###################           
i_alloc(&fwidth,maxfronts);                        //                      ###################    
i_alloc(&fheight,maxfronts);                       //                      ###################      
extern void xinitialize(int shi, int nola);
extern int  xcategorize(int xmin, int xmax, int ymax, int num_soil_layers, int *soil_type_by_layer,
                        double *cum_layer_thickness_cm, struct soil_properties_ *soil_properties, 
                        int *topleftx, int *toplefty, int *fwidth, int *fheight);
extern void xplot(int shade, int num_fronts_to_plot, int xorg, int yorg, int *topleftx, int *toplefty,
                  int *fwidth, int *fheight, int refresh_flag);
#endif

/* END OF X-WINDOWS Initialization */
/*#################################*/



//###############################
// Allocate array memory
i_alloc(&soil_type_by_layer,MAX_NUM_SOIL_LAYERS);
d_alloc(&layer_thickness_cm,MAX_NUM_SOIL_LAYERS);
d_alloc(&cum_layer_thickness_cm,MAX_NUM_SOIL_LAYERS);


//###############################
// run-time options

//debug_flag=FALSE;                 // set to TRUE to print lots of debugging output to the screen.
 int num_time_steps= 90001;//68804; //41020; //24700; //90001;           // total simulation time is this number times the time_step_s
// int num_time_steps= 2;//288;          // total simulation time is this number times the time_step_s
nint=120;                         // the number of trapezoids used in integrating the Geff function
strcpy(vG_param_file_name,"data/vG_default_params.dat");  // the name of the soil parameter input file
 //double forcing_resolution_s=300.0;  // 15 minute data
double forcing_resolution_s=3600.0;  // 15 minute data

//model time step:  30.0000 s, forcing data time interval: 30 timesteps.
//forcing data interval =  30, 30,  30.0000
//###################################################################
// INITIALIZATIONS
time_step_s=300.0;
// time_step_s=3600.0;
double time_step_h=time_step_s/3600;

num_soil_types=15;          // must be less than or equal to MAX_NUM_SOIL_TYPES 
num_soil_layers=3;          // must be less than or equal to MAX_NUM_SOIL_LAYERS 
num_fronts=0;               // must be less than MAX_NUM_FRONTS
root_zone_bottom_layer = 2; // example, iff set to 2, then the root zone is in layer 1 and 2

// **** PICK TEST *******//
 bool sim_test = 1; // 0: for synthetic and 1 for real(phillipsburg)
 if (sim_test == 0)
   initial_psi_cm= 50.0;        // cm
 else
   initial_psi_cm=2000;       // cm - phillipsburg
 
AET_thresh_Theta = 0.85;    // scaled soil moisture (0-1) above which AET=PET
AET_expon = 1.0;            // exponent that allows curvature of the rising portion of the Budyko curve
// check? wilting_point_psi_cm=initial_psi_cm*101325.0/9810.0*100.0;  // assumed 15 atm.
wilting_point_psi_cm = 154950/10.; // AJ-correction
// printf("wilting::: %lf \n", wilting_point_psi_cm);
load_ICs=FALSE;      // set to false iff initial conditions are not specified TODO fixeme iff TRUE

/*## EXAMPLE ##*/
 if (sim_test == 0) {
   layer_thickness_cm[1]=10.0;
   layer_thickness_cm[2]=10.0;
   layer_thickness_cm[3]=10.0;
 }
 else {
   // Philiphsburg example
   layer_thickness_cm[1]=44.0;
   layer_thickness_cm[2]=131.0;
   layer_thickness_cm[3]=25.0;
 }

soil_type_by_layer[1]=13;  // sand    
soil_type_by_layer[2]=14;  // sandy loam
soil_type_by_layer[3]=15;  // clay



//########################################################################################
// option sanity checks.  Detect nonsensical or infeasible options or initializations here
if(num_soil_types > (MAX_NUM_SOIL_TYPES -1))
  {printf("Too many soil types specified.  Increase MAX_NUM_SOIL_TYPES in all.h .  Program stopped.\n");exit(0);}
if(num_soil_layers > (MAX_NUM_SOIL_LAYERS -1))
  {printf("Too many soil types specified.  Increase MAX_NUM_SOIL_LAYERS in code.  Program stopped.\n");exit(0);}


//##########################################################################################
// calculate the cumulative (absolute) depth from land surface to bottom of each soil layer
cum_layer_thickness_cm[0]=0.0;
// note, since only one command associated with this loop, don't need curly braces {}
for(layer=1;layer<=num_soil_layers;layer++) 
  cum_layer_thickness_cm[layer]=layer_thickness_cm[layer]+cum_layer_thickness_cm[layer-1]; 
      

//allocate memory to create an array of structures to hold the soils properties data.
struct soil_properties_ *soil_properties = (struct soil_properties_*) 
                              malloc((num_soil_types+1)*sizeof(struct soil_properties_));

lgar_read_vG_param_file(vG_param_file_name, num_soil_types, wilting_point_psi_cm, soil_properties);


//###################################################################
// calculate the initial theta and wilting point moisture content
// ffor the soil in each layer from initial_psi and assumed wilting point psi
// and create initial fronts and include them in each of the soil layers

int soil;
front=0;
for(layer=1;layer<=num_soil_layers;layer++)
  {
  front++;
  soil = soil_type_by_layer[layer];
  theta_init = calc_theta_from_h(initial_psi_cm,soil_properties[soil].vg_alpha_per_cm,
                                 soil_properties[soil].vg_m,soil_properties[soil].vg_n,
                                 soil_properties[soil].theta_e,soil_properties[soil].theta_r);

  // the next lines create the initial moisture profile
  bottom_flag=TRUE;  // all intial wetting fronts are in contact with the bottom of the layer they exist in
  // NOTE: The listInsertFront function does lots of stuff.
  current = listInsertFront(cum_layer_thickness_cm[layer],theta_init,front,layer,bottom_flag);
  current->psi_cm = initial_psi_cm;
  Se = calc_Se_from_theta(current->theta,soil_properties[soil].theta_e,soil_properties[soil].theta_r);
  current->K_cm_per_s = calc_K_from_Se(Se, soil_properties[soil].Ksat_cm_per_s, soil_properties[soil].vg_m);  // cm/s
  
  }

//listPrint();
// abort();
//###################################################################

#ifdef XWINDOWS   
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for(layer=1;layer<=num_soil_layers;layer++)
  {
  soil=soil_type_by_layer[layer];
  xmin=10.0;
  xmax=-1.0;

  if(soil_properties[soil].theta_r < xmin) xmin = soil_properties[soil].theta_r;
  if(soil_properties[soil].theta_e > xmax) xmax = soil_properties[soil].theta_e;
  }
ymax=cum_layer_thickness_cm[num_soil_layers];
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif 

if(debug_flag==TRUE)   // this stuff is all demo/test code.
  {
    printf("Initial linked list of wetting fronts:\n"); 
    //printf("depth,theta,layer,front,to_bottom,dzdt ");
    //listPrint();
  // note in the next line, using the file print function (fprintf), with arg. "stdout" causes write to screen
  fprintf(stdout,"Initial depth of water stored in soil: %lf cm\n",volstart);
  fprintf(stdout,"Should be equal to:                    %lf cm\n",39.1658320808);  // a kind of unit test.
  fprintf(stdout,"difference:                            %e cm\n",volstart-39.1658320808);
  printf("Initial linked list of wetting fronts:\n");
  printf("after calculating dzdt for %d fronts.\n",num_dzdt_calculated); 
  printf("depth,theta,layer,front,to_bottom,dzdt ");
  listPrint();
  }

if(illiterate_flag==FALSE)
  {
  // open the output file ffor writing
  if((outl_fptr=fopen("data_layers.out","w"))==NULL)
    {printf("Problem opening output file. Program stopped.\n"); exit(0);}
  if((outd_fptr=fopen("data_variables.out","w"))==NULL)
    {printf("Problem opening output file. Program stopped.\n"); exit(0);}
  // add headings to variables forcing file
  fprintf(outd_fptr,"precip(mm),AET(mm),runoff(mm),ponded_depth(mm),storage(mm),bottom_flux(mm),mass_bal_err(mm)\n");
  }

// open the forcing file ffor reading
//if((in_forcing_fptr=fopen("forcing_data_syn_case3.txt","r"))==NULL)
// if((in_forcing_fptr=fopen("forcing_data_syn_case4.txt","r"))==NULL)
//if((in_forcing_fptr=fopen("Phillipsburg_data1.txt","r"))==NULL)
 if((in_forcing_fptr=fopen("forcing/Phillipsburg_data2_PET.txt","r"))==NULL)

  {printf("Problem opening forcing_data_syn.txt input file. Program stopped.\n"); exit(0);}

forcing_interval=(int) (forcing_resolution_s/(time_step_s)+1.0e-08); // add 1.0e-08 to prevent truncation error
 if(debug_flag)
   {
     fprintf(stdout,"model time step: %8.4lf s, forcing data time interval: %d timesteps.\n",time_step_s,forcing_interval);
   }
 
 if( ((forcing_interval) % (int)(time_step_s+1.0e-08)) != 0) // here % is the remainder operator
  {
    //printf("forcing data interval (s) must be integer divisible by the time step (s).  Program stopped.\n");
    //exit(-1);  // -1 is arbitrary, and can be targeted.  Usually int functions return 0 if successful.
  }

fgets(jstr,255,in_forcing_fptr);  // read the header line and ignore.

#ifdef XWINDOWS
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(PLOT)
  {
  int j;
  /* initialize the xlib display functions */

  xinitialize(winw,winh);

  /* plot the initial conditions */
  num_fronts_to_plot=xcategorize(xmin, xmax, ymax, num_soil_layers, soil_type_by_layer,
                                 cum_layer_thickness_cm, soil_properties, 
                                 topleftx, toplefty, fwidth, fheight);  // find the corners of the rectangles to draw

  refresh_flag=TRUE; //erases window before drawing.    
  
  // plot the wetting fronts in an x-window
  xplot(chrome, num_fronts_to_plot, xorg,yorg, topleftx,toplefty,fwidth,fheight,refresh_flag); 
  }
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
  
// if(yes_stop) getchar();  //uncomment to have the code stop and wait ffor someone to press enter.


// REMINDER OF MASS BALANCE TOTALIZERS
// volprecip
// volstart 
// volend   
// volin    
// volrunoff
// volon    
// volAET    
// volrech 
// note: at the start and end of each time step the following must be true (assuming volon=0 at t=0):
// volstart+volprecip-volAET-volin-volrunoff-volon-volrech-volend = 0.0 


// sum up the initial depth of water stored in the soil
volstart=lgar_calc_mass_bal(num_soil_layers,cum_layer_thickness_cm);

ponded_depth_cm=0.0;   // always start simulation with no ponded depth

if(debug_flag)
  printf("\n ************ ************ TIME STEPPING NOW ****************** \n");
double mm_s_to_cm_hr = 1. * (1.0/10.); // mm per sec to cm per hour conversion

double volend_timestep_cm;
double volstart_timestep_cm;
double volrech_timestep_cm;
double volprecip_timestep_cm;
double volin_timestep_cm = 0.0;
double volrunoff_timestep_cm = 0.0;

time_t result = time(NULL);
clock_t start_time, end_time;
double elapsed;
start_time = clock();

for(time_step_num=0;time_step_num<num_time_steps;time_step_num++)
  {

    if(debug_flag) {
      printf("\n-------------------------------TS start----------------------------------------------- \n");
      printf(" ************ time step %d ******** \n",time_step_num);
    }
    
    state_previous = NULL;
    state_previous = listCopy(head);
    
    if(debug_flag) {
      //printf("Previous state : \n");
      //listPrintC(*state_previous);
    }
     
    // ----------------- read forcing data (precip and PET) -----------------------
    if(time_step_num % forcing_interval == 0)  // necessary because model and forcing timesteps differ
      {
	fscanf(in_forcing_fptr,"%d %d %d %d %d %lf %lf",&yr,&mo,&da,&hr,&mi,
	       &precip_mm_per_15min,&PET_mm_per_15min);
	if(debug_flag)
	  fprintf(stderr,"Time = %d %d %d %d %d %lf %lf\n",yr,mo,da,hr,mi,precip_mm_per_15min,PET_mm_per_15min);
	 
      }
     
    if(!is_epsilon_less_than(precip_mm_per_15min, 1.0e-06))
      {
	error=0;
      }
    
    precip_timestep_cm = precip_mm_per_15min/forcing_interval * mm_s_to_cm_hr;  // 10 converts from mm to cm
    PET_timestep_cm    = PET_mm_per_15min/forcing_interval * mm_s_to_cm_hr;     // 10 converts from mm to cm

     
    if (debug_flag)
      printf("per timestep (precip, pet) = (%lf, %lf) \n ", precip_timestep_cm,PET_timestep_cm );
    double AET_timestep_temp_cm = 0.0;
    
    // ------------------ iff raining, take some or all of the PET from rainfall depending on rain rate --------------
    /*
     if(!is_epsilon_less_than(precip_mm_per_15min, 1.0e-06))  // in essence, iff rainrate != 0, note ! negates return val.
       {
	 //abort();
	 if(precip_timestep_cm < PET_timestep_cm )
	   {
	     volET += precip_timestep_cm;  // in this ccase all rainfall evaporates but leaves some PET demand
	     PET_timestep_cm -= precip_timestep_cm;  // note that in this ccase some PET remains unsatisfied
	     precip_timestep_cm = 0.0;                    // but all precip is used to partially satisfy PET
	     // AET_timestep_temp_cm = 
	   }
	 else
	   {
	     volET += PET_timestep_cm;
	     precip_timestep_cm -= PET_timestep_cm;  // reduce rainfall by PET
	     PET_timestep_cm = 0.0;                       // all PET is satisfied
	   }
       }
     */
     
    if (PET_timestep_cm>0) {
      // Calculate AET from PET and root zone soil moisture.  Note PET was reduced iff raining
      volPET+= PET_timestep_cm;
      
      AET_timestep_cm = calc_aet(PET_timestep_cm, time_step_h, wilting_point_psi_cm, soil_properties, soil_type_by_layer, AET_thresh_Theta, AET_expon);
      
      // this part needs to be tested....
     //printf("AET = %lf %lf %0.6e \n", PET_mm_per_15min, PET_timestep_cm, AET_timestep_cm);
     //abort();
    }
    else
      AET_timestep_cm = 0.0;

    // printf("AET = %lf %lf %6.10f \n", PET_mm_per_15min, PET_timestep_cm, AET_timestep_cm);
    // mass balance
    volprecip += precip_timestep_cm * time_step_h;
    volstart_timestep_cm = lgar_calc_mass_bal(num_soil_layers,cum_layer_thickness_cm);
    volprecip_timestep_cm = precip_timestep_cm * time_step_h;
    
    // -------------------- add incoming precip to ponded_depth_cm ---------------------------
    //    printf("---- Ponded depth = %lf %lf \n ", ponded_depth_cm, precip_timestep_cm * time_step_h);
    //if (ponded_depth_cm>0) abort();
    ponded_depth_cm += precip_timestep_cm * time_step_h;
     
    int wf_free_drainage_demand = wetting_front_free_drainage();

    //printf("wf_free_drainage_demand = %d \n", wf_free_drainage_demand);
    //bool surficial_wf_created = false;
    
    volin_timestep_cm =0.0;
    //total_mass_added += precip_timestep_cm * time_step_h;
    
    int soil_num = soil_type_by_layer[head->layer_num];
    
    double theta_e = soil_properties[soil_num].theta_e;  // rhs of the new front, assumes theta_e as per Peter
    
    bool is_top_wf_saturated = head->theta >= theta_e ? true : false;
    bool create_surficial_front = (precip_previous_timestep_cm == 0.0 && precip_timestep_cm >0.0);
    
    double mass_source_to_soil_timestep = 0.0;
    //printf ("__**___Surficial front create? = %d %d %lf %lf \n ",  create_surficial_front, is_top_wf_saturated, head->theta, theta_e);
    
    if(create_surficial_front && !is_top_wf_saturated)  
      //  This means that there is no wetting front in the top layer to accept the water, must create one.
      {
	double temp_pd = 0.0; // necessary to assign zero precip due to the creation of new wetting front; AET will still be taken out of the layers
	 
	lgar_move_wetting_fronts(&temp_pd,time_step_h, wf_free_drainage_demand, volend_timestep_cm, num_soil_layers, &AET_timestep_cm, cum_layer_thickness_cm, soil_type_by_layer, soil_properties);
	//listPrint();
	dry_depth = lgar_calc_dry_depth(nint, time_step_h, soil_type_by_layer, soil_properties, cum_layer_thickness_cm,&delta_theta);
	//printf("SF = %lf %lf ",AET_timestep_cm, dry_depth);
	theta1 = head->theta;
	lgar_create_surfacial_front(&ponded_depth_cm, &volin_timestep_cm, dry_depth, theta1, soil_type_by_layer, soil_properties, cum_layer_thickness_cm, nint, time_step_h);
	
	//listPrint();
	state_previous = NULL;
	state_previous = listCopy(head);
	
	volin += volin_timestep_cm;
	//volrunoff += ponded_depth_cm;
	//printf("surfical front created and ponded depth = %lf %lf %lf \n", ponded_depth_cm, volrunoff, volin_timestep_cm);
	//if (ponded_depth_cm >0)
	//listPrint();
	// abort();
      }


    //     if (time_step_num >= 0 && precip_this_timestep_cm >0 && precip_previous_timestep_cm >0) {
    if (ponded_depth_cm > 0 && !create_surficial_front) {
      //printf("--- LGAR MAIN **** insert water........ %lf, %d \n", ponded_depth_cm, create_surficial_front);
      volrunoff_timestep_cm = lgar_insert_water(&ponded_depth_cm, &volin_timestep_cm, precip_timestep_cm, dry_depth, nint, time_step_h, wf_free_drainage_demand, soil_type_by_layer, soil_properties, cum_layer_thickness_cm);
      
      //volin_this_timestep = ponded_depth_cm;
      volin += volin_timestep_cm;
      volrunoff += volrunoff_timestep_cm;
      volrech_timestep_cm = volin_timestep_cm; // water leaving thru the bottom
      volon = ponded_depth_cm;
      //printf("Mass in = %lf %lf %lf \n", volin, ponded_depth_cm, runoff_timestep);
      if (volrunoff_timestep_cm < 0)
	abort();
      
    }
    else {
      //printf("wetting front created = %lf %d \n", ponded_depth_cm ,!create_surficial_front );
      double hp_cm_max = 0.0; //h_p_max = 0.0;
      
      if (ponded_depth_cm < hp_cm_max) {
	volrunoff += 0.0;
	volon = ponded_depth_cm;
	ponded_depth_cm = 0.0;
	volrunoff_timestep_cm=0.0;
      }
      else {
	volrunoff_timestep_cm=(ponded_depth_cm - hp_cm_max);
	volrunoff += (ponded_depth_cm - hp_cm_max);
	volon = hp_cm_max;
	ponded_depth_cm = hp_cm_max;
	
      }
     }
    //printf("Runoff = %6.12f \n",volrunoff_timestep_cm);
    //printf("LGAR MAIN ******************* potential infiltration ........%lf \n", ponded_depth_cm);
    //lgar_potential_infiltration(time_step_s, cum_layer_thickness_cm, soil_type_by_layer, soil_properties);
    
    if (!create_surficial_front) {
      //printf("LGAR MAIN ******************* move wetting front ........%lf %lf \n", ponded_depth_cm, volin_this_timestep);
      //lgar_move_fronts(&ponded_depth_cm,time_step_h, wf_free_drainage_demand, volend_timestep, cum_layer_thickness_cm, soil_type_by_layer, soil_properties);
      
      lgar_move_wetting_fronts(&volin_timestep_cm,time_step_h, wf_free_drainage_demand, volend_timestep_cm, num_soil_layers, &AET_timestep_cm, cum_layer_thickness_cm, soil_type_by_layer, soil_properties);
      
      // this is the volume of water leaving through the bottom
      volrech_timestep_cm = volin_timestep_cm;
    }
     
#ifdef XWINDOWS
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if(PLOT)
       {
	 
	 /* plot conditions */
	 num_fronts_to_plot=xcategorize(xmin, xmax, ymax, num_soil_layers, soil_type_by_layer,
					cum_layer_thickness_cm, soil_properties,
					topleftx, toplefty, fwidth, fheight);  // find the corners of the rectangles to draw
	 
	 refresh_flag=TRUE; //erases window before drawing.    
	 
	 // plot the wetting fronts in an x-window
	 xplot(chrome, num_fronts_to_plot, xorg,yorg, topleftx,toplefty,fwidth,fheight,refresh_flag); 
       }
     //getchar();  // rigged demo
     //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif

     //printf("Calling dz/dt from main.... \n");
     // calculate dzdt for each front : derivs in python
     
     //listPrint();
     int num_dzdt_calculated = lgar_dzdt_calc(nint, soil_type_by_layer, soil_properties, cum_layer_thickness_cm,
					      ponded_depth_cm);
     if(debug_flag)
       listPrint();
     
     
     // keep a copy of the previous LL (needed?)
     previous = head;
     
     // move fronts
     //     collisions with layers
     //     mergers of overtaking fronts
     
     volAET += AET_timestep_cm;
     
     volrech += volrech_timestep_cm; 
       
     volend_timestep_cm = lgar_calc_mass_bal(num_soil_layers,cum_layer_thickness_cm);
     // volstart+volprecip-volAET-volin-volrunoff-volon-volrech-volend = 0.0
     //volend+= volend_timestep;
     //printf("mass balance = %lf %lf %lf %lf %lf \n", volstart, volin, volin_timestep_cm, volstart+volin, volend_timestep_cm);
     //printf ("mass balance = %lf %lf \n ", volrech, volrech_timestep_cm);
     precip_previous_timestep_cm = precip_timestep_cm;
     //fprintf(out_fptr,"#TS = %5d %7.6f %7.6f\n",time_step_num, precip_mm_per_15min, PET_mm_per_15min);
     
     if(illiterate_flag==FALSE) {
       fprintf(outl_fptr,"#TS = %5d %7.6f %7.6f\n",time_step_num, precip_timestep_cm, AET_timestep_cm);
       
       write_state(outl_fptr);
     }
     
     //volstart+volprecip-volET-volin-volrunoff-volon-volrech-volend = 0.0 
	 //double mass_added_timestep = precip_this_timestep_cm * time_step_h;
	 double local_mb = volstart_timestep_cm + volprecip_timestep_cm - volin_timestep_cm*0 - volrunoff_timestep_cm - AET_timestep_cm - volon -volrech_timestep_cm -volend_timestep_cm;

     if(debug_flag) {
       printf("local mass balance = %0.10e %0.10e %0.10e %0.10e %0.10e %0.10e \n", local_mb, volstart_timestep_cm,volprecip_timestep_cm, volrunoff_timestep_cm,AET_timestep_cm, volend_timestep_cm);
     }
     if (fabs(local_mb) >1e-7) {
       printf("Local mass balance (in this timestep) is %.6e, larger than expected, needs some debugging...\n ",local_mb);
       abort();
     }
     //if (volrech_timestep_cm > 0)
     //  abort();
     if(illiterate_flag==FALSE)
       fprintf(outd_fptr,"%6.10f,%6.10f,%6.10f,%6.10f,%6.10f,%6.10f,%6.20f\n",precip_timestep_cm*10, AET_timestep_cm*10, volrunoff_timestep_cm*10, volon*10, volend_timestep_cm*10 , volrech_timestep_cm*10, local_mb*10);
     
   }

 end_time = clock();

 elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC;
 
 volend = lgar_calc_mass_bal(num_soil_layers,cum_layer_thickness_cm);

 
 //############  END TIME LOOP
  
 fclose(in_forcing_fptr);
 
 if(!illiterate_flag)
   {
     fclose(outl_fptr);
     fclose(outd_fptr);
   }

 // volin is the volume infiltrates into the soil; volrunoff + volin = total volume injected/added
 double global_error_cm = volstart + volprecip - volrunoff - volAET - volon -volrech -volend;


 printf("\n---------------------- Simulation Summary  ------------------------ \n");
 printf("Time (sec)                 = %6.10f \n", elapsed);
 printf("-------------------------- Mass balance ------------------- \n");
 printf("initial water in soil      = %14.10f cm\n", volstart);
 printf("total precipitation input  = %14.10f cm\n", volprecip);
 printf("total infiltration         = %14.10f cm\n", volin);
 printf("final water in soil        = %14.10f cm\n", volend);
 printf("water remaining on surface = %14.10f cm\n", volon);
 printf("surface runoff             = %14.10f cm\n", volrunoff);
 printf("total percolation          = %14.10f cm\n", volrech);
 printf("total AET                  = %14.10f cm\n", volAET);
 printf("total PET                  = %14.10f cm\n", volPET);
 printf("global balance             =   %.6e cm\n", global_error_cm);
 /*
 printf("\n---------------------- Simulation Summary  ------------------------ \n");
 printf("Time (sec)               = %6.10f \n", elapsed);
 printf("-------------------------- Mass balance ------------------- \n");
 printf("vol start                = %6.10f \n", volstart);
 printf("total_mass_added         = %6.10f \n", volprecip);
 printf("vol added (to soil)      = %6.10f \n", volin);
 printf("vol end   (in soil)      = %6.10f \n", volend);
 printf("vol ponded (on surface)  = %6.10f \n", volon);
 printf("vol runoff (surface)     = %6.10f \n", volrunoff);
 printf("vol runoff (bottom)      = %6.10f \n",  volrech);
 printf("vol AET                  = %6.10f \n", volET);
 printf("Global balance           = %.6e \n", global_mb);
 */
 // double global_mb = volstart + (volin+volrunoff) - volET - volon-volrech-volend;





 
 
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// The rest of the code in main() below this line provides examples on usage of linked list functions.
//    it is all commented out.
//-----------------------------------
// EXAMPLE how to delete all fronts
//
//while(!listIsEmpty())
//  {            
//  struct wetting_front *temp = listDeleteFirst();
//  printf("\nDeleted value:");
//  printf("(%d,%d) ",temp->front_num, temp->layer_num);
//  }  
//
//printf("\nList after deleting all items: ");
//listPrint();

//-------------------------------------------------
// EXAMPLE how to find a front and print its values
//
//struct wetting_front *foundLink = listFindFront(3); // find front 3
//
//if(foundLink != NULL)
//  {
//  printf("Element found: ");
//  printf("(front: %d in layer: %d  at Z=%8.4lf with theta=%8.7lf) ",foundLink->front_num,foundLink->layer_num,
//                                                                    foundLink->depth,foundLink->theta);
//  printf("\n");  
//  }
//else
//  {
//  printf("Element not found.");
//  }

//----------------------------------
// EXAMPLE how to delete a front 
//
//listDeleteFront(3);  // argument here is wetting front number
//printf("List after deleting wetting front 3: ");
//listPrint();
//printf("\n");

//--------------------------------------------------------------------------------------------
// EXAMPLE how to sort a linked list of fronts by depth from surface down to the wetting front
//
// printf("\n");
// listSortFrontsByDepth();
//
// printf("List after sorting the layer depth: ");
// listPrint();

//-----------------------------------------------------------------
// EXAMPLE how to sort a linked list in terms of decreasing depth
//
// listReverseOrder(&head);
// printf("\nList after reversing the list: ");
// listPrint();

// EXAMPLE how to delete all fronts 
//
// delete all fronts
// while(!listIsEmpty())
//   {            
//   struct wetting_front *temp = listDeleteFirst();
//   printf("\nDeleted value:");
//   printf("(%d,%d) ",temp->front_num,temp->layer_num);
//   }  
//
// printf("\nList after deleting all items: ");
// listPrint();   

//--------------------------------------
// EXAMPLE how to add a front to the beginning of the list
//
// bottom_flag=FALSE;
// listInsertFirst(10.2,0.28,1,1,bottom_flag);
// after inserting a new link at position 3 */
// listPrint();
}

/*################################################*/
/*################################################*/
/*################################################*/
/* END OF MAIN PROGRAM                            */
/*################################################*/
/*################################################*/
/*################################################*/

extern void write_state(FILE *out){

  struct wetting_front *current = head;

  fprintf(out, "[");
  while(current != NULL)
  {
    if (current == head)
      fprintf(out,"(%lf,%lf,%d,%d,%lf)",current->depth_cm*10., current->theta, current->layer_num,current->front_num, current->psi_cm*10.);
    else
      fprintf(out,"|(%lf,%lf,%d,%d,%lf)",current->depth_cm*10., current->theta, current->layer_num,current->front_num, current->psi_cm*10.);
  current = current->next;
  }
  fprintf(out, "]\n");
  //fprintf(out, "],");

}

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
  /*
do
  {
  if((current->theta > wettest_theta_in_root_zone) && (current->layer_num <= root_zone_bottom_layer) ) 
    {
    wettest_theta_in_root_zone = current->theta;
    layer = current->layer_num;
    if(layer<=0) {printf("problem in calc_aet.  Program stopped.\n"); exit(-1);}  // should never happen
    soil = soil_type[layer];
    }
  current = current->next;
  } while (current != NULL);
 printf("wettest theta = %lf %d \n", wettest_theta_in_root_zone, root_zone_bottom_layer);
// calculate AET from PET and soil moisture
Se    = calc_Se_from_theta(wettest_theta_in_root_zone, soil_props[soil].theta_e, soil_props[soil].theta_r);

Se_wp = calc_Se_from_theta(soil_props[soil].theta_wp, soil_props[soil].theta_e, soil_props[soil].theta_r);

if( Se > AET_thresh_Theta)
  {
  return (PET_timestep_cm);
  }
if( Se <= Se_wp)
  {
  return 0.0;
  }
return (PET_timestep_cm * pow((Se - Se_wp)/(AET_thresh_Theta - Se_wp),AET_expon));
*/
}


