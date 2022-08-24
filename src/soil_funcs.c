#include "../include/all.h"

/*##################################################*/
/*##################################################*/
/*##################################################*/
//--------------------------------------------------
//
//        SSS         OOO     LL    LL  
//      SSS  SSS    OOO OOO   TT    LL
//     SS      SS  OO     OO  TT    LL
//     SS          OO     OO  TT    LL
//      SSSS       OO     OO  TT    LL
//          SSS    OO     OO  TT    LL
//            SSS  OO     OO  TT    LL
//     SS      SS  OO     OO  TT    LL 
//      SSS  SSS    OOO OOO   TT    LL
//        SSSS        OOO     TT    LLLLLLLL
//
//  Property Calculations from Constitutive Relations



extern double calc_Geff(double theta1, double theta2, double theta_e, double theta_r,
                        double vg_alpha, double vg_n, double vg_m, double h_min, double Ksat, int nint)
{
/***********************************************************************************************/
/* This function calculates the unsaturated capillary drive Geff(0i,0o)  (note "0" are thetas) */
/* for the Green and Ampt redistribution function.                                             */
/* used for redistribution following the equation published by Ogden and Saghafian 1995.       */
/* to compile: "gcc one_block.c -lm"                                                           */
/*                                                                                             */
/* author: Fred Ogden, June, 2021,                                                             */
/***********************************************************************************************/

// local variables
// note: units of h in cm.  units of K in cm/s
  
  int i,local_debug_flag=FALSE;
  double h_i,h_f,Se_i,Se_f,K_i,K_f;  // variables to store initial and final values
  double Se;
  double h2;      // the head at the right-hand side of the trapezoid being integrated [m]
  double dh;      // the delta h over which integration is performed [m]
  double Se1,Se2; // the scaled moisture content on left- and right-hand side of trapezoid
  double K1,K2;   // the K(h) values on the left and right of the region dh integrated [m]
  double Geff;    // this is the result to be returned.
  
  if(theta2<=theta1 && false)
    {
      printf("in function calc_Geff(), theta2 <=theta 1, when theta1 must be less than theta2.\n");
      printf("theta1=%lf, theta2=%lf\n",theta1, theta2);
      printf("returning 0.\n");
      //return(0.0); AJ
    }
  
  Se_i = calc_Se_from_theta(theta1,theta_e,theta_r);    // scaled initial water content (0-1) [-]
  Se_f = calc_Se_from_theta(theta2,theta_e,theta_r);    // scaled final water content (0-1) [-]
  
  //printf("Se_ : %lf %lf %lf %lf \n", Se_i, Se_f, theta1, theta2);
  //printf("SeX_ : %lf %lf %lf %lf \n", vg_alpha, vg_n, vg_m, Ksat);
  
  h_i = calc_h_from_Se(Se_i,vg_alpha,vg_m,vg_n);  // capillary head associated with Se_i [cm]
  h_f = calc_h_from_Se(Se_f,vg_alpha,vg_m,vg_n);  // capillary head associated with Se_f [cm]
  //printf("Inside G limits: %lf %lf %lf \n", h_i, h_f, h_min);
  if(h_i < h_min) return h_min;   /* if the lower limit of integration is less than h_min FIXME */
  
  if(local_debug_flag)
    {
      // set debug_flag to TRUE to see iff calc_Se_from_h function is working properly
      Se = calc_Se_from_h(h_i,vg_alpha,vg_m,vg_n);
      //printf("Se_i=%8.6lf  Se_inverse=%8.6lf\n",Se_i,Se);
      Se = calc_Se_from_h(h_f,vg_alpha,vg_m,vg_n);
      //printf("Se_f=%8.6lf  Se_inverse=%8.6lf\n",Se_f,Se);
    }
  
  // nint = number of "dh" intervals to integrate over using trapezoidal rule
  dh = (h_f-h_i)/(double)nint;
  Geff = 0.0;
  
  // integrate k(h) dh from h_i to h_f, using trapezoidal rule, with subscript
  // 1 denoting the left-hand side of the trapezoid, and 2 denoting the right-hand side
  
  Se1= Se_i;  // could just use Se_i in next statement.  Done 4 completeness.
  K1 = calc_K_from_Se(Se1,Ksat,vg_m);
  h2 = h_i+dh;
  //printf("G params : %lf %lf %lf %lf %lf \n", vg_alpha, vg_n, vg_m, Se1, K1*3600*10);
  for(i=1;i<=nint;i++)
    {
      
      Se2 = calc_Se_from_h(h2,vg_alpha, vg_m, vg_n);
      K2 = calc_K_from_Se(Se2,Ksat,vg_m);
      Geff += (K1+K2)*dh/2.0;                  // trapezoidal rule
      // reset ffor next time through loop
      //printf("integ: %d %6.8f %6.8f %6.8f %6.8f\n", i, Geff*10, K1*10, K2*10, dh);
      K1 = K2;
      h2 += dh;
    }
  //printf("Integral: %lf %lf %lf \n", Geff, Ksat, Geff/Ksat);
  Geff = fabs(Geff/Ksat);       // by convention Geff is a positive quantity
  //if(Geff < h_min) Geff = h_min; AJ // as per Morel-Seytoux and Khanji
  //printf ("Inside G DONE = %lf %lf \n", Geff, h_min);
  
  return(Geff);
  
}


double calc_theta_from_h(double h,double alpha, double m, double n, double theta_e, double theta_r)
{
/**************************************/
/* function to calculate theta from h */
/**************************************/

return(1.0/(pow(1.0+pow(alpha*h,n),m))*(theta_e-theta_r)+theta_r);
}

double calc_Se_from_h(double h,double alpha, double m, double n)
{
/***********************************/
/* function to calculate Se from h */
/***********************************/
if(is_epsilon_less_than(h,1.0e-01)) return 1.0;  // this function doesn't work well ffor tiny h
else return(1.0/(pow(1.0+pow(alpha*h,n),m)));
}

double calc_K_from_Se(double Se, double Ksat, double m)
{
/***********************************/
/* function to calculate K from Se */
/***********************************/
return(Ksat*sqrt(Se)*pow( 1.0-    pow(1.0-   pow(Se,1.0/m)  ,m)     ,2.0));  // same units as Ksat 
}

double calc_h_from_Se(double Se, double alpha, double m, double n)
{
/***********************************/
/* function to calculate h from Se */
/***********************************/
return(1.0/alpha*pow(pow(Se,-1.0/m)-1.0,1.0/n));
}

double calc_Se_from_theta(double theta,double e,double r)
{
/***************************************/
/* function to calculate Se from theta */
/***************************************/
return((theta-r)/(e-r));
}


