#include "../include/all.hxx"
#include "iostream"
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



/***********************************************************************************************/
/* This function calculates the unsaturated capillary drive Geff(0i,0o)  (note "0" are thetas) */
/* for the Green and Ampt redistribution function.                                             */
/* used for redistribution following the equation published by Ogden and Saghafian 1995.       */
/* to compile: "gcc one_block.c -lm"                                                           */
/*                                                                                             */
/* author: Fred Ogden, June, 2021, edited by Peter La Follette in 2023                         */
//updated to calculate G using an adaptive integral method, rather than always using a fixed number of trapezoids.
//This really helps when the limits of integration in terms of psi are quite far apart, and also helps to save runtime.
/***********************************************************************************************/

// extern double calc_Geff(bool use_closed_form_G, double theta1, double theta2, double theta_e, double theta_r,
//                         double vg_alpha, double vg_n, double vg_m, double h_min, double Ksat, int nint,
// 			double lambda, double bc_psib_cm)
// {
//   double Geff;       // this is the result to be returned.

//   if (use_closed_form_G==false) {
//     // local variables
//     // note: units of h in cm.  units of K in cm/s
//     double h_i,h_f,Se_i,Se_f;  // variables to store initial and final values
//     double Se;
//     double h2;         // the head at the right-hand side of the trapezoid being integrated [m]
//     double dh;         // the delta h over which integration is performed [m]
//     double Se1,Se2;    // the scaled moisture content on left- and right-hand side of trapezoid
//     double K1,K2;      // the K(h) values on the left and right of the region dh integrated [m]


//     Se_i = calc_Se_from_theta(theta1,theta_e,theta_r);    // scaled initial water content (0-1) [-]
//     Se_f = calc_Se_from_theta(theta2,theta_e,theta_r);    // scaled final water content (0-1) [-]

//     h_i = calc_h_from_Se(Se_i,vg_alpha,vg_m,vg_n);  // capillary head associated with Se_i [cm]
//     h_f = calc_h_from_Se(Se_f,vg_alpha,vg_m,vg_n);  // capillary head associated with Se_f [cm]

//     if(h_i < h_min) {/* if the lower limit of integration is less than h_min FIXME */
//       //return h_min; // commenting out as this is not used in the Python version
//     }

//     if (verbosity.compare("high") == 0) {
//       // debug statements to see if calc_Se_from_h function is working properly
//       Se = calc_Se_from_h(h_i,vg_alpha,vg_m,vg_n);
//       printf("Se_i = %8.6lf,  Se_inverse = %8.6lf\n", Se_i, Se);

//       Se = calc_Se_from_h(h_f,vg_alpha,vg_m,vg_n);
//       printf("Se_f = %8.6lf,  Se_inverse = %8.6lf\n", Se_f, Se);
//     }
    
//     // nint = number of "dh" intervals to integrate over using trapezoidal rule
//     dh = (h_f-h_i)/(double)nint;
//     Geff = 0.0;
    
//     // integrate k(h) dh from h_i to h_f, using trapezoidal rule, with subscript
//     // 1 denoting the left-hand side of the trapezoid, and 2 denoting the right-hand side

//     Se1 = Se_i;  // could just use Se_i in next statement.  Done 4 completeness.
//     K1  = calc_K_from_Se(Se1, Ksat, vg_m);
//     h2  = h_i + dh;
    
//     for(int i=1; i<=nint; i++) {

//       Se2 = calc_Se_from_h(h2, vg_alpha, vg_m, vg_n);
//       K2  = calc_K_from_Se(Se2, Ksat, vg_m);
//       Geff += (K1+K2)*dh/2.0;                  // trapezoidal rule
//       // reset for next time through loop
//       K1 = K2;
//       h2 += dh;
//     }
    
//     Geff = fabs(Geff/Ksat);       // by convention Geff is a positive quantity

//     if (verbosity.compare("high") == 0){
//       printf ("Capillary suction (G) = %8.6lf \n", Geff);
//     }

//     return Geff;

//   }
//   else {

//     double Se_f = calc_Se_from_theta(theta1,theta_e,theta_r);    // the scaled moisture content of the wetting front
//     double Se_i = calc_Se_from_theta(theta2,theta_e,theta_r);    // the scaled moisture content below the wetting front
//     double H_c  = bc_psib_cm*((2+3*lambda) / (1+3*lambda));         // Green ampt capillary drive parameter, which can be used in the approximation of G with the Brooks-Corey model (See Ogden and Saghafian, 1997)

//     Geff = H_c*(pow(Se_i,(3+1/lambda))-pow(Se_f,(3+1/lambda)))/(1-pow(Se_f,(3+1/lambda)));
    
//     if (isinf(Geff)) {
//       Geff = H_c;
//     }
//     if (isnan(Geff)) {
//       Geff = H_c;
//     }
    
//     return Geff;
    
//   }
  
// }



extern double calc_Geff(bool use_closed_form_G, double theta1, double theta2, double theta_e, double theta_r,
                        double vg_alpha, double vg_n, double vg_m, double h_min, double Ksat, int nint, double lambda, double bc_psib_cm)

{
  double Geff;       // this is the result to be returned.

  if (use_closed_form_G==false){
    // local variables
    // note: units of h in cm.  units of K in cm/s
    double h_i,h_f,Se_i,Se_f;  // variables to store initial and final values
    double Se;
    double h2;         // the head at the right-hand side of the trapezoid being integrated [m]
    double dh;         // the delta h over which integration is performed [m]
    double Se1,Se2;    // the scaled moisture content on left- and right-hand side of trapezoid
    double K1,K2;      // the K(h) values on the left and right of the region dh integrated [m]

    Se_i = calc_Se_from_theta(theta1,theta_e,theta_r);    // scaled initial water content (0-1) [-]
    Se_f = calc_Se_from_theta(theta2,theta_e,theta_r);    // scaled final water content (0-1) [-]

    h_i = calc_h_from_Se(Se_i,vg_alpha,vg_m,vg_n);  // capillary head associated with Se_i [cm]
    h_f = calc_h_from_Se(Se_f,vg_alpha,vg_m,vg_n);  // capillary head associated with Se_f [cm]

    if(h_i < h_min) {/* if the lower limit of integration is less than h_min FIXME */
      //return h_min; // commenting out as this is not used in the Python version
    }

    if (verbosity.compare("high") == 0) {
      // debug statements to see if calc_Se_from_h function is working properly
      Se = calc_Se_from_h(h_i,vg_alpha,vg_m,vg_n);
      printf("Se_i = %8.6lf,  Se_inverse = %8.6lf\n", Se_i, Se);

      Se = calc_Se_from_h(h_f,vg_alpha,vg_m,vg_n);
      printf("Se_f = %8.6lf,  Se_inverse = %8.6lf\n", Se_f, Se);
    }

    // nint = number of "dh" intervals to integrate over using trapezoidal rule
    dh = -1*(h_f-h_i)/(double)nint;
    dh = dh*0.01;

    Geff = 0.0;

    // integrate k(h) dh from h_i to h_f, using trapezoidal rule, with subscript
    // 1 denoting the left-hand side of the trapezoid, and 2 denoting the right-hand side

    Se1 = Se_i;  // could just use Se_i in next statement.  Done 4 completeness.
    K1  = calc_K_from_Se(Se1, Ksat, vg_m);
    h2  = h_f + dh;

    while(h2<h_i) {

      double prior_h2 = h2;

      Se2 = calc_Se_from_h(h2, vg_alpha, vg_m, vg_n);
      K2  = calc_K_from_Se(Se2, Ksat, vg_m);

      if ( (K1-K2)/K2 > 0.02 ){//if K1 disagrees with K2 by more than this fraction, then dh is made smaller
        dh = dh*0.5;
      }
      else if ( (K1-K2)/K2 < 0.03 ){//but if K1 and K2 are within a certain fraction of each other, then dh is made larger
        dh = dh*10.0;
      }
      else{
        dh = dh*2.0; //but if K1 and K2 are within a certain fraction of each other, then dh is made larger
      }

      if (h2<h_i){
        Geff += (K1+K2)*dh/2.0;                  // trapezoidal rule
      }
      else{
        dh = h_i - prior_h2;
        Se2 = calc_Se_from_h(h_i, vg_alpha, vg_m, vg_n);
        K2  = calc_K_from_Se(Se2, Ksat, vg_m);
        Geff += (K1+K2)*dh/2.0;  
      }

      // reset for next time through loop
      K1 = K2;
      h2 += dh;

    }

    //std::cerr<<"Integral = "<< Geff<<" "<<Ksat<<"\n";
    Geff = fabs(Geff/Ksat);       // by convention Geff is a positive quantity

    if (verbosity.compare("high") == 0){
      printf ("Capillary suction (G) = %8.6lf \n", Geff);
    }

    return Geff;

    //if(Geff < h_min) Geff = h_min; AJ // as per Morel-Seytoux and Khanji

    //return Geff;

  }
   else {

     double Se_f = calc_Se_from_theta(theta1,theta_e,theta_r);    // the scaled moisture content of the wetting front
     double Se_i = calc_Se_from_theta(theta2,theta_e,theta_r);    // the scaled moisture content below the wetting front
     double H_c = bc_psib_cm*((2+3*lambda)/(1+3*lambda));            // Green ampt capillary drive parameter, which can be used in the approximation of G with the Brooks-Corey model (See Ogden and Saghafian, 1997)
     Geff = H_c*(pow(Se_i,(3+1/lambda))-pow(Se_f,(3+1/lambda)))/(1-pow(Se_f,(3+1/lambda)));
     if (isinf(Geff)){
       Geff = H_c;
     }
     if (isnan(Geff)){
       Geff = H_c;
     }

     return Geff;

  }

}



/**************************************/
/* function to calculate theta from h */
/**************************************/
double calc_theta_from_h(double h,double alpha, double m, double n, double theta_e, double theta_r)
{
  return(1.0/(pow(1.0+pow(alpha*h,n),m))*(theta_e-theta_r)+theta_r);
}

/***********************************/
/* function to calculate Se from h */
/***********************************/
double calc_Se_from_h(double h,double alpha, double m, double n)
{
  if(is_epsilon_less_than(h,1.0e-01)) return 1.0;  // this function doesn't work well ffor tiny h
  else return(1.0/(pow(1.0+pow(alpha*h,n),m)));
}

/***********************************/
/* function to calculate K from Se */
/***********************************/
double calc_K_from_Se(double Se, double Ksat, double m)
{
  return (Ksat * sqrt(Se) * pow(1.0 - pow(1.0 - pow(Se,1.0/m), m), 2.0));  // same units as Ksat
}

/***********************************/
/* function to calculate h from Se */
/***********************************/
double calc_h_from_Se(double Se, double alpha, double m, double n)
{
  return(1.0/alpha*pow(pow(Se,-1.0/m)-1.0,1.0/n));
}

/***************************************/
/* function to calculate Se from theta */
/***************************************/
double calc_Se_from_theta(double theta,double e,double r)
{
  return((theta-r)/(e-r));
}
