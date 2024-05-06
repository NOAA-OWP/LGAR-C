#ifndef _GIUH_C

#define _GIUH_C

#include "giuh.h"
#include <stdio.h>
#include <math.h>


//##############################################################
//############### GIUH CONVOLUTION INTEGRAL   ##################
//##############################################################
extern double giuh_convolution_integral(int adaptive_timestep, double subtimestep_h, double runoff_m,int num_giuh_ordinates, 
					double *giuh_ordinates, double *runoff_queue_m_per_timestep)
{
  //##############################################################
  // This function solves the convolution integral involving N
  //  GIUH ordinates.
  //##############################################################

  if (adaptive_timestep){//adaptive time step case

    double runoff_m_now;
    int N,i,j;
    
    N=num_giuh_ordinates;
    runoff_queue_m_per_timestep[N]=0.0;
    
    if (subtimestep_h==1.0/12.0){
      for(i=0;i<N;i++)
        {
          runoff_queue_m_per_timestep[i]+=giuh_ordinates[i]*runoff_m;
        }
      
      runoff_m_now=runoff_queue_m_per_timestep[0];
      
      for(i=1;i<N;i++)  // shift all the entries in preperation ffor the next timestep
        {
          runoff_queue_m_per_timestep[i-1]=runoff_queue_m_per_timestep[i];
        }
      runoff_queue_m_per_timestep[N-1]=0.0;
    }
    else {//The time step is not equal to the smallest internal time step. In this case, the giuh code must be repeated j times because the time between each entry in runoff_queue_m_per_timestep is fixed.
      double factor = (1.0/12.0)/subtimestep_h;
      runoff_m = runoff_m*factor;
      int resample_limit = fmax(1,(int) 1.0/factor);
      for (j=0;j<resample_limit;j++){
        for(i=0;i<N;i++)
          {
            runoff_queue_m_per_timestep[i]+=giuh_ordinates[i]*runoff_m;
          }
        
        runoff_m_now+=runoff_queue_m_per_timestep[0];
        
        for(i=1;i<N;i++)  // shift all the entries in preperation ffor the next timestep
          {
            runoff_queue_m_per_timestep[i-1]=runoff_queue_m_per_timestep[i];
          }
        runoff_queue_m_per_timestep[N-1]=0.0;
      }

    }
    
    return runoff_m_now;
  }

  else {//fixed time step case
    double runoff_m_now;
    int N,i;
    
    N=num_giuh_ordinates;
    runoff_queue_m_per_timestep[N]=0.0;
    
    for(i=0;i<N;i++)
      {
        runoff_queue_m_per_timestep[i]+=giuh_ordinates[i]*runoff_m;
      }
    
    runoff_m_now=runoff_queue_m_per_timestep[0];
    
    for(i=1;i<N;i++)  // shift all the entries in preperation ffor the next timestep
      {
        runoff_queue_m_per_timestep[i-1]=runoff_queue_m_per_timestep[i];
      }
    runoff_queue_m_per_timestep[N-1]=0.0;
    
    return runoff_m_now;
  }
}


#endif
