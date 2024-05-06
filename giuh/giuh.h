#ifndef _GIUH_H

#define _GIUH_H

extern double giuh_convolution_integral(int adaptive_timestep, double subtimestep_h, double runoff_m, int num_giuh_ordinates, 
                                   double *giuh_ordinates, double *runoff_queue_m_per_timestep);

#endif
