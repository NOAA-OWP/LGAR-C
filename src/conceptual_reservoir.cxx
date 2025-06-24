#ifndef CR_CXX_INCLUDE
#define CR_CXX_INCLUDE

#include "../include/all.hxx"

//A single nonlinear reservoir is used to represent water stored in the catchment that contributes to streamflow, which lumps both quickflow and groundwater contributions to streamflow
extern double calc_CR_Q(double subtimestep_h, double a, double b, double precip_for_CR_subtimestep_cm, double *CR_storage_cm){
    double QCR = subtimestep_h * (a * pow(*CR_storage_cm, b));
    if (*CR_storage_cm < 0.01){ // idea here is that we do want the contribution to streamflw to actually become 0 in arid or semi arid environments
        QCR = 0.0;
    }
    if (*CR_storage_cm + (subtimestep_h * precip_for_CR_subtimestep_cm - QCR) > 0.0){
        *CR_storage_cm += (subtimestep_h * precip_for_CR_subtimestep_cm - QCR);
    }
    else {
        QCR = *CR_storage_cm + precip_for_CR_subtimestep_cm;
        *CR_storage_cm = 0.0;
    }
    return(QCR);
}

#endif

