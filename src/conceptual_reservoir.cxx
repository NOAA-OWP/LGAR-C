#ifndef CR_CXX_INCLUDE
#define CR_CXX_INCLUDE

#include "../include/all.hxx"


extern double calc_CR_Q(
    double subtimestep_h,
    double a_fast, double a_slow,
    double b_fast, double b_slow,
    double frac_slow,  // fraction (0 - 1) of recharge going to slow reservoir
    double precip_for_CR_subtimestep_cm_per_h,
    double *CR_fast_storage_cm,
    double *CR_slow_storage_cm)
{
    // Partition recharge between fast and slow reservoirs
    double input_slow = precip_for_CR_subtimestep_cm_per_h * frac_slow;
    double input_fast = precip_for_CR_subtimestep_cm_per_h - input_slow; // implicit (1 - frac_slow)

    // === FAST reservoir outflow ===
    double Q_fast = subtimestep_h * (a_fast * pow(*CR_fast_storage_cm, b_fast));
    if (*CR_fast_storage_cm < 0.01) Q_fast = 0.0;

    double delta_fast = subtimestep_h * input_fast - Q_fast;
    if (*CR_fast_storage_cm + delta_fast > 0.0) {
        *CR_fast_storage_cm += delta_fast;
    } else {
        Q_fast = *CR_fast_storage_cm + subtimestep_h * input_fast;
        *CR_fast_storage_cm = 0.0;
    }

    // === SLOW reservoir outflow ===
    double Q_slow = subtimestep_h * (a_slow * pow(*CR_slow_storage_cm, b_slow));
    if (*CR_slow_storage_cm < 0.01) Q_slow = 0.0;

    double delta_slow = subtimestep_h * input_slow - Q_slow;
    if (*CR_slow_storage_cm + delta_slow > 0.0) {
        *CR_slow_storage_cm += delta_slow;
    } else {
        Q_slow = *CR_slow_storage_cm + subtimestep_h * input_slow;
        *CR_slow_storage_cm = 0.0;
    }

    return Q_fast + Q_slow;
}

#endif
