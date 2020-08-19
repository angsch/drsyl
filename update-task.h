#ifndef UPDATE_TASK_H
#define UPDATE_TASK_H

#include <omp.h>


//  C := C - alpha * A  * B.
void update(
    int m, int n, int k,
    omp_lock_t *lock,
    const double alpha, double *restrict const Ain, int ldAin,
    const double ainnorm, const scaling_t ainscale,
    double *restrict const B, int ldB, const double bnorm,
    double *restrict const C, int ldC, double *restrict const cnorm, 
    scaling_t *restrict const cscale);

#endif
