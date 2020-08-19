#ifndef SYL_H
#define SYL_H

#include "typedefs.h"

/**
 * Solves AX + XB = scale * C.
 * */
void syl(
    double sgn,
    int m, int n,
    const double *restrict A, int ldA, const double anorm,
    const double *restrict B, int ldB, const double bnorm,
    double *restrict C, int ldC, double *restrict cnorm,
    scaling_t *restrict scale,
    const double smin);


void blocked_syl(
    double sgn,
    int m, int n,
    const double *restrict A, int ldA, const double anorm,
    const double *restrict B, int ldB, const double bnorm,
    double *restrict C, int ldC, double *restrict cnorm,
    scaling_t *restrict scale,
    const double smin);

#endif
