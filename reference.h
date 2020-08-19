#ifndef REFERENCE_H
#define REFERENCE_H

#include "typedefs.h"

void solve_sylvester_dtrsyl(
    const int transA, const int transB,
    double sgn,
    double ***A_tiles, int ldA,
    double ***B_tiles, int ldB,
    double ***C_tiles, int ldC,
    partitioning_t *p);

#endif
