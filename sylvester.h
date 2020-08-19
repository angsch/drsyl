#ifndef SYLVESTER_H
#define SYLVESTER_H

#include "typedefs.h"

void solve_tiled_sylvester(
    double sgn,
    double ***A_tiles, int ldA,
    double ***B_tiles, int ldB,
    double ***C_tiles, int ldC,
    partitioning_t *p,
    scaling_t *scale,
    memory_layout_t mem_layout);

#endif
