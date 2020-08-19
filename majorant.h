#ifndef MAJORANT_H
#define MAJORANT_H

#include "typedefs.h"

void bound_triangular_matrix(
    double ***A_tiles, int ldA,
    double *restrict A_norms,
    int num_tiles, int *p,
    matrix_desc_t type,
    memory_layout_t mem_layout);


// Compute upper bound for each rhs.
void compute_column_majorants(
    int m, int n,
    const double *restrict C, int ldC,
    double *restrict C_norms);

#endif
