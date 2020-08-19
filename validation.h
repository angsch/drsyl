#ifndef VALIDATION_H
#define VALIDATION_H

#include "typedefs.h"


void validate(
    double sgn,
    int m, int n,
    const double *restrict const A, int ldA,
    const double *restrict const B, int ldB,
    double *restrict C, int ldC,
    double *restrict X, int ldX,
    const scaling_t scale);


void validate_tiled(
    double sgn,
    double ***A_tiles, int ldA, partitioning_t *part_A,
    double ***B_tiles, int ldB, partitioning_t *part_B,
    double ***C_tiles, int ldC, partitioning_t *part_C,
    double ***X_tiles, int ldX,
    const scaling_t scale,
    memory_layout_t mem_layout);


int validate_quasi_triangular_shape(
    int n, const double *restrict const A, int ldA, matrix_desc_t type);


int validate_spectra(double sgn,
    int m, double *lambda_A, int *lambda_type_A,
    int n, double *lambda_B, int *lambda_type_B);

#endif
