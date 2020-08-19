#ifndef RANDOM_H
#define RANDOM_H

#include "defines.h"
#include "typedefs.h"

int random_integer(int low, int high);
double random_double(double low, double high);

void generate_eigenvalues(const int n, double complex_ratio,
    double *restrict const lambda, int *restrict const lambda_type);

void generate_upper_quasi_triangular_matrix(
    double *** T_blocks, int ld,
    partitioning_t *p,
    const double *restrict const lambda, const int *restrict const lambda_type,
    memory_layout_t mem_layout);

void generate_dense_matrix(
    double *** A_blocks, int ld,
    partitioning_t *p,
    memory_layout_t mem_layout);
#endif
