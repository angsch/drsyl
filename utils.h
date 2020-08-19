#ifndef UTILS_H
#define UTILS_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "typedefs.h"

int get_size_with_padding(const int n);

static inline int is_aligned(const void *restrict const ptr)
{
    return ((uintptr_t)ptr % ALIGNMENT == 0);
}

static inline int min(int a, int b)
{
    return a < b ? a : b;
}

static inline int max(int a, int b)
{
    return a > b ? a : b;
}

static inline double minf(double x, double y)
{
    return x < y ? x : y;
}

static inline double maxf(double x, double y)
{
    return x > y ? x : y;
}

static inline double find_absmax_in_2x2(const double *restrict const A, int ldA)
{
    double largest = maxf(fabs(A[0]), fabs(A[1]));
    largest = maxf(largest, fabs(A[ldA]));
    largest = maxf(largest, fabs(A[ldA + 1]));

    return largest;
}

static inline scaling_t min_element(int n, const scaling_t *restrict const a)
{
    scaling_t smallest;

#ifdef INTSCALING
    smallest = 0;
#else
    smallest = INFINITY;
#endif

    for (int i = 0; i < n; i++)
        if (a[i] < smallest)
            smallest = a[i];

    return smallest;
}


static inline void set_zero(int n, int m, double *restrict const A, int ldA)
{
    for (int j = 0; j < m; j++)
        for (int i = 0; i < n; i++)
            A[(size_t) j * ldA + i] = 0.0;

}

static inline void copy_block(int n, int m,
    const double *restrict const Ain, int ldA,
    double *restrict const Bout, int ldB)
{
#define Ain(i,j) Ain[(i) + ldA * (j)]
#define Bout(i,j) Bout[(i) + ldB * (j)]

    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            Bout(i,j) = Ain(i,j);
        }
    }

#undef Ain
#undef Bout
}


void copy_matrix(
    double ***in_blocks, int ldin,
    double ***out_blocks, int ldout,
    partitioning_t *p,
    memory_layout_t layout);


void scale(
    int n, double *restrict const x, const scaling_t *beta);

void scale_tile(int m, int n, 
    double *restrict const X, int ldX, const scaling_t *beta);

void scale_excluding(int m, int n, int ilo, int ihi, int jlo, int jhi,
    double *restrict const X, int ldX, const scaling_t *beta);

double convert_scaling(scaling_t alpha);

void init_scaling_factor(
    int n, scaling_t *restrict const alpha);

double compute_upscaling(
    scaling_t alpha_min, scaling_t alpha);

double compute_combined_upscaling(
    scaling_t alpha_min, scaling_t alpha, scaling_t beta);

void update_global_scaling(
    scaling_t *global, scaling_t phi);

void update_norm(double *norm, scaling_t phi);

/// Copy all selected eigenvalues to a compact memory representation.
void compact_eigenvalues(int n, const int *restrict selected,
    const double *restrict lambda, const int *restrict lambda_type,
    double *restrict compact_lambda, int *restrict compact_lambda_type);


void dgemm(
    const char transa, const char transb,
    const int m, const int n, const int k,
    const double alpha, const double *restrict const A, const int ldA,
    const double *restrict const B, const int ldB,
    const double beta, double *restrict C, const int ldC);


double dlange(const char norm, const int m, const int n,
    const double *restrict const A, int ldA);

#endif
