#include "utils.h"

#include <stdio.h>
#include <math.h>


int get_size_with_padding(const int n)
{
#if defined(__AVX512F__)
    // Round up n to the next multiple of 8.
    return 8 * ((n + 7) / 8);
#else
    // Round up n to the next multiple of 4.
    return 4 * ((n + 3) / 4);
#endif
}


void copy_matrix(
    double ***in_blocks, int ldin,
    double ***out_blocks, int ldout,
    partitioning_t *p,
    memory_layout_t layout)
{
    int num_blk_rows = p->num_blk_rows;
    int num_blk_cols = p->num_blk_cols;
    int *first_row = p->first_row;
    int *first_col = p->first_col;

    // Copy out_blocks := in_blocks.
    for (int j = 0; j < num_blk_cols; j++) {
        for (int i = 0; i < num_blk_rows; i++) {
            const int num_rows = first_row[i + 1] - first_row[i];
            const int num_cols = first_col[j + 1] - first_col[j];

            // Copy block out_i,j := in_i,j.
            if (layout == COLUMN_MAJOR) {
                copy_block(num_rows, num_cols, in_blocks[i][j], ldin,
                    out_blocks[i][j], ldout);
            }
            else { //     TILE_LAYOUT
                copy_block(num_rows, num_cols, in_blocks[i][j], num_rows,
                    out_blocks[i][j], num_rows);
            }
        }
    }
}




void scale(
    int n, double *restrict const x, const scaling_t *beta)
{
#ifdef INTSCALING
    double alpha = ldexp(1.0, beta[0]);
#else
    double alpha = beta[0];
#endif

    // Scale vector, if necessary.
    if (alpha != 1.0) {
        for (int i = 0; i < n; i++) {
            x[i] = alpha * x[i];
        }
    }
}


void scale_tile(int m, int n, 
    double *restrict const X, int ldX, const scaling_t *beta)
{
#ifdef INTSCALING
    double alpha = ldexp(1.0, beta[0]);
#else
    double alpha = beta[0];
#endif

#define X(i,j) X[(i) + (j) * ldX]

    // Scale vector, if necessary.
    if (alpha != 1.0) {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                X(i,j) = alpha * X(i,j);
            }
        }
    }

#undef X
}


// Scale everything that is not enclosed in (ilo:ihi,jlo:jhi)
void scale_excluding(int m, int n, int ilo, int ihi, int jlo, int jhi,
    double *restrict const X, int ldX, const scaling_t *beta)
{
#ifdef INTSCALING
    double alpha = ldexp(1.0, beta[0]);
#else
    double alpha = beta[0];
#endif

#define X(i,j) X[(i) + (j) * ldX]

    // Scale vector, if necessary.
    if (alpha != 1.0) {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if ((ilo <= i && i <= ihi) && (jlo <= j && j <= jhi)) {
                    // Skip. We are in the area that is to be spared.
                }
                else {
                    X(i,j) = alpha * X(i,j);
                }
            }
        }
    }

#undef X
}


/// Convert int scaling factor to double-precision scaling, if necessary.
double convert_scaling(scaling_t alpha)
{
#ifdef INTSCALING
    double scaling = ldexp(1.0, alpha);
#else
    double scaling = alpha;
#endif

    return scaling;
}


void init_scaling_factor(int n, scaling_t *restrict const alpha)
{
#ifdef INTSCALING
    for (int i = 0; i < n; i++)
        alpha[i] = 0;
#else
    for (int i = 0; i < n; i++)
        alpha[i] = 1.0;
#endif
}


/**
 * @brief Compute common scaling factor alpha_min / alpha.
 */
double compute_upscaling(scaling_t alpha_min, scaling_t alpha)
{
    double scaling;

#ifdef INTSCALING
    // Common scaling is 2^alpha_min / 2^alpha.
    scaling_t exp = alpha_min - alpha;
    scaling = ldexp(1.0, exp);
#else
    scaling = alpha_min / alpha;
#endif

    return scaling;
}


/**
 * @brief Compute common scaling factor (alpha_min / alpha) * beta.
 */
double compute_combined_upscaling(
    scaling_t alpha_min, scaling_t alpha, scaling_t beta)
{
    double scaling;

#ifdef INTSCALING
    // Common scaling is (2^alpha_min / 2^alpha) * 2^beta.
    scaling_t exp = alpha_min - alpha + beta;
    scaling = ldexp(1.0, exp);
#else
    scaling = (alpha_min / alpha) * beta;
#endif

    return scaling;
}


void update_global_scaling(scaling_t *global, scaling_t phi)
{
#ifdef INTSCALING
    *global = phi + (*global);
#else
    *global = phi * (*global);
#endif
}


void update_norm(double *norm, scaling_t phi)
{
#ifdef INTSCALING
    *norm = ldexp(1.0, phi) * (*norm);
#else
    *norm = phi * (*norm);
#endif
}


void dgemm(
    const char transa, const char transb,
    const int m, const int n, const int k,
    const double alpha, const double *restrict const A, const int ldA,
    const double *restrict const B, const int ldB,
    const double beta, double *restrict const C, const int ldC)
{
    extern void dgemm_(
        const char *transa, const char *transb,
        const int *m, const int *n, const int *k,
        const double *alpha, const double *a, const int *lda,
        const double *b, const int *ldb,
        const double *beta, double *c, const int *ldc);

    dgemm_(&transa, &transb,
           &m, &n, &k,
           &alpha, A, &ldA,
           B, &ldB,
           &beta, C, &ldC);

}


double dlange(const char norm, const int m, const int n,
    const double *restrict const A, const int ldA)
{
    double work;
    extern double dlange_(const char *norm, const int *m, const int *n,
        const double *a, const int *ldA, double *work);

    double res = dlange_(&norm, &m, &n, A, &ldA, &work);

    return res;
}

