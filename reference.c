#include "reference.h"
#include "partition.h"
#include "utils.h"
#include "defines.h"
#include "timing.h"
#include "validation.h"
#include "norm.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>



void dtrsyl(
    const char transA, const char transB, const double sgn,
    const int m, const int n,
    const double *restrict const A, int ldA,
    const double *restrict const B, int ldB,
    double *restrict const C, int ldC,
    double *restrict const scale)
{
    extern void dtrsyl_(
        const char *trana, const char *tranb, const int *isgn,
        const int *m, const int *n,
        const double *A, const int *ldA,
        const double *B, const int *ldB,
        double *C, const int *ldC, double *scale,
        int *info);

    assert(sgn == 1.0 || sgn == -1.0);

    int info;
    int isgn;
    if (sgn == 1.0)
        isgn = 1;
    else
        isgn = -1;

    dtrsyl_(&transA, &transB, &isgn,
        &m, &n,
        A, &ldA,
        B, &ldB,
        C, &ldC, scale,
        &info);

    if (info == 1) {
        printf("WARNING: A and B have common or very close eigenvalues.\n");
    }
    else if (info < 0) {
        printf("ERROR: The %d-th argument to DTRSYL was invalid\n", -info);
    }
}



void solve_sylvester_dtrsyl(
    const int transA, const int transB,
    double sgn,
    double ***A_tiles, int ldA,
    double ***B_tiles, int ldB,
    double ***C_tiles, int ldC,
    partitioning_t *p)
{
    int num_tile_rows = p->num_blk_rows;
    int num_tile_cols = p->num_blk_cols;
    int *first_row = p->first_row;
    int *first_col = p->first_col;

    // Extract matrix dimensions.
    const int m = first_row[num_tile_rows];
    const int n = first_col[num_tile_cols];

    double tm_start, tm_end;

    double scale = 1.0;

    // Duplicate C.
    double *X = (double *) malloc((size_t)ldC * n * sizeof(double));
    memcpy(X, C_tiles[0][0], (size_t)ldC * n * sizeof(double));

    // Compute solution with LAPACK.
    tm_start = get_time();
    dtrsyl(transA, transB, sgn, m, n,
        A_tiles[0][0], ldA, B_tiles[0][0], ldB,
        X, ldC, &scale);
    tm_end = get_time();

#if INTSCALING
    validate(sgn, m, n, A_tiles[0][0], ldA, B_tiles[0][0], ldB,
        C_tiles[0][0], ldC, X, ldC, ilogb(scale));

#else
    validate(sgn, m, n, A_tiles[0][0], ldA, B_tiles[0][0], ldB,
        C_tiles[0][0], ldC, X, ldC, scale);
#endif
    printf("...scale = %.6e\n", scale);

    // Clean up.
    free(X);


    printf("...LAPACK time = %.2f s.\n", tm_end - tm_start);
}

