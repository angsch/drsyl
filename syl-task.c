#include "syl-task.h"
#include "majorant.h"
#include "robust.h"
#include "utils.h"
#include "norm.h"
#include "defines.h"
#include "partition.h"

#include <assert.h>
#include <math.h>
#include <mm_malloc.h>



// ||A^T||_oo
static double bound_row_vector(int n, const double *A, int ldA)
{
    double ub = 0;
    for (int i = 0; i < n; i++)
        ub = fmax(ub, fabs(A[i * ldA]));

    return ub;
}



void syl(
    double sgn,
    int m, int n,
    const double *restrict A, int ldA, const double anorm,
    const double *restrict B, int ldB, const double bnorm,
    double *restrict C, int ldC, double *restrict cnorm,
    scaling_t *restrict scale,
    const double smin)
{
    // Validate inputs.
    assert(sgn == 1.0 || sgn == -1.0);
    assert(ldA >= m);
    assert(ldC >= m);
    assert(ldB >= n);

    if (m == 0 || n == 0)
        return;

    // Init scaling factor.
    init_scaling_factor(1, scale);

#define A(k,i) A[(k) + (i) * ldA]
#define B(j,l) B[(j) + (l) * ldB]
#define C(k,l) C[(k) + (l) * ldC]

    // Solve A * X + sgn * X * B = scale * C.
    //
    // The (k,l)-th block of X is determined starting from the bottom-left
    // corner column by column by
    //
    // A(k,k) * X(k,l) + sgn * X(k,l) * B(l,l) = C(k,l) - R(k,l)
    //
    // where
    //          m-1                           l-1
    // R(k,l) = SUM [A(k,i) * X(i,l)] + sgn * SUM [X(k,j) * B(j,l)].
    //         i=k+1                          j=0

    // l: column index
    // k: row index

    int Bll_is_1x1, Akk_is_1x1;

    // Loop over columns from left to right.
    for (int l = 0; l < n; l++) {
        // Check if the next block B(l,l) is 1-by-1 or 2-by-2.
        if (l < n - 1 && B(l + 1, l) != 0.0) {
            // 2-by-2 block spans B(l:l+1,l:l+1).
            Bll_is_1x1 = 0;
        }
        else {
            // 1-by-1 block spans B(l,l).
            Bll_is_1x1 = 1;
        }


        // Loop over rows bottom up.
        for (int k = m - 1; k >= 0; k--) {
            // Next block A(k,k) is 1-by-1 or 2-by-2.
            if (k > 0 && A(k, k - 1) != 0.0) {
                // 2-by-2 block spans A(k-1:k,k-1:k).
                Akk_is_1x1 = 0;
            }
            else {
                // 1-by-1 block spans A(k,k).
                Akk_is_1x1 = 1;
            }


            if (Bll_is_1x1 && Akk_is_1x1) {
                ////////////////////////////////////////////////////////////////
                // Solve C(k,l) = C(k,l) / (A(k,k) + sgn * B(l,l)).
                ////////////////////////////////////////////////////////////////
                scaling_t phi;
                init_scaling_factor(1, &phi);
                solve_a1x1_b1x1(A(k,k), sgn, B(l,l), smin, &C(k,l), &phi);
                update_global_scaling(scale, phi);

                // Scale remaining tile.
                scale_excluding(m, n, k, k, l, l, C, ldC, &phi);


                ////////////////////////////////////////////////////////////////
                // Column update to the top.
                ////////////////////////////////////////////////////////////////
                double colbound = vector_infnorm(k-1, &C(0,l));
                phi = protect_update(anorm, fabs(C(k,l)), colbound);
                update_global_scaling(scale, phi);

                // Apply the scaling to the whole tile.
                scale_tile(m, n, C, ldC, &phi);

                // Now it is safe to execute the linear update.
                for (int i = 0; i < k; i++)
                    C(i,l) = C(i,l) - A(i,k) * C(k,l);


                ////////////////////////////////////////////////////////////////
                // Row update to the right.
                ////////////////////////////////////////////////////////////////
                double rowbound = bound_row_vector(n-(l+1), &C(k,l+1), ldC);
                phi = protect_update(bnorm, fabs(C(k,l)), rowbound);
                update_global_scaling(scale, phi);

                // Apply the scaling to the whole tile.
                scale_tile(m, n, C, ldC, &phi);

                // Now it is safe to execute the linear update.
                for (int j = l + 1; j < n; j++)
                    C(k,j) = C(k,j) - sgn * C(k,l) * B(l,j);

            }
            else if (Bll_is_1x1 && !Akk_is_1x1) {
                // B(l,l) is 1-by-1, A(k:k-1,k:k-1) is 2-by-2.

                ////////////////////////////////////////////////////////////////
                // Solve.
                ////////////////////////////////////////////////////////////////
                scaling_t phi;
                init_scaling_factor(1, &phi);
                solve_a2x2_b1x1(&A(k-1,k-1), ldA, sgn, B(l,l),
                    smin, &C(k-1,l), &phi);
                update_global_scaling(scale, phi);

                // Scale remaining tile.
                scale_excluding(m, n, k - 1, k, l, l, C, ldC, &phi);


                ////////////////////////////////////////////////////////////////
                // Column update to the top.
                ////////////////////////////////////////////////////////////////
                //for (int i = 0; i < k - 1; i++) {
                //    C(i,l) = C(i,l) - A(i,k-1) * C(k-1,l);
                //    C(i,l) = C(i,l) - A(i,k) * C(k,l);
                //}
                // Determine bounds for the first linear update - note that the
                // 2 linear updates cannot be joined.
                double colbound = vector_infnorm(k-2, &C(0,l));
                phi = protect_update(anorm, fabs(C(k-1,l)), colbound);
                update_global_scaling(scale, phi);

                // Apply the scaling to the whole tile.
                scale_tile(m, n, C, ldC, &phi);

                // Now it is safe to execute the first linear update.
                for (int i = 0; i < k - 1; i++)
                    C(i,l) = C(i,l) - A(i,k-1) * C(k-1,l);

                // Determine bounds for the second linear update.
                colbound = vector_infnorm(k-2, &C(0,l));
                phi = protect_update(anorm, fabs(C(k,l)), colbound);
                update_global_scaling(scale, phi);

                // Apply the scaling to the whole tile.
                scale_tile(m, n, C, ldC, &phi);

                // Now it is safe to execute the second linear update.
                for (int i = 0; i < k - 1; i++)
                    C(i,l) = C(i,l) - A(i,k) * C(k,l);


                ////////////////////////////////////////////////////////////////
                // Row update to the right.
                ////////////////////////////////////////////////////////////////
                // Determine bounds so that the 2 linear updates can be joined.
                double rowbound = maxf(bound_row_vector(n-l-1, &C(k-1,l+1), ldC),
                                       bound_row_vector(n-l-1, &C(k,l+1), ldC));
                double lhsbound = maxf(fabs(C(k-1,l)), fabs(C(k,l)));
                phi = protect_update(bnorm, lhsbound, rowbound);
                update_global_scaling(scale, phi);

                // Apply the scaling to the whole tile.
                scale_tile(m, n, C, ldC, &phi);

                // Now it is safe to execute the linear update.
                for (int j = l + 1; j < n; j++) {
                    C(k-1,j) = C(k-1,j) - sgn * C(k-1,l) * B(l,j);
                    C(k,j) = C(k,j) - sgn * C(k,l) * B(l,j);
                }
            }
            else if (!Bll_is_1x1 && Akk_is_1x1) {
                // B(l:l+1,l:l+1) is 2-by-2, A(k,k) is 1-by-1.

                scaling_t phi;
                init_scaling_factor(1, &phi);
                solve_a1x1_b2x2(A(k,k), sgn, &B(l,l), ldB,
                    smin, &C(k,l), ldC, &phi);
                update_global_scaling(scale, phi);

                // Scale remaining tile.
                scale_excluding(m, n, k, k, l, l + 1, C, ldC, &phi);

                // Column update to the top.
                // Determine bounds so that the 2 linear updates can be joined.
                double colbound = maxf(vector_infnorm(k-1, &C(0,l)),
                                       vector_infnorm(k-1, &C(0,l+1)));
                double rhsbound = maxf(fabs(C(k,l)), fabs(C(k,l+1)));
                phi = protect_update(anorm, rhsbound, colbound);
                update_global_scaling(scale, phi);

                // Apply the scaling to the whole tile.
                scale_tile(m, n, C, ldC, &phi);

                // Now it is safe to execute the linear update jointly.
                for (int i = 0; i < k; i++) {
                    C(i,l) = C(i,l) - A(i,k) * C(k,l);
                    C(i,l+1) = C(i,l+1) - A(i,k) * C(k,l+1);
                }


                // Row update to the right.
                //for (int j = l + 2; j < n; j++) {
                //    C(k,j) = C(k,j) - sgn * C(k,l) * B(l,j);
                //    C(k,j) = C(k,j) - sgn * C(k,l+1) * B(l+1,j);
                //}
                // Determine bounds for the first linear update - note that the
                // 2 linear updates cannot be joined.
                double rowbound = bound_row_vector(n-(l+2), &C(k,l+2), ldC);
                phi = protect_update(bnorm, fabs(C(k,l)), rowbound);
                update_global_scaling(scale, phi);

                // Apply the scaling to the whole tile.
                scale_tile(m, n, C, ldC, &phi);

                // Now it is safe to execute the first linear update.
                for (int j = l + 2; j < n; j++)
                    C(k,j) = C(k,j) - sgn * C(k,l) * B(l,j);

                // Determine bounds for the second linear update.
                // recomputation is necessary because first linear update overwrites this
                rowbound = bound_row_vector(n-(l+2), &C(k,l+2), ldC);
                phi = protect_update(bnorm, fabs(C(k,l+1)), rowbound);
                update_global_scaling(scale, phi);

                // Apply the scaling to the whole tile.
                scale_tile(m, n, C, ldC, &phi);

                // Now it is safe to execute the second linear update.
                for (int j = l + 2; j < n; j++)
                    C(k,j) = C(k,j) - sgn * C(k,l+1) * B(l+1,j);

            }
            else {
                // B(l:l+1,l:l+1) and A(k-1:k,k-1:k) are 2-by-2.

                ////////////////////////////////////////////////////////////////
                // Solve.
                ////////////////////////////////////////////////////////////////
                // DLASY2
                scaling_t phi;
                init_scaling_factor(1, &phi);
                solve_a2x2_b2x2(&A(k-1,k-1), ldA, sgn, &B(l,l), ldB,
                    &C(k-1,l), ldC, &phi);
                update_global_scaling(scale, phi);

                // Scale remaining tile.
                scale_excluding(m, n, k-1, k, l, l+1, C, ldC, &phi);

                ////////////////////////////////////////////////////////////////
                // Column update to the top.
                ////////////////////////////////////////////////////////////////
                // Split the 4 linear updates
                //for (int i = 0; i < k - 1; i++) {
                //    C(i,l) = C(i,l) - A(i,k-1) * C(k-1,l);            // (1)
                //    C(i,l+1) = C(i,l+1) - A(i,k-1) * C(k-1,l+1);      // (2)
                //    C(i,l) = C(i,l) - A(i,k) * C(k,l);                // (3)
                //    C(i,l+1) = C(i,l+1) - A(i,k) * C(k,l+1);          // (4)
                //}
                // into two jointly executed linear updates (1)+(2) and (3)+(4).
                // Determine bounds for (1)+(2).
                double colbound = maxf(vector_infnorm(k-2, &C(0,l)),
                                       vector_infnorm(k-2, &C(0,l+1)));
                double rhsbound = maxf(fabs(C(k-1,l)), fabs(C(k-1,l+1)));
                phi = protect_update(anorm, rhsbound, colbound);
                update_global_scaling(scale, phi);

                // Apply the scaling to the whole tile.
                scale_tile(m, n, C, ldC, &phi);

                // Now it is safe to execute the linear updates (1)+(2).
                for (int i = 0; i < k - 1; i++) {
                    C(i,l) = C(i,l) - A(i,k-1) * C(k-1,l);
                    C(i,l+1) = C(i,l+1) - A(i,k-1) * C(k-1,l+1);
                }

                // Determine bounds for (3)+(4).
                colbound = maxf(vector_infnorm(k-2, &C(0,l)),
                                vector_infnorm(k-2, &C(0,l+1)));
                rhsbound = maxf(fabs(C(k,l)), fabs(C(k,l+1)));
                phi = protect_update(anorm, rhsbound, colbound);
                update_global_scaling(scale, phi);

                // Apply the scaling to the whole tile.
                scale_tile(m, n, C, ldC, &phi);

                // Now it is safe to execute the linear updates (3)+(4).
                for (int i = 0; i < k - 1; i++) {
                    C(i,l) = C(i,l) - A(i,k) * C(k,l);
                    C(i,l+1) = C(i,l+1) - A(i,k) * C(k,l+1);
                }


                ////////////////////////////////////////////////////////////////
                // Row update to the right.
                ////////////////////////////////////////////////////////////////
                // Split the 4 linear updates
                //for (int j = l + 2; j < n; j++) {
                //    C(k-1,j) = C(k-1,j) - sgn * C(k-1,l) * B(l,j);     // (1)
                //    C(k,j) = C(k,j) - sgn * C(k,l) * B(l,j);           // (2)
                //    C(k-1,j) = C(k-1,j) - sgn * C(k-1,l+1) * B(l+1,j); // (3)
                //    C(k,j) = C(k,j) - sgn * C(k,l+1) * B(l+1,j);       // (4)
                //}
                // into two jointly executed linear updates (1)+(2) and (3)+(4).
                // Determine bounds for (1)+(2).
                double rowbound = maxf(bound_row_vector(n-l-2, &C(k-1,l+2), ldC),
                                       bound_row_vector(n-l-2, &C(k,l+2), ldC));
                double lhsbound = maxf(fabs(C(k-1,l)),fabs(C(k,l)));
                phi = protect_update(bnorm, lhsbound, rowbound);
                update_global_scaling(scale, phi);

                // Apply the scaling to the whole tile.
                scale_tile(m, n, C, ldC, &phi);

                // Now it is safe to execute the linear updates (1)+(2).
                for (int j = l + 2; j < n; j++) {
                    C(k-1,j) = C(k-1,j) - sgn * C(k-1,l) * B(l,j);
                    C(k,j) = C(k,j) - sgn * C(k,l) * B(l,j);
                }

                // Determine bounds for (3)+(4).
                rowbound = maxf(bound_row_vector(n-l-2, &C(k-1,l+2), ldC),
                                bound_row_vector(n-l-2, &C(k,l+2), ldC));
                lhsbound = maxf(fabs(C(k-1,l+1)), fabs(C(k,l+1)));
                phi = protect_update(bnorm, lhsbound, rowbound);
                update_global_scaling(scale, phi);

                // Apply the scaling to the whole tile.
                scale_tile(m, n, C, ldC, &phi);

                // Now it is safe to execute the linear update (3)+(4).
                for (int j = l + 2; j < n; j++) {
                    C(k-1,j) = C(k-1,j) - sgn * C(k-1,l+1) * B(l+1,j);
                    C(k,j) = C(k,j) - sgn * C(k,l+1) * B(l+1,j);
                }

            }

            if (!Akk_is_1x1) {
                // We processed a 2-by-2 block in A, so skip the next row.
                k--;
            }
        }

        if (!Bll_is_1x1) {
            // We processed a 2-by-2 block in B, so skip the next column.
            l++;
        }
    }

    // Compute upper bound of C.
    *cnorm = matrix_infnorm(m, n, C, ldC);

#undef A
#undef B
#undef C
}


// This is a light-weight version of the linear update
// C := C - T * X.
// (m x n)  (m x k) (k x n)
static void update_top(
    int m, int n, int k,
    const double *restrict const T, int ldT, const double tnorm,
    double *restrict const X, int ldX, double *restrict const xnorm,
    scaling_t *restrict const xscale,
    double *restrict const C, int ldC, double *restrict const cnorm,
    scaling_t *restrict const cscale)
{
    // Scaling of X, C.
    scaling_t xscaling = *xscale;
    scaling_t cscaling = *cscale;

    // Local scaling factor.
    scaling_t zeta;

    // Status flag if C or X have to be rescaled.
    int rescale_C = 0;
    int rescale_X = 0;

    // Bound right-hand side C.
    *cnorm = matrix_infnorm(m, n, C, ldC);

    ////////////////////////////////////////////////////////////////////////////
    // Compute scaling factor.
    ////////////////////////////////////////////////////////////////////////////

    // Simulate consistent scaling.
    if (cscaling < xscaling) {
        // The common scaling factor is cscale.
        const double s = compute_upscaling(cscaling, xscaling);

        // Mark X for scaling. Physical rescaling is deferred.
        rescale_X = 1;

        // Update norm.
        *xnorm = s * (*xnorm);
    }
    else if (xscaling < cscaling) {
        // The common scaling factor is tscale.
        const double s = compute_upscaling(xscaling, cscaling);

        // Mark C for scaling. Physical rescaling is deferred.
        rescale_C = 1;

        // Update norm.
        *cnorm = s * (*cnorm);
    }

    // Compute scaling factor needed to survive the linear update.
    zeta = protect_update(tnorm, *xnorm, *cnorm);

#ifdef INTSCALING
    if (zeta != 0) {
        rescale_X = 1;
        rescale_C = 1;
    }
#else
    if (zeta != 1.0) {
        rescale_X = 1;
        rescale_C = 1;
    }
#endif

    // Physically rescale X, if necessary.
    if (rescale_X) { // k x n
        if (cscaling < xscaling) {
            // The common scaling factor is cscaling. Combine with robust update
            // scaling factor.
            const double s = compute_combined_upscaling(cscaling, xscaling, zeta);
            for (int j = 0; j < n; j++)
                for (int i = 0; i < k; i++)
                    X[i + ldX * j] = s * X[i + ldX * j];
        }
        else if (xscaling < cscaling) {
            // The common scaling factor is xscaling. Rescale with
            // robust update factor.
            const double s = convert_scaling(zeta);
            for (int j = 0; j < n; j++)
                for (int i = 0; i < k; i++)
                    X[i + ldX * j] = s * X[i + ldX * j];
        }
        else {
            // X and C are consistently scaled.

            // Rescale with robust update factor.
            const double s = convert_scaling(zeta);
            for (int j = 0; j < n; j++)
                for (int i = 0; i < k; i++)
                    X[i + ldX * j] = s * X[i + ldX * j];
        }
    }


    // Physically rescale C, if necessary.
    if (rescale_C) { // m x n
        if (cscaling < xscaling) {
            // The common scaling factor is cscaling. Rescale with
            // robust update factor.
            const double s = convert_scaling(zeta);
            for (int j = 0; j < n; j++)
                for (int i = 0; i < m; i++)
                    C[i + j * ldC] = s * C[i + j * ldC];
        }
        else if (xscaling < cscaling) {
            // The common scaling factor is xscaling. Combine with robust
            // scaling factor.
            const double s = compute_combined_upscaling(xscaling, cscaling, zeta);
            for (int j = 0; j < n; j++)
                for (int i = 0; i < m; i++)
                    C[i + j * ldC] = s * C[i + j * ldC];
        }
        else {
            // X and C are consistently scaled.

            // Rescale with robust update factor.
            const double s = convert_scaling(zeta);
            for (int j = 0; j < n; j++)
                for (int i = 0; i < m; i++)
                    C[i + j * ldC] = s * C[i + j * ldC];
        }
    }

    // Update global scaling of C.
#ifdef INTSCALING
    *cscale = min(cscaling, xscaling) + zeta;
#else
    *cscale = minf(cscaling, xscaling) * zeta;
#endif


    ////////////////////////////////////////////////////////////////////////////
    // Compute update.
    ////////////////////////////////////////////////////////////////////////////

    //  C := C - T  * X.
    //(mxn)    (mxk)(kxn)
    dgemm('N', 'N',
          m, n, k,
          -1.0, T, ldT,
          X, ldX,
          1.0, C, ldC);


    ////////////////////////////////////////////////////////////////////////////
    // Recompute norm.
    ////////////////////////////////////////////////////////////////////////////
    *cnorm = matrix_infnorm(m, n, C, ldC);
}



// This is a light-weight version of the linear update
// C := C - sgn * X * B.
// (m x n)  (m x k) (k x n)
static void update_right(
    int m, int n, int k, const double sgn,
    double *restrict const X, int ldX, double *restrict const xnorm,
    scaling_t *restrict const xscale,
    const double *restrict const B, int ldB, const double bnorm,
    double *restrict const C, int ldC, double *restrict const cnorm,
    scaling_t *restrict const cscale)
{
    // Scaling of X, C.
    scaling_t xscaling = *xscale;
    scaling_t cscaling = *cscale;

    // Local scaling factor.
    scaling_t zeta;

    // Status flag if C or X have to be rescaled.
    int rescale_C = 0;
    int rescale_X = 0;

    // Bound right-hand side C.
    *cnorm = matrix_infnorm(m, n, C, ldC);


    ////////////////////////////////////////////////////////////////////////////
    // Compute scaling factor.
    ////////////////////////////////////////////////////////////////////////////

    // Simulate consistent scaling.
    if (cscaling < xscaling) {
        // The common scaling factor is cscale.
        const double s = compute_upscaling(cscaling, xscaling);

        // Mark X for scaling. Physical rescaling is deferred.
        rescale_X = 1;

        // Update norm.
        *xnorm = s * (*xnorm);
    }
    else if (xscaling < cscaling) {
        // The common scaling factor is tscale.
        const double s = compute_upscaling(xscaling, cscaling);

        // Mark C for scaling. Physical rescaling is deferred.
        rescale_C = 1;

        // Update norm.
        *cnorm = s * (*cnorm);
    }

    // Compute scaling factor needed to survive the linear update.
    zeta = protect_update(*xnorm, bnorm, *cnorm);

#ifdef INTSCALING
    if (zeta != 0) {
        rescale_X = 1;
        rescale_C = 1;
    }
#else
    if (zeta != 1.0) {
        rescale_X = 1;
        rescale_C = 1;
    }
#endif

    // Physically rescale X, if necessary.
    if (rescale_X) { // (m x k)
        if (cscaling < xscaling) {
            // The common scaling factor is cscaling. Combine with robust update
            // scaling factor.
            const double s = compute_combined_upscaling(cscaling, xscaling, zeta);
            for (int j = 0; j < k; j++)
                for (int i = 0; i < m; i++)
                    X[i + ldX * j] = s * X[i + ldX * j];
        }
        else if (xscaling < cscaling) {
            // The common scaling factor is xscaling. Rescale with
            // robust update factor.
            const double s = convert_scaling(zeta);
            for (int j = 0; j < k; j++)
                for (int i = 0; i < m; i++)
                    X[i + ldX * j] = s * X[i + ldX * j];
        }
        else {
            // X and C are consistently scaled.

            // Rescale with robust update factor.
            const double s = convert_scaling(zeta);
            for (int j = 0; j < k; j++)
                for (int i = 0; i < m; i++)
                    X[i + ldX * j] = s * X[i + ldX * j];
        }
    }


    // Physically rescale C, if necessary.
    if (rescale_C) { // m x n
        if (cscaling < xscaling) {
            // The common scaling factor is cscaling. Rescale with
            // robust update factor.
            const double s = convert_scaling(zeta);
            for (int j = 0; j < n; j++)
                for (int i = 0; i < m; i++)
                    C[i + j * ldC] = s * C[i + j * ldC];
        }
        else if (xscaling < cscaling) {
            // The common scaling factor is xscaling. Combine with robust
            // scaling factor.
            const double s = compute_combined_upscaling(xscaling, cscaling, zeta);
            for (int j = 0; j < n; j++)
                for (int i = 0; i < m; i++)
                    C[i + j * ldC] = s * C[i + j * ldC];
        }
        else {
            // X and C are consistently scaled.

            // Rescale with robust update factor.
            const double s = convert_scaling(zeta);
            for (int j = 0; j < n; j++)
                for (int i = 0; i < m; i++)
                    C[i + j * ldC] = s * C[i + j * ldC];
        }
    }

    // Update global scaling of C.
#ifdef INTSCALING
    *cscale = min(cscaling, xscaling) + zeta;
#else
    *cscale = minf(cscaling, xscaling) * zeta;
#endif

    ////////////////////////////////////////////////////////////////////////////
    // Compute update.
    ////////////////////////////////////////////////////////////////////////////

    //  C := C - sgn * X  * B.
    //(mxn)          (mxk)(kxn)
    dgemm('N', 'N',
          m, n, k,
          -sgn, X, ldX,
          B, ldB,
          1.0, C, ldC);


    ////////////////////////////////////////////////////////////////////////////
    // Recompute norm.
    ////////////////////////////////////////////////////////////////////////////
    *cnorm = matrix_infnorm(m, n, C, ldC);
}



void blocked_syl(
    double sgn,
    int m, int n,
    const double *restrict A, int ldA, const double anorm,
    const double *restrict B, int ldB, const double bnorm,
    double *restrict C, int ldC, double *restrict cnorm,
    scaling_t *restrict scale,
    const double smin)
{
    // Validate inputs.
    assert(sgn == 1.0 || sgn == -1.0);
    assert(ldA >= m);
    assert(ldC >= m);
    assert(ldB >= n);

    if (m <= 0 || n <= 0)
        return;

#define A(k,i) A[(k) + (i) * ldA]
#define B(j,l) B[(j) + (l) * ldB]
#define C(k,l) C[(k) + (l) * ldC]


    ////////////////////////////////////////////////////////////////////////////
    // Compute internal blocking.
    ////////////////////////////////////////////////////////////////////////////

    // Partition C internally even further. This partitioning induces a
    // finer partitioning of A and B.
    int blksz = 24;
    int num_blk_rows = (m + blksz - 1) / blksz;
    int num_blk_cols = (n + blksz - 1) / blksz;
    int first_row[num_blk_rows + 1];
    int first_col[num_blk_cols + 1];


    // Extract the lambda type of A.
    int lambda_type_A[m];
    lambda_type_A[m - 1] = REAL;
    for (int i = 0; i < m - 1; i++) {
        if (A(i+1,i) != 0.0) {
            lambda_type_A[i] = CMPLX;
            lambda_type_A[i + 1] = CMPLX;
            i++;
        }
        else {
            lambda_type_A[i] = REAL;
        }
    }


    // Extract the lambda type of B.
    int lambda_type_B[n];
    lambda_type_B[n - 1] = REAL;
    for (int i = 0; i < n; i++) {
        if (B(i+1,i) != 0.0) {
            lambda_type_B[i] = CMPLX;
            lambda_type_B[i + 1] = CMPLX;
            i++;
        }
        else {
            lambda_type_B[i] = REAL;
        }
    }

    partition(m, lambda_type_A, num_blk_rows, blksz, first_row);
    partition(n, lambda_type_B, num_blk_cols, blksz, first_col);

    // Prepare internal scaling factors.
    double Cnorms_internal[num_blk_rows * num_blk_cols];
    scaling_t scales_internal[num_blk_rows * num_blk_cols];
    init_scaling_factor(num_blk_rows * num_blk_cols, scales_internal);
#define scales_internal(row,col) scales_internal[(row) + (col) * num_blk_rows]
#define Cnorms_internal(row,col) Cnorms_internal[(row) + (col) * num_blk_rows]



    ////////////////////////////////////////////////////////////////////////////
    // Solve Sylvester equation.
    ////////////////////////////////////////////////////////////////////////////

    for (int k = num_blk_rows - 1; k >= 0; k--) {
        for (int l = 0; l < num_blk_cols; l++) {
            // Solve A(k,k) * X(k,l) + sgn * X(k,l) * B(l,l) = C(k,l).
            {
                // Compute dimensions of C(k,l).
                int num_rows = first_row[k + 1] - first_row[k];
                int num_cols = first_col[l + 1] - first_col[l];

                syl(sgn, num_rows, num_cols,
                    &A(first_row[k], first_row[k]), ldA, anorm,
                    &B(first_col[l], first_col[l]), ldB, bnorm,
                    &C(first_row[k], first_col[l]), ldC, &Cnorms_internal(k,l),
                    &scales_internal(k,l), smin);
            }

            // Update tiles to the top.
            for (int i = k - 1; i >= 0; i--) {
                // Dimensions of C(i,l).
                const int num_rows = first_row[i + 1] - first_row[i];
                const int num_cols = first_col[l + 1] - first_col[l];

                // Number of columns in A(i,k)/rows in C(k,l).
                int num_inner = first_row[k + 1] - first_row[k];

                // C(i,l) := C(i,l) - A(i,k) * C(k,l).
                update_top(num_rows, num_cols, num_inner,
                    &A(first_row[i], first_row[k]), ldA, anorm,
                    &C(first_row[k], first_col[l]), ldC, &Cnorms_internal(k,l), &scales_internal(k,l),
                    &C(first_row[i], first_col[l]), ldC, &Cnorms_internal(i,l), &scales_internal(i,l));
            }

            // Update tiles to the right.
            for (int j = l + 1; j < num_blk_cols; j++) {
                // Dimensions of C(k,j).
                const int num_rows = first_row[k + 1] - first_row[k];
                const int num_cols = first_col[j + 1] - first_col[j];

                // Number of columns in C(k,l)/rows in B(l,j).
                const int num_inner = first_col[l + 1] - first_col[l];

                // C(k,j) := C(k,j) - sgn * C(k,l) * B(l,j).
                update_right(num_rows, num_cols, num_inner, sgn,
                    &C(first_row[k], first_col[l]), ldC, &Cnorms_internal(k,l), &scales_internal(k,l),
                    &B(first_col[l], first_col[j]), ldB, bnorm,
                    &C(first_row[k], first_col[j]), ldC, &Cnorms_internal(k,j), &scales_internal(k,j));
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Consolidate scaling factors and scale consistently.
    ////////////////////////////////////////////////////////////////////////////

    scaling_t scale_min = min_element(num_blk_rows * num_blk_cols, scales_internal);
    for (int k = num_blk_rows - 1; k >= 0; k--) {
        for (int l = 0; l < num_blk_cols; l++) {
            // Compute dimensions of C(k,l).
            int num_rows = first_row[k + 1] - first_row[k];
            int num_cols = first_col[l + 1] - first_col[l];

            scaling_t ratio;
            #ifdef INTSCALING
                ratio = scale_min - scales_internal(k,l);
            #else
                ratio = scale_min / scales_internal(k,l);
            #endif

            scale_tile(num_rows, num_cols, &C(first_row[k], first_col[l]), ldC, &ratio);
        }
    }

    update_global_scaling(scale, scale_min);

    // Compute upper bound of C.
    *cnorm = matrix_infnorm(m, n, C, ldC);


#undef A
#undef B
#undef C
}


