#include "utils.h"
#include "defines.h"
#include "norm.h"
#include "update-task.h"
#include "robust.h"

#include <stdlib.h>
#include <stdio.h>
#include <mm_malloc.h>
#include <string.h>
#include <math.h>



void update(
    int m, int n, int k,
    omp_lock_t *lock,
    const double alpha, double *restrict const Ain, int ldAin,
    const double ainnorm, const scaling_t ainscale,
    double *restrict const B, int ldB, const double bnorm,
    double *restrict const C, int ldC, double *restrict const cnorm, 
    scaling_t *restrict const cscale)
{
#ifndef NDEBUG
    printf("update (%dx%d)(%dx%d)\n", m, k, k, n);
    printf("update thread id = %d\n", omp_get_thread_num());
#endif

    // Scaling of C.
    scaling_t cscaling = *cscale;

    // Local scaling factor.
    scaling_t zeta;

    // Pointer to A - either a copy or the original memory.
    double *A;

    // Local norm of A.
    double anorm;

    // Status flag if A or C have to be rescaled.
    int rescale_A = 0;
    int rescale_C = 0;


    ////////////////////////////////////////////////////////////////////////////
    // Compute scaling factor.
    ////////////////////////////////////////////////////////////////////////////

    while (!omp_test_lock(lock)) {
        #pragma omp taskyield
        ;
    }

    // Copy norm.
    anorm = ainnorm;

    // Bound right-hand side C.
    *cnorm = matrix_infnorm(m, n, C, ldC);

    // Simulate consistent scaling.
    if (cscaling < ainscale) {
        // The common scaling factor is cscale.
        const double s = compute_upscaling(cscaling, ainscale);

        // Mark A for scaling. Physical rescaling is deferred.
        rescale_A = 1;

        // Update norm.
        anorm = s * ainnorm;
    }
    else if (ainscale < cscaling) {
        // The common scaling factor is ascale.
        const double s = compute_upscaling(ainscale, cscaling);

        // Mark C for scaling. Physical rescaling is deferred.
        rescale_C = 1;

        // Update norm.
        *cnorm = s * (*cnorm);
    }
    else {
        // Nothing to do. C and A are consistently scaled.
    }

    // Compute scaling factor needed to survive the linear update.
    zeta = protect_update(anorm, bnorm, *cnorm);

#ifdef INTSCALING
    if (zeta != 0) {
        rescale_A = 1;
        rescale_C = 1;
    }
#else
    if (zeta != 1.0) {
        rescale_A = 1;
        rescale_C = 1;
    }
#endif

    // If A has to be rescaled, take a copy and do the scaling on the copy.
    if (rescale_A) {
        A = (double *) _mm_malloc((size_t)ldAin * k * sizeof(double), ALIGNMENT);
        if (cscaling < ainscale) {
            // Copy A and simultaneously rescale.
            double s = compute_combined_upscaling(cscaling, ainscale, zeta);
            for (int j = 0; j < k; j++)
                for (int i = 0; i < m; i++)
                    A[i + ldAin * j] = s * Ain[i + ldAin * j];
        }
        else if (ainscale < cscaling) {
            // Copy A and simultaneously rescale with robust update factor.
            double s = convert_scaling(zeta);
            for (int j = 0; j < k; j++)
                for (int i = 0; i < m; i++)
                    A[i + ldAin * j] = s * Ain[i + ldAin * j];
        }
    }
    else {
        // A does not have to be rescaled. Operate on original memory.
        A = Ain;
    }

    // If C has to be rescaled, directly modify C.
    if (rescale_C) {
        if (cscaling < ainscale) {
            const double s = convert_scaling(zeta);
            for (int j = 0; j < n; j++)
                for (int i = 0; i < m; i++)
                    C[i + j * ldC] = s * C[i + j * ldC];
        }
        else if (ainscale < cscaling) {
            const double s = compute_combined_upscaling(ainscale, cscaling, zeta);
            for (int j = 0; j < n; j++)
                for (int i = 0; i < m; i++)
                    C[i + j * ldC] = s * C[i + j * ldC];
        }
        else {
            // A and C are consistently scaled.
            const double s = convert_scaling(zeta);
            for (int j = 0; j < n; j++)
                for (int i = 0; i < m; i++)
                    C[i + j * ldC] = s * C[i + j * ldC];
        }
    }

    // Update global scaling of Y.
#ifdef INTSCALING
    *cscale = min(cscaling, ainscale) + zeta;
#else
    *cscale = minf(cscaling, ainscale) * zeta;
#endif


    ////////////////////////////////////////////////////////////////////////////
    // Compute update.
    ////////////////////////////////////////////////////////////////////////////

    //  C := C - alpha * A  * B.
    //(mxn)            (mxk)(kxn)
    dgemm('N', 'N',
          m, n, k,
          -alpha, A, ldAin,
          B, ldB,
          1.0, C, ldC);


    omp_unset_lock(lock);

    ////////////////////////////////////////////////////////////////////////////
    // Free workspace.
    ////////////////////////////////////////////////////////////////////////////

    if (rescale_A) {
        _mm_free(A);
    }
}
