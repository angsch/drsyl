#include "robust.h"
#include "utils.h"
#include "defines.h"
#include "norm.h"

#include <math.h>
#include <stdio.h>
#include <float.h>
#include <math.h>


int MIN_EXP = DBL_MIN_EXP - 1; // -1022
int MAX_EXP = DBL_MAX_EXP - 1; //  1023

// Overflow threshold.
const double g_omega = 1.e+307;
const double g_omega_inv = 1.e-307;


////////////////////////////////////////////////////////////////////////////////
// protect real division
////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Compute scaling such that the division b / t cannot overflow
 * where b, t are real-valued.
 *
 * If the return type is double-prevision, this routine returns a scaling alpha
 * such that x = (alpha * b) / t cannot overflow.
 *
 * If the return type is int, this routine returns a scaling alpha such that
 * x = (2^alpha * b) / t cannot overflow.
 *
 * Assume |b|, |t| are bounded by Omega.
 *
 * Credits: Carl Christian Kjelgaard Mikkelsen.
 */
static double protect_real_division(double b, double t)
{
    // Initialize scaling factor.
    double scale = 1.0;

    // Find scaling alpha such that x = (alpha * b) / t cannot overflow.
    if (fabs(t) < g_omega_inv) {
        if (fabs(b) > fabs(t) * g_omega) {
            // Please observe that scales will be strictly less than 1.
            scale = (fabs(t) * g_omega) / fabs(b);
        }
    }
    else { // fabs(t) >= g_omega_inv
        // Exploit short circuiting, i.e., the left side is evaluated first.
        // If 1.0 > abs(t) holds, then it is safe to compute
        // fabs(t) * g_omega.
        if (1.0 > fabs(t) && fabs(b) > fabs(t) * g_omega) {
            scale = 1.0 / fabs(b);
        }
    }

    return scale;
}

////////////////////////////////////////////////////////////////////////////////
// protect sum
////////////////////////////////////////////////////////////////////////////////

// Returns scaling such that sum := (alpha * x) + (alpha * y) cannot overflow.
double protect_sum(double x, double y)
{
    double scale = 1.0;

    // Protect against overflow if x and y have the same sign.
    if ((x > 0 && y > 0) || (x < 0 && y < 0))
        if (fabs(x) > g_omega - fabs(y))
            scale = 0.5;

    return scale;
}



////////////////////////////////////////////////////////////////////////////////
// protect multiplication (internal)
////////////////////////////////////////////////////////////////////////////////

// Returns scaling alpha such that y := t * (alpha * x) cannot overflow.
static double protect_mul(double tnorm, double xnorm)
{
    // Initialize scaling factor.
    double scale = 1.0;

    // Process simplified decision tree of protect_update().
    if (fabs(xnorm) <= 1.0) {
        if (fabs(tnorm) * fabs(xnorm) > g_omega) {
            scale = 0.5;
        }
    }
    else { // xnorm > 1.0
        if (fabs(tnorm) > g_omega / fabs(xnorm)) {
            scale = 0.5 / fabs(xnorm);
        }
    }

    return scale;
}


////////////////////////////////////////////////////////////////////////////////
// protect update
////////////////////////////////////////////////////////////////////////////////

// Returns scaling alpha such that y := (alpha * y) - t * (alpha * x) cannot
// overflow.
scaling_t protect_update(double tnorm, double xnorm, double ynorm)
{
    // Initialize scaling factor.
    double scale = 1.0;

    // Process decision tree.
    if (xnorm <= 1.0) {
        if (tnorm * xnorm > g_omega - ynorm) {
            scale = 0.5;
        }
    }
    else { // xnorm > 1.0
        if (tnorm > (g_omega - ynorm) / xnorm) {
            scale = 0.5 / xnorm;
        }
    }

#ifdef INTSCALING
    return ilogb(scale);
#else
    return scale;
#endif
}



////////////////////////////////////////////////////////////////////////////////
// protect update scalar
////////////////////////////////////////////////////////////////////////////////

static double protect_update_scalar(double t, double x, double y)
{
    double scale = 1.0;

    // Protect p = x * y.
    double alpha1 = protect_mul(x, t);
    double p = t * (alpha1 * x);
    if (fabs(ilogb(y) - ilogb(p)) > 52) {
        // The factors are far apart. Either y or p is the final result.
        if (ilogb(p) > ilogb(y))
            scale = alpha1;
    }
    else {
        // Scale y consistently.
        y = y / alpha1;
        double alpha2 = protect_sum(y, -p);
        scale = alpha1 * alpha2;
    }

    return scale;
}



////////////////////////////////////////////////////////////////////////////////
// [ a ] * [ x ] + sgn * [ x ] * [ b ] = scale * [ c ]
////////////////////////////////////////////////////////////////////////////////


void solve_a1x1_b1x1(const double a, const double sgn, const double b,
    const double smin, double *c, scaling_t /* == int*/ *scale)
{
    // Compute a + sgn * b robustly. Note that the scaling contributes as
    // reciprocal to the global scaling.
    double s = protect_sum(a, sgn * b);
    double t = (s * a) + sgn * (s * b);

    // Replace entries with too small magnitude. The problem is ill-conditioned.
    if (fabs(t) < smin) {
        printf("WARNING: The eigenvalues are very close.\n");
        t = copysign(smin, t);
    }

    // Compute a scaling to survive the real-valued division.
    double alpha = protect_real_division(c[0], t);

    // Execute the division safely.
    c[0] = (alpha * c[0]) / t;

    // Return scaling factor.
#ifdef INTSCALING
    scale[0] = ilogb(alpha / s);
#else
    scale[0] = alpha / s;
#endif
}



////////////////////////////////////////////////////////////////////////////////
// [ a11 a12 ] * [ x11 ] + sgn * [ x11 ] * b = [ c11 ]
// [ a21 a22 ]   [ x21 ]         [ x21 ]       [ c21 ]
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// solve 2x2 real system
////////////////////////////////////////////////////////////////////////////////


// Credits: Carl Christian Kjelgaard Mikkelsen
static double backsolve_real_2x2_system(
    double *T, int ldT, double *b, const double smin)
{
#define T(i,j) T[(i) + (j) * ldT]

    // Global scaling factor.
    double alpha = 1.0;

    double xnorm = max(fabs(b[0]), fabs(b[1]));

    if (fabs(T(1,1)) < smin) {
        // Perturb to have a sufficiently large divisor.
        T(1,1) = smin;
    }

    double s = protect_real_division(b[1], T(1,1));
    if (s != 1.0) {
        // Apply scaling to right-hand side.
        b[0] = s * b[0];
        b[1] = s * b[1];

        // Update global scaling.
        alpha = s * alpha;

        // Update the infinity norm of the solution.
        xnorm = s * xnorm;
    }

    // Execute the division.
    b[1] = b[1] / T(1,1);

#ifdef INTSCALING
    s = ldexp(1.0, protect_update(fabs(T(0,1)), fabs(b[1]), xnorm));
#else
    s = protect_update(fabs(T(0,1)), fabs(b[1]), xnorm);
#endif

    if (s != 1.0) {
        // Apply scaling to right-hand side.
        b[0] = s * b[0];
        b[1] = s * b[1];

        // Update global scaling.
        alpha = s * alpha;
    }

    // Execute the linear update.
    b[0] = b[0] - b[1] * T(0,1);

    // Recompute norm.
    xnorm = max(fabs(b[0]), fabs(b[1]));

    if (fabs(T(0,0)) < smin) {
        // Perturb to have a sufficiently large divisor.
        T(0,0) = smin;
    }

    s = protect_real_division(b[0], T(0,0));
    if (s != 1.0) {
        // Apply scaling to right-hand side.
        b[0] = s * b[0];
        b[1] = s * b[1];

        // Update global scaling.
        alpha = s * alpha;

        // Update the infinity norm of the solution.
        xnorm = s * xnorm;
    }

    // Execute the division.
    b[0] = b[0] / T(0,0);

    return alpha;

#undef T
}

// Swap row 0 and row 1.
static void swap_rows(int n, double *C)
{
#define C(i,j) C[(i) + (j) * 2]

    // Swap row 0 and row 1.
    for (int j = 0; j < n; j++) {
        double swap = C(0,j);
        C(0,j) = C(1,j);
        C(1,j) = swap;
    }

#undef C
}


static void find_real_pivot(double *C, int *pivot_row, int *pivot_col)
{
#define C(i,j) C[(i) + (j) * 2]

    // Find the coordinates of the pivot element.
    int row = 0;
    int col = 0;
    double cmax = 0.0;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            double lmax = fabs(C(i,j));
            if (lmax > cmax) {
                row = i;
                col = j;
                cmax = lmax;
            }
        }
    }

    *pivot_row = row;
    *pivot_col = col;

#undef C
}


// Complete pivoting.
static void solve_2x2_real_system_internal(
    const double *restrict const T, int ldT, 
    double lambda,
    double *restrict const b, double *restrict const scale,
    const double smin)
{
#define T(i,j) T[(i) + (j) * ldT]
#define C(i,j) C[(i) + (j) * 2]

    // Solve
    // (T - lambda I)x = b.

    // C = [(T - lambda * I) | b]
    double C[2 * 3];

    // Compute t + (-lambda) robustly. Recall the diagonals in the the 2-by-2
    // T block are equal, so that we have to protect only one subtraction.
    double s = protect_sum(T(0,0), -lambda);
    double csr = (s * T(0,0)) - (s * lambda);

    // Apply scaling to T. Note that scaling of b is not safe. Therefore s is
    // incorporated into the global scaling at the very end of this routine.
    // C := [s * (T - lambda I) | b].
    C(0,0) = csr;         C(0,1) = s * T(0,1);  C(0,2) = b[0];
    C(1,0) = s * T(1,0);  C(1,1) = csr;         C(1,2) = b[1];

    ////////////////////////////////////////////////////////////////////////////
    // Transform A to echelon form with complete pivoting.
    ////////////////////////////////////////////////////////////////////////////

    // Find pivot element in entire matrix.
    int pivot_row = 0, pivot_col = 0;
    find_real_pivot(C, &pivot_row, &pivot_col);

    // Permute pivot to the top-left corner.
    if (pivot_row == 1) {
        // Swap row 0 and row 1.
        swap_rows(3, C);
    }
    if (pivot_col == 1) {
        // Swap column 0 and column 1.
        for (int i = 0; i < 2; i++) {
            double swap = C(i,0);
            C(i,0) = C(i,1);
            C(i,1) = swap;
        }
    }

    if (fabs(C(0,0)) < smin) {
        // Perturb the pivot element to be sufficiently large.
        C(0,0) = smin;
    }

    // Compute multiplier, the reciprocal of the pivot.
    double ur11r = 1.0 / C(0,0);

    // Multiply first row with reciprocal of C(0,0).
    {
    C(0,0) = 1.0;
    C(0,1) = C(0,1) * ur11r;

    // Treat rhs.
    double beta = protect_mul(C(0,2), ur11r);
    *scale = beta;
    C(0,2) = C(0,2) * beta;
    C(1,2) = C(1,2) * beta;
    C(0,2) = C(0,2) * ur11r;
    }

    // Second row - CR(1,0) * first_row.
    {
    C(1,1) = C(1,1) - C(1,0) * C(0,1);

    // Treat rhs.
    double beta = protect_update_scalar(C(1,0), C(0,2), C(1,2));
    *scale = (*scale) * beta;
    C(0,2) = C(0,2) * beta;
    C(1,2) = C(1,2) * beta;
    C(1,2) = C(1,2) - C(1,0) * C(0,2);

    // (1,0) has been annihilated.
    C(1,0) = 0.0;
    }

    // The system is now in upper triangular form.

    ////////////////////////////////////////////////////////////////////////////
    // Backward substitution.
    ////////////////////////////////////////////////////////////////////////////

    double alpha = backsolve_real_2x2_system(&C(0,0), 2, &C(0,2),smin);
    *scale = (*scale) * alpha;

    // Copy the solution back.
    if (pivot_col == 1) {
        b[0] = C(1,2);
        b[1] = C(0,2);
    }
    else {
        b[0] = C(0,2);
        b[1] = C(1,2);
    }

#undef T
#undef C
}



void solve_a2x2_b1x1(
    const double *restrict const A, int ldA,
    const double sgn,
    const double b,
    const double smin,
    double *restrict const C,
    scaling_t *restrict const scale)
{
    // Solve system
    // ([ a11 a12 ] - (-sgn * b) I ) = [ c1 ]
    //  [ a21 a22 ]                    [ c2 ]

#ifdef INTSCALING
    // Local scaling factor.
    double phi = 1.0;

    const double lambda = -sgn * b;
    solve_2x2_real_system_internal(A, ldA, lambda, C, &phi, smin);

    // Convert double-precision scaling factor to int scaling factor.
    *scale = ilogb(phi);
#else
    const double lambda = -sgn * b;
    solve_2x2_real_system_internal(A, ldA, lambda, C, scale, smin);
#endif
}



////////////////////////////////////////////////////////////////////////////////
// a * [ x1 x2 ] + sgn * [ x1 x2 ] * [ b11 b12 ] = scale * [ c1 ]
//                                   [ b21 b22 ]           [ c2 ]
////////////////////////////////////////////////////////////////////////////////

void solve_a1x1_b2x2(
    const double a,
    const double sgn,
    const double *restrict const B, int ldB,
    const double smin,
    double *restrict const c, int ldc,
    scaling_t *restrict const scale)
{
#define B(i,j) B[(i) + (j) * ldB]

    // Construct system
    // [ sgn * b11 - a     sgn * b21     ] * [ x1 ] = [ c1 ]
    // [ sgn * b12         sgn * b22 - a ]   [ x2 ]   [ c2 ].
    // aka
    // (D - a * I) x = y.
    double D[2*2]; double y[2];
    D[0] = sgn * B(0,0);
    D[1] = sgn * B(0,1);
    D[2] = sgn * B(1,0);
    D[3] = sgn * B(1,1);
    y[0] = c[0];
    y[1] = c[ldc];

    // Internal scaling factor.
    double s;

    solve_2x2_real_system_internal(D, 2, -a, y, &s, smin);

#ifdef INTSCALING
    scale[0] = ilogb(s);
#else
    scale[0] = s;
#endif

    // Transpose solution vector.
    c[0]   = y[0];
    c[ldc] = y[1];

#undef B
}


////////////////////////////////////////////////////////////////////////////////
// [ a11 a12 ] * [ x11 x12 ] + sgn * [ x11 x12 ] * [ b11 b12 ] = [ c11 c12 ]
// [ a21 a22 ]   [ x21 x22 ]         [ x21 x22 ]   [ b21 b22 ] = [ c21 c22 ]
////////////////////////////////////////////////////////////////////////////////

static void find_pivot(int n, double *C, int ld, int *pivot_row, int *pivot_col)
{
#define C(i,j) C[(i) + (j) * ld]

    // Find the coordinates of the pivot element.
    int row = 0;
    int col = 0;
    double cmax = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double lmax = fabs(C(i,j));
            if (lmax > cmax) {
                row = i;
                col = j;
                cmax = lmax;
            }
        }
    }

    *pivot_row = row;
    *pivot_col = col;

#undef C
}


static void swap_rows_4(int row1, int row2, int n, double *C, int ld)
{
#define C(i,j) C[(i) + (j) * ld]

    if (row1 != row2) {
        // Swap C(row1,:) and C(row2,:).
        for (int j = 0; j < n; j++) {
            double tmp = C(row1,j);
            C(row1,j) = C(row2,j);
            C(row2,j) = tmp;
        }
    }

#undef C
}


static void swap_cols_4(int col1, int col2, int n, double *C, int ld)
{
#define C(i,j) C[(i) + (j) * ld]

    if (col1 != col2) {
        // Swap C(:,col1) and C(:,col2).
        for (int i = 0; i < n; i++) {
            double tmp = C(i,col1);
            C(i,col1) = C(i,col2);
            C(i,col2) = tmp;
        }
    }

#undef C
}




// Credits: Carl Christian Kjelgaard Mikkelsen
static double backsolve_real_4x4_system(double *T, int ldT, double *b)
{
#define T(i,j) T[(i) + (j) * ldT]

    // Global scaling factor.
    double alpha = 1.0;

    // Compute upper bound of T.
    const double tnorm = matrix_infnorm(4, 4, T, ldT);

    // Upper bound of the solution.
    double xnorm = vector_infnorm(4, b);

    // n = 4.
    for (int i = 3; i >= 0; i--) {
        double s = protect_real_division(b[i], T(i,i));
        if (s != 1.0) {
            // Apply scaling to right-hand side.
            b[0] = s * b[0];
            b[1] = s * b[1];
            b[2] = s * b[2];
            b[3] = s * b[3];

            // Update global scaling.
            alpha = s * alpha;

            // Update the infinity norm of the solution.
            xnorm = s * xnorm;
        }

        // Execute the division.
        b[i] = b[i] / T(i,i);

        // Are there more columns to process?
        if (i > 0) {
            s = convert_scaling(protect_update(tnorm, fabs(b[i]), xnorm));
            if (s != 1.0) {
                // Apply scaling to right-hand side.
                b[0] = s * b[0];
                b[1] = s * b[1];
                b[2] = s * b[2];
                b[3] = s * b[3];

                // Update global scaling.
                alpha = s * alpha;
            }

            // Execute the linear update.
            for (int j = 0; j < i; j++) {
                b[j] = b[j] - T(j,i) * b[i];
            }

            // Recompute norm.
            xnorm = vector_infnorm(i, b);
        }
    }

    return alpha;

#undef T
}





// A * X + sgn * X * B = scale * C.
void solve_a2x2_b2x2(
    const double *restrict const A, int ldA,
    double sgn,
    const double *restrict const B, int ldB,
    double *restrict const C, int ldC,
    scaling_t *restrict const scale)
{
#define T(i,j) T[(i) + (j) * 4]
#define A(i,j) A[(i) + (j) * ldA]
#define B(i,j) B[(i) + (j) * ldB]
#define C(i,j) C[(i) + (j) * ldC]

    // LAPACK DLASY2 does not magic to compute smin in this routine.
    // Machine epsilon according to Demmel and as used in LAPACK. Note that
    // this is half of the machine epsilon defined in the ISO C standard.
    const double eps = DBL_EPSILON / 2;
    double smin = maxf(find_absmax_in_2x2(A, ldA), find_absmax_in_2x2(B, ldB));
    smin = maxf(eps * smin, DBL_MIN);


    double T[4 * 4 + 4];

    // Construct 4-by-4 system matrix.
    T(0,0)=A(0,0)+sgn*B(0,0); T(0,1)=A(0,1);            T(0,2)=sgn*B(1,0);        T(0,3)=0.0;
    T(1,0)=A(1,0);            T(1,1)=A(1,1)+sgn*B(0,0); T(1,2)=0.0;               T(1,3)=sgn*B(1,0);
    T(2,0)=sgn*B(0,1);        T(2,1)=0.0;               T(2,2)=A(0,0)+sgn*B(1,1); T(2,3)=A(0,1);
    T(3,0)=0.0;               T(3,1)=sgn*B(0,1);        T(3,2)=A(1,0);            T(3,3)=A(1,1)+sgn*B(1,1);

    // Construct right-hand side.
    T(0,4) = C(0,0);
    T(1,4) = C(1,0);
    T(2,4) = C(0,1);
    T(3,4) = C(1,1);

    // Initialize local scaling factor.
    double s = 1.0;

    ////////////////////////////////////////////////////////////////////////////
    // Transform T to echelon form with complete pivoting.
    ////////////////////////////////////////////////////////////////////////////

    int pivot_row = 0, pivot_col = 0;
    int col_pivots[4];
    for (int i = 0; i < 4; i++) {
        // Find pivot element in T(i:n,i:n).
        find_pivot(4 - i, &T(i,i), 4, &pivot_row, &pivot_col);
        // Correct offset of pivot element.
        pivot_row += i; pivot_col += i;
        col_pivots[i] = pivot_col;

        // Permute pivot to position T(i,i). The permutation encloses the rhs.
        swap_rows_4(pivot_row, i, 5, T, 4);
        swap_cols_4(pivot_col, i, 4, T, 4);

        if (fabs(T(i,i)) < smin) {
            // Perturb the pivot element to be sufficiently large.
            T(i,i) = smin;
        }

        // Compute multiplier, the reciprocal of the pivot.
        double tiir = 1.0 / T(i,i);

        // Multiply i-th row with reciprocal of T(i,i).
        {
            // Treat matrix.
            T(i,i) = 1.0;
            for (int j = i + 1; j < 4; j++)
                T(i,j) = T(i,j) * tiir;

            // Treat rhs.
            double beta = protect_mul(T(i,4), tiir);
            s = s * beta;
            for (int ii = 0; ii < 4; ii++)
                T(ii,4) = T(ii,4) * beta;
            T(i,4) = T(i,4) * tiir;
        }

        // Annihilate T(k,i) for k = i + 1, ..., 4 via T(k,:) - T(k,i) * T(i,:).
        {
            for (int k = i + 1; k < 4; k++) {
                for (int jj = i + 1; jj < 4; jj++)
                    T(k,jj) = T(k,jj) - T(k,i) * T(i,jj);

                // Treat rhs.
                double beta = protect_update_scalar(T(k,i), T(i,4), T(k,4));
                s = s * beta;
                for (int ii = 0; ii < 4; ii++)
                    T(ii,4) = T(ii,4) * beta;
                T(k,4) = T(k,4) - T(k,i) * T(i,4);

                // (k,i) has been annihilated.
                T(k,i) = 0.0;
            }
        }
    }

    // The system is now in upper triangular form.

    ////////////////////////////////////////////////////////////////////////////
    // Backward substitution.
    ////////////////////////////////////////////////////////////////////////////

    double alpha = backsolve_real_4x4_system(&T(0,0), 4, &T(0,4));
    s = s * alpha;

    // Apply column pivoting.
    for (int i = 0; i < 4; i++) {
        if (col_pivots[3 - i] != 3 - i) {
            double tmp = T(3-i,4);
            T(3-i,4) = T(col_pivots[3 - i],4);
            T(col_pivots[3 - i],4) = tmp;
        }
    }

    // Copy the solution back.
    C(0,0) = T(0,4);
    C(1,0) = T(1,4);
    C(0,1) = T(2,4);
    C(1,1) = T(3,4);



#ifdef INTSCALING
    *scale = ilogb(s);
#else
    *scale = s;
#endif

#undef T
#undef A
#undef B
#undef C
}

