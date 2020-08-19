#ifndef ROBUST_H
#define ROBUST_H

#include "typedefs.h"


/**
 * @brief Computes a scaling such that x + y does not overflow.
 *
 * The routine computes a scaling such that (scale * x) + (scale * y)
 * does not overflow.
 *
 * @param[in]  x    Scalar x.
 * @param[in]  y    Scalar y.
 * @return A scaling factor.
 */
double protect_sum(double x, double y);



/**
 * @brief Computes scaling such that the update y := y - t x cannot overflow.
 *
 * If the return type is of type double, this routine
 * returns a scaling alpha such that y := (alpha * y) - t * (alpha * x)
 * cannot overflow.
 *
 * If the return type is of type int, this routine
 * returns a scaling alpha such that y := (2^alpha * y) - t * (2^alpha * x)
 * cannot overflow.
 *
 * Assume 0 <= t, x, y <= Omega.
 *
 * Credits: Carl Christian Kjelgaard Mikkelsen.
 */
scaling_t protect_update(double tnorm, double xnorm, double ynorm);



/**
 * @brief Solves a * x + sgn * x * b = c robustly.
 *
 * If the type of scale is double, the routine solves (scale * c) / (a + sgn * b)
 * whereas, if the type of scale is int, the routine solves
 * (2^scale * c) / (a + sgn * b) such that no overflow occurs.
 *
 * @param[in]      a       Real scalar a.
 * @param[in]      sgn     +1.0 or -1.0
 * @param[in]      b       Real scalar b.
 * @param[in]      smin    Safe minimum.
 * @param[in, out] c       On entry, the right-hand side c. On exit, the real
 *                         solution (scale * x) respectively (2^scale * x).
 * @param[out]     scale   Scalar scaling factor of x.
 */
void solve_a1x1_b1x1(const double a, const double sgn, const double b,
    const double smin, double *c, scaling_t *scale);



/**
 * @brief Solve A * x + sgn * x * b = c where A is 2-by-2 and b is 1-by-1.
 *
 * Solves the real-valued system
 * [ a11 a12 ] * [ x1 ] + sgn * [ x1 ] * b = [ c1 ]
 * [ a21 a22 ]   [ x2 ]         [ x2 ]       [ c2 ]
 * such that it cannot overflow.
 * Internally, the system is converted into the eigenvalue problem
 * [ a11-lambda  a12        ] * [ x1 ] = [ c1 ]
 * [ a21         a22-lambda ]   [ x2 ]   [ c2 ]
 * where lambda := - sgn * b. The system is converted into echelon form through
 * complete pivoting. The resulting triangular system is solved through backward
 * substitution. All operations are robust.
 *
 * @param[in]      A       Real 2-by-2 matrix A.
 * @param[in]      ldA     The leading dimension of A. ldA >= 2.
 * @param[in]      sgn     +1.0 or -1.0.
 * @param[in]      b       The scalar b.
 * @param[in]      smin    Safe minimum.
 * @param[in, out] c       Real vector of length 2. On entry, the right-hand
 *                         side. On exit, the solution.
 * @param[out]     scale   Scalar scaling factor of the solution x.
 * 
 * Credits: Carl Christian Kjelgaard Mikkelsen.
 */
void solve_a2x2_b1x1(
    const double *restrict const A, int ldA,
    const double sgn,
    const double b,
    const double smin,
    double *restrict const c, scaling_t *restrict const scale);


/**
 * @brief Solve a * x + sgn * x * B = c where a is 1-by-1 and B is 2-by-2.
 *
 * Solves the real-valued system
 * a * [ x1 x2 ] + sgn * [ x1 x2 ] * [ b11 b12 ] = [ c1 ]
 *                                   [ b21 b22 ]   [ c2 ]
 * such that if cannot overflow.
 * Internally, the transposed system is interpreted as the eigenvalue problem
 * [ sgn * b11 - a     sgn * b21     ] * [ x1 ] = [ c1 ]
 * [ sgn * b12         sgn * b22 - a ]   [ x2 ]   [ c2 ].
 * This system is transformed into echelon form through complete pivoting. The
 * resulting triangular system is solved through backward substitution. All
 * operations are robust.
 *
 * @param[in]      a       The scalar a.
 * @param[in]      sgn     +1.0 or -1.0.
 * @param[in]      B       Real 2-by-2 matrix B.
 * @param[in]      ldB     The leading dimension of B. ldB >= 2.
 * @param[in]      smin    Safe minimum.
 * @param[in, out] c       Row vector of length 2. On entry, the right-hand
 *                         side. On exit, the solution.
 * @param[in]      ldc     The leading dimension of c. ldc >= 1.
 * @param[out]     scale   Scalar scaling factor of the solution x.
 */
void solve_a1x1_b2x2(
    const double a,
    const double sgn,
    const double *restrict const B, int ldB,
    const double smin,
    double *restrict const c, int ldc,
    scaling_t *restrict const scale);




/**
 * @brief Solve A * X + sgn * X * B = C where all matrices are 2-by-2.
 *
 * Solves the real-valued system
 * [ a11 a12 ] * [ x11 x12 ] + sgn * [ x11 x12 ] * [ b11 b12 ] = [ c11 c12 ]
 * [ a21 a22 ]   [ x21 x22 ]         [ x21 x22 ]   [ b21 b22 ] = [ c21 c22 ]
 * such that overflow cannot occur.
 *
 * @param[in]      A       Real 2-by-2 matrix A.
 * @param[in]      ldA     The leading dimension of A. ldA >= 2.
 * @param[in]      sgn     The sign.
 * @param[in]      B       Real 2-by-2 matrix B.
 * @param[in]      ldB     The leading dimension of B. ldB >= 2.
 * @param[in, out] C       Real 2-by-2 matrix. On entry, C. On exit, the
 *                         solution X.
 * @param[in]      ldC     The leading dimension of C. ldC >= 2.
 * @param[out]     scale   Scalar scaling factor of the solution X.
 */
void solve_a2x2_b2x2(
    const double *restrict const A, int ldA,
    double sgn,
    const double *restrict const B, int ldB,
    double *restrict const C, int ldC,
    scaling_t *restrict const scale);

#endif
