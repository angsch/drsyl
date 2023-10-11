# Robust Task-Parallel Solution of the Triangular Sylvester Equation 

The robust solver `drsyl` solves the scaled triangular Sylvester equation $AX + XB = \alpha C$ in a tiled and task-parallel fashion. The matrices $A$ and $B$ are assumed to be upper quasi-triangular; the right-hand side matrix $C$ and the solution matrix $X$ are dense. The scalar $\alpha$ is computed alongside with the solution $X$ such that overflow is avoided and the solution matrix $X$ at no point during the computation contains infinities representing overflow. This property qualifies `drsyl` as robust. The name `drsyl` is inspired by the LAPACK naming convention and stands for *d*ouble precision *r*obust triangular *Syl*vester equation solver.

## Prerequisites

* A compiler that supports OpenMP. The compiler optimization level must be chosen such that associative math is disabled.
* An efficient BLAS implementation.
* An implementation of the LAPACK routines. Only needed for the validation.

## Building and executing

### Makefile
A Makefile for the GNU compiler linked against OpenBLAS and the Intel compiler linked against MKL is included. Try `make`. 

### CMake

Navigate into the directory and run the following.
```
mkdir build
cd build
cmake ..
make -j
```
### Execution

The executable `drsyl` takes 8 input parameters.

* `m`. The matrix size of A and the row count of C and X.
* `n`. The matrix size of B and the column count of C and X.
* `tlsz`. The tile size. A good starting value is 400. For a good performance, the tile size needs to be tuned.
* `cmplx-ratio-A`. The proportion of 2-by-2 blocks on the diagonal of A.
* `cmplx-ratio-B`. The proportion of 2-by-2 blocks on the diagonal of B.
* `sign`. The sign of the matrix B in the Sylvester equation.
* `mem-layout`. The storage format of the matrices, 0=column major or 1=tile layout.
* `seed`. The seed for the random number generator to reproduce runs.

For a quick test, type `./drsyl 5000 5000 400 0.5 0.5 1 0 0`. The solver uses task parallelism. The number of threads is controlled with `OMP_NUM_THREADS`. For best performance, it is recommended to pin the threads.

## Remarks

By default, a double-precision number is used for the scaling factor $\alpha$. Sometimes the systems grow so quickly that overflow protection with a double-precision scaling factor does not suffice. Then setting `-DINTSCALING` during the build process activates integer scaling factors and allows for solving systems that are not solvable by a double-precision scaling factor. This change requires a complete rebuild (`make clean`, `make`).

More information can be found [here](https://people.cs.umu.se/angies/sylvester).

## Publication Reference

Schwarz, A., Kjelgaard Mikkelsen, C.C. (2020). Robust Task-Parallel Solution of the Triangular Sylvester Equation. In: Wyrzykowski, R., Deelman, E., Dongarra, J., Karczewski, K. (eds) Parallel Processing and Applied Mathematics. PPAM 2019. Lecture Notes in Computer Science(), vol 12043. Springer, Cham.  https://doi.org/10.1007/978-3-030-43229-4_8
