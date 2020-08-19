#include "norm.h"
#include "defines.h"
#include "utils.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>


// Credits: Björn Adlerborn
double vector_infnorm(int n, const double *x)
{
    double norm = 0;
    for (int i = 0; i < n; ++i) {
        double abs =  fabs(x[i]);
        if (abs > norm) {
            norm = abs;
        }
    }
    return norm;
}


double vector_cmplx_infnorm(int n, const double *x_re, const double *x_im)
{
    double norm = 0;
    for (int i = 0; i < n; ++i) {
        // Compute len = sqrt(x_re[i] * x_re[i] + x_im[i] * x_im[i]) robustly.
        double maxabs = maxf(fabs(x_re[i]), fabs(x_im[i]));
        double len = maxabs * sqrt(  (x_re[i] / maxabs) * (x_re[i] / maxabs)
                                   + (x_im[i] / maxabs) * (x_im[i] / maxabs));

        if (len > norm) {
            norm = len;
        }
    }
    return norm;
}


double vector_2norm(int n, const double *x)
{
    double norm = 0;
    for (int i = 0; i < n; ++i) {
        norm += x[i] * x[i];
    }
    return sqrt(norm);
}


// Credits: Björn Adlerborn, slightly modified.
double matrix_infnorm(int n, int m, const double *A, int ldA)
{
#define A(i,j) A[(i) + (j) * ldA]

    double *rowsums = calloc(n, sizeof(double));

    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; i++) {
            rowsums[i] += fabs(A(i,j)); 
        }
    }

    double norm = rowsums[0];
    for (int i = 1; i < n; i++) {
        if (rowsums[i] > norm) {
            norm = rowsums[i];
        }
    }

    free(rowsums);
    return norm;

#undef A
}






double matrix_tile_layout_infnorm(double ***A_tiles,
    int num_row_blks, int num_col_blks,
    const int *restrict const first_row,
    const int *restrict const first_col)
{
    // Find the row count of the matrix.
    const int n = first_row[num_row_blks];

    double *rowsums = calloc(n, sizeof(double));

    for (int j = 0; j < num_col_blks; j++) {
        for (int i = 0; i < num_row_blks; i++) {
            // Find block properties.
            double *A = A_tiles[i][j];
            const int num_rows = first_row[i + 1] - first_row[i];
            const int num_cols = first_col[j + 1] - first_col[j];

            // Loop over entries in block.
            for (int jj = 0; jj < num_cols; jj++) {
                for (int ii = 0; ii < num_rows; ii++) {
                    const int row = first_row[i] + ii;
                    rowsums[row] += fabs(A[ii + num_rows * jj]);
                }
            }
        }
    }

    double norm = rowsums[0];
    for (int i = 1; i < n; i++) {
        if (rowsums[i] > norm) {
            norm = rowsums[i];
        }
    }

    free(rowsums);
    return norm;

}


double tiled_matrix_frobeniusnorm(
    double ***X_tiles, int ldX,
    partitioning_t *p,
    memory_layout_t layout)
{
    double norm = 0.0;

    int *first_row = p->first_row;
    int *first_col = p->first_col;
    int num_blk_cols = p->num_blk_cols;
    int num_blk_rows = p->num_blk_rows;

    if (layout == COLUMN_MAJOR) {
        const int m = first_row[num_blk_rows];
        const int n = first_col[num_blk_cols];
        double norm = dlange('F', m, n, X_tiles[0][0], ldX);

        return norm;
    }
    else { //     TILE_LAYOUT
        // Loop over tiles.
        for (int blkj = 0; blkj < num_blk_cols; blkj++) {
            for (int blki = 0; blki < num_blk_rows; blki++) {
                // Locate base pointer of current block.
                double *X = X_tiles[blki][blkj];

                // Compute actual number of rows and columns.
                const int m = first_row[blki + 1] - first_row[blki];
                const int n = first_col[blkj + 1] - first_col[blkj];

                // Add squared Frobenius norm of Xij.
                double normF = dlange('F', m, n, X, m);
                norm += normF * normF;
            }
        }

        return sqrt(norm);
    }
}

