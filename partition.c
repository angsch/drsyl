#include "partition.h"
#include "defines.h"
#include "utils.h"

#include <assert.h>
#include <stdio.h>



void partition(
    int n, const int *restrict const lambda_type,
    int num_tiles, int tlsz, int *restrict const first_row)
{
    // Compute candidate partitioning.
    for (int i = 0; i < num_tiles; i++)
        first_row[i] = i * tlsz;

    // Fill pad so that #rows = first_row[i + 1] - first_row[i].
    first_row[num_tiles] = n;

    // Count complex eigenvalues per tile.
    int num_cmplx = 0;

    // Absolute indices.
    int first_idx = 0;
    int last_idx = min(first_idx + tlsz, n);

    for (int k = 0; k < num_tiles; k++) {
        // Count complex eigenvalues in the candidate tile.
        for (int i = first_idx; i < last_idx; i++)
            if (lambda_type[i] == CMPLX)
                num_cmplx++;

        if ((num_cmplx % 2) == 0) {
            // Accept candidate partitioning.
            first_idx = last_idx;
            last_idx = min(first_idx + tlsz, n);
            num_cmplx = 0;
        }
        else {
            // (num_cmplx % 2) == 1

            // 2-by-2 block is split across tiles. Adapt candidate partitioning.
            {
                // If this tile is the last one, skip the remaining computation.
                if (k == num_tiles - 1) {
                    continue;
                }

                // Let the next tile start one entry later.
                first_row[k + 1]++;
            }

            // Advance to the next tile.
            first_idx = last_idx + 1;
            last_idx = min(first_idx + tlsz - 1, n);
            num_cmplx = 0;
        }
    }

    // Crude fix if n % tlsz == 1 AND last block is 2-by-2.
    if (first_row[num_tiles - 1] > first_row[num_tiles])
        first_row[num_tiles - 1] = first_row[num_tiles];

}


void partition_matrix(
    double *restrict const A, int ldA,
    memory_layout_t layout,
    const partitioning_t *restrict const p,
    double *** A_tiles)
{
    // Extract row and column partitioning.
    const int *first_row = p->first_row;
    const int *first_col = p->first_col;
    const int num_blk_rows = p->num_blk_rows;
    const int num_blk_cols = p->num_blk_cols;

    switch (layout) {
    case COLUMN_MAJOR:
    {
        #define A(i,j) A[(i) + (j) * ldA]
        for (int i = 0; i < num_blk_rows; i++) {
            for (int j = 0; j < num_blk_cols; j++) {
                A_tiles[i][j] = &A(first_row[i], first_col[j]);
            }
        }
        #undef A
    }
    break;

    case TILE_LAYOUT: {
        // Use column major order to store blocks.
        for (int i = 0; i < num_blk_rows; i++) {
            for (int j = 0; j < num_blk_cols; j++) {
                A_tiles[i][j]
                    = A + first_row[num_blk_rows] * first_col[j] // Full block cols to the left of us.
                        + first_row[i] * (first_col[j+1] - first_col[j]); // Offset in our block column.
            }
        }
    }
    break;

    default:
    {
        assert(0);
    }
    break;
    }
}
