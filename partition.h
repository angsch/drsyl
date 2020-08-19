#ifndef PARTITION_H
#define PARTITION_H

#include "typedefs.h"

/**
 * Generate a partitioning that does not split 2x2 blocks.
 *
 * @param[in]  n           Dimension of the system
 * @param[in]  lambda_type Array of length n. Describes if the eigenvalue at
 *                         (i,i) is a real or part of a complex eigenvalue.
 * @param[in]  num_tiles   Number of tiles.
 * @param[in]  tlsz        Tile size.
 * @param[out] first_row   Partitioning of an array with n entries into segments
 *                         of approximately tlsz size such that complex eigen-
 *                         values are not split. On exit, first_row[i] holds
 *                         complex eigenvalues. On exit, first_row[i] holds
 *                         the start row index of the i-th tile row.
 */
void partition(
    int n, const int *restrict const lambda_type,
    int num_tiles, int tlsz, int *restrict const first_row);



/**
 * Apply partitioning to a matrix.
 * 
 * @param[in]  A         The matrix to be partitioned.
 * @param[in]  ldA       The leading dimension of A.
 * @param[in]  layout    COLUMN_MAJOR or TILE_LAYOUT.
 * @param[in]  p         Partitioning of A that contains row and column indices
 *                       of all blocks.
 * @param[out] A_tiles   On entry, a num_blks-by-num_blks array of pointers. On
 *                       exit, A_tiles[i][j] holds the base pointer to the
 *                       tile (i,j) according to the partitioning.
 */
void partition_matrix(
    double *restrict const A, int ldA,
    memory_layout_t layout,
    const partitioning_t *restrict const p,
    double ***A_tiles);


#endif
