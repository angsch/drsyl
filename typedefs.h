#ifndef TYPEDEF_H
#define TYPEDEF_H

enum memory_layout {COLUMN_MAJOR, TILE_LAYOUT};
typedef enum memory_layout memory_layout_t;

enum matrix_desc {UPPER_TRIANGULAR, LOWER_TRIANGULAR};
typedef enum matrix_desc matrix_desc_t;

typedef struct {
    int num_blk_rows; // num_blk_rows + 1 exists
    int num_blk_cols; // num_blk_cols + 1 exists
    int *first_row;
    int *first_col;
} partitioning_t;



#ifdef INTSCALING
typedef int scaling_t;
#else
typedef double scaling_t;
#endif

#endif
