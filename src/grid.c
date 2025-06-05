// src/grid.C
#include "grid.h"
#include <fftw3.h> // For fftwf_malloc and fftw_free
#include <string.h> // For memset 
#include <stdlib.h>

/**
 * @brief A function for allocating and initializing a Grid struct.
 * 
 * The function allocates memory for a Grid struct and its associated arrays
 * using fftwf_malloc(). The arrays are padded with two extra rows, so that
 * they support FFTW's real-to-complex transformations required in the FFTs.
 * After allocation the grid's arrays are initialized to all zeroes using grid_clear().
 *
 * @param N The logical size of the square array to be created (N * N). 
 * 
 * @return A pointer to the newly allocated and initialized Grid struct.
 */
Grid* grid_create(int N) {
    Grid *grid = (Grid* )malloc(sizeof(Grid));
    grid->N = N;
    const int arraysize = N * (N + 2); // fftw uses two extra rows

    // Allocate arrays using fftwf_malloc
    grid->u = (float *)(fftwf_malloc(sizeof(float) * arraysize));
    grid->v = (float *)(fftwf_malloc(sizeof(float) * arraysize));

    grid->f_x = (float *)(fftwf_malloc(sizeof(float) * arraysize));
    grid->f_y = (float *)(fftwf_malloc(sizeof(float) * arraysize));

    // No FFTs are applied on these but in case the user for some reason
    // desires to do so, I still allocate them with fftwf_malloc
    grid->D = (float *)(fftwf_malloc(sizeof(float) * arraysize));
    grid->D0 = (float *)(fftwf_malloc(sizeof(float) * arraysize));

    grid_clear(grid); // Ensures created grid is empty

    return grid;
}


/**
 * @brief A function for clearing a grid.
 * 
 * The function takes in a pointer to a Grid struct and if the grid exists, it 
 * uses memset() to set the N * (N + 2) block of floating point numbers associated 
 * with each of the grid's arrays to all zeroes.
 *
 * @param grid A pointer to the grid for which the arrays are to be cleared.
 */
void grid_clear(Grid *grid) {
    int N = grid->N;
    if (grid) {
        const int arraysize = N * (N + 2); // fftw uses two extra rows

        // Initialize arrays to zero (fftwf_malloc does not initialize memory)
        memset(grid->u, 0, sizeof(float) * arraysize);
        memset(grid->v, 0, sizeof(float) * arraysize);
        memset(grid->f_x, 0, sizeof(float) * arraysize);
        memset(grid->f_y, 0, sizeof(float) * arraysize);
        memset(grid->D, 0, sizeof(float) * arraysize); // Dye field likely won't be FT'd but maybe? Also consistency
        memset(grid->D0, 0, sizeof(float) * arraysize);
    }
}


/**
 * @brief A function for freeing the memory associated with a given grid.
 *
 * The function takes in a pointer to a Grid struct and if that grid exists,
 * it uses fftwf_free() to free the memory blocks created for each of the
 * grid's arrays. In general you want to use this function at the very end
 * of your simulation loop. This will help prevent any memory issues.
 *
 * @param grid A pointer to the grid for which the memory is to be freed.
 */
void grid_free(Grid *grid) {
    // For freeing memory at the end
    if (grid) {
        fftwf_free(grid->u);
        fftwf_free(grid->v);
        fftwf_free(grid->f_x);
        fftwf_free(grid->f_y);
        fftwf_free(grid->D);
        fftwf_free(grid->D0);
        free(grid);
    }
}