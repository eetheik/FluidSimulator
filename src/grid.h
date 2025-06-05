// /src/grid.h
#ifndef GRID_H
#define GRID_H

typedef struct {
    int N; // Size of grid N x N
    float *u; // x-velocity
    float *v; // y-velocity
    float *f_x; // x-force
    float *f_y; // y-force
    float *D; // dye density
    float *D0; // temporary array for dye density
} Grid;

Grid* grid_create(int N);

void grid_clear(Grid *grid);

void grid_free(Grid *grid);

#endif