#include "../src/grid.h"
#include "../src/fluid_solver.h"
#include "../src/io.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/**
 * @brief A function defining the main simulation loop for the kelvin_helmholtz_instability simulation
 *
 * This simulation initializes a dye field with value 1.0 for the upper half of the simulation grid and
 * a value of 0.0 for the lower half of the grid. A small disturbance is applied at the boundary for the 
 * first 10 timesteps, resulting in the formation of swirls during the following evolution of the fluid. 
 * This simulation demonstrates the Kelvin- Helmholtz instability.
 * 
 * Exact details for the disturbance applied can be consulted from section 3.1 of the report.
 * 
 * Output binary files are stored in the /outputs directory with the file basename "kelvin_helmholtz_instability", 
 * followed by the date and time of the simulation, e.g. "kelvin_helmholtz_instability_new_19_05_2025_18h02m26s.bin"
 *
 */
int main() {
    // I am not sure how well this works on other terminals besides my own
    int render_to_terminal = 1; // Set this to 1 and it will render to terminal (zoom out in terminal) 
    
    // Specify simulation parameters (task 4)
    const int N = 500;
    const int total_steps = 700;
    float dt = 0.01;
    float visc = 0.001; 
    float U0 = 5.0;
    float delta = 0.025;
    float A = 1.0;
    float sigma = 0.02;
    int k = 4;

    char* filename = NULL;
    if (render_to_terminal != 1) {
        char* filename = init_binary_output_file(N, N, total_steps, "kelvin_helmholtz_instability");
    }

    Grid* grid = grid_create(N); // From grid we have access to u, v, u0, v0, D, D0

    // Necessary for initializing solver (could abstract away to some VelocitySolver object, which needs init.?)
    // Also would benefit from renaming to something like init_velocity_solver(N). I am torn between keeping the
    // name as is or changing it. The name describes it very explicitly, but it might not be obvious to a user
    // that the simulation needs a FFT even, so why would they initialize it. In any case the initialization should
    // not be done inside the velocity solver function since then we would reinitialize ever loop iteration
    init_FFT(N); 

    // Create requested x, y arrays:
    float *x_i = (float *)(malloc(sizeof(float) * N));
    float *y_j = (float *)(malloc(sizeof(float) * N));
    for (int i = 0; i < N; i++) {
        x_i[i] = (i + 0.5) / N;
        y_j[i] = (i + 0.5) / N;
    }

    // if y_j >= 1/2 then D_ij = 1.0, 
    // if y_j < 1/2 then D_ij = 0.0
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (y_j[j] <= 0.5) {
                grid->D[i + N*j] = 1.0;
            } else {
                grid->D[i + N*j] = 0.0;
            }
        }
    }

    // Simulation loop
    for (int step = 0; step < total_steps; step++) {
        // For first ten steps apply given force
        if (step < 10) {
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    grid->f_x[i + N * j] = U0 * tanh((y_j[j] - 0.5) / delta);
                    grid->f_y[i + N * j] = A * sin(2 * M_PI * k * x_i[i]) * exp(-(y_j[j] - 0.5) * (y_j[j] - 0.5) / (2*sigma*sigma));
                }
            }
        } else {
            // No force
            memset(grid->f_x, 0, sizeof(float) * N * (N + 2));
            memset(grid->f_y, 0, sizeof(float) * N * (N + 2));
        }

        // Solve velocity field
        evolve_velocity_field(grid->N, grid->u, grid->v, grid->f_x, grid->f_y, visc, dt);
        // Then dye field from velocity field
        evolve_dye_field(grid->N, grid->u, grid->v, grid->D, grid->D0, dt);
        
        if (render_to_terminal == 1) {
            render_colored_field(grid->N, grid->D, 4);
            usleep(11000);
        } else { 
            write_array_to_binary(grid->D, step, N, N, filename);
            if (step + 1 == total_steps) {
                printf("Filename stored to %s\n", filename);
            }
        }
        
    }

    grid_free(grid);
    deinit_FFT();
    return 0;
}