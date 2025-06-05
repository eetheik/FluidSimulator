#include "../src/grid.h"
#include "../src/fluid_solver.h"
#include "../src/io.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> 

// In this simulation we add four Gaussians to the density field, 
// and input a force that sends them towards eachother. Their interaction
// is investigated under different viscosity conditions.

/**
 * @brief Function adds a Gaussian to a given field @p field
 *
 * This function adds a two dimensional Gaussian blob to a field @p field
 * The Gaussian's location is given by @p x_center, @p y_center, its amplitude
 * by @p A and its standard deviation i.e. width by @p sigma. 
 *
 * NOTE THAT THIS FUNCTION IS ADDITIVE, meaning that it does not override anything
 * previously stored in the field @p field.
 *
 * @param field The field to which the Gaussian is added.
 * @param N The logical grid size (N * N)
 * @param x_center Location of the x-center of the Gaussian blob on [0,1]
 * @param y_center Location of the y-center of the Gaussian blob on [0,1]
 * @param A Amplitude of the Gaussian (determines e.g. dye field density or force strength)
 * @param sigma The standard deviation for the Gaussian. Taken as a float, not cov matrix.
 * 
 */
void add_gaussian(float *field, int N, float x_center, float y_center, float A, float sigma) {
    for (int i = 0; i < N; i++) {
        float x = (i + 0.5) / N;
        for (int j = 0; j < N; j++) {
            float y = (j + 0.5) / N;
            float r2 = (x-x_center)*(x-x_center) + (y-y_center)*(y-y_center);
            field[i + N * j] += A * expf(-r2 / (2 * sigma * sigma));
        }
    }
}

/**
 * @brief Function defining the main simulation loop for the viscosity investigation.
 *
 * This function takes in the viscosity @p visc as an argument, and simulates the interaction
 * of four Gaussian blobs under the behavior of a Gaussian force applied to each blob which sends
 * the blobs towards eachother at the start of the simulation, in a fluid with the given viscosity.
 *
 * Simulation specific parameters e.g. the timestep or grid size can be assigned inside of the
 * simulation function itself.
 *
 * @param visc A floating point number defining the kinematic viscosity used in the simulation
 * @param base_filename A pointer to a character array storing the base of the filename for storing the simulation output
 * 
 * Outputs are stored in the /outputs directory with the base filename, followed by the date and time of the simulation,
 * e.g. "viscosity_investigation_visc_0.0100_19_05_2025_17h42m10s.bin"
 */
int simulation(float visc, char *base_filename) {
    // Specify simulation parameters 
    const int N = 500;
    const int total_steps = 200;
    float dt = 0.01;
    // visc defined in function call

    // Blob parameter
    const float blob1_x = 0.2;
    const float blob1_y = 0.2;

    const float blob2_x = 0.2;
    const float blob2_y = 0.8;

    const float blob3_x = 0.8;
    const float blob3_y = 0.8;

    const float blob4_x = 0.8;
    const float blob4_y = 0.2;

    const float blob_sigma = 0.05;

    // Set this to 1 and it will render to terminal (zoom out in terminal)
    // Otherwise an output binary file is created which can be read in
    // from external software for animation.
    int render_to_terminal = 1; 
    
    char* filename = NULL;
    if (render_to_terminal != 1) {
        char* filename = init_binary_output_file(N, N, total_steps, base_filename);
    }

    Grid* grid = grid_create(N);
    // From grid we have access to
    // u, v, u0, v0, D, D0

    // Necessary for initializing solver (could abstract away to some VelocitySolver object, which needs init.?)
    // Also would benefit from renaming to something like init_velocity_solver(N)
    init_FFT(N); 

    // Add Gaussian blobs to density field
    add_gaussian(grid->D, N, blob1_x, blob1_y, 0.5, blob_sigma);
    add_gaussian(grid->D, N, blob2_x, blob2_y, 0.5, blob_sigma);
    add_gaussian(grid->D, N, blob3_x, blob3_y, 0.5, blob_sigma);
    add_gaussian(grid->D, N, blob4_x, blob4_y, 0.5, blob_sigma);

    // Simulation loop
    for (int step = 0; step < total_steps; step++) {
        // if (step % 100 == 0) printf("Simulation step %i completed\n", step);
        if (step < 100) {
            // We need to zero the forces because the add_gaussian function is additive
            // so we would be adding to the previous velocities (since the force arrays
            // are repurposed for velocity arrays in the solver)
            memset(grid->f_x, 0, sizeof(float) * N * (N + 2));
            memset(grid->f_y, 0, sizeof(float) * N * (N + 2));

            add_gaussian(grid->f_x, N, blob1_x, blob1_y, 5.0, 1.5*blob_sigma);
            add_gaussian(grid->f_x, N, blob2_x, blob2_y, 5.0, 1.5*blob_sigma);
            add_gaussian(grid->f_x, N, blob3_x, blob3_y, -5.0, 1.5*blob_sigma);
            add_gaussian(grid->f_x, N, blob4_x, blob4_y, -5.0, 1.5*blob_sigma);

            add_gaussian(grid->f_y, N, blob1_x, blob1_y, 5.0, 1.5*blob_sigma);
            add_gaussian(grid->f_y, N, blob2_x, blob2_y, -5.0, 1.5*blob_sigma);
            add_gaussian(grid->f_y, N, blob3_x, blob3_y, -5.0, 1.5*blob_sigma);
            add_gaussian(grid->f_y, N, blob4_x, blob4_y, 5.0, 1.5*blob_sigma);
            } else {
            // No force
            memset(grid->f_x, 0, sizeof(float) * N * (N + 2));
            memset(grid->f_y, 0, sizeof(float) * N * (N + 2));
        }

        // Solve velocity field, then dye field from 
        // velocity field, then write dye field to file
        evolve_velocity_field(grid->N, grid->u, grid->v, grid->f_x, grid->f_y, visc, dt);
        evolve_dye_field(grid->N, grid->u, grid->v, grid->D, grid->D0, dt);

        if (render_to_terminal == 1) {
            render_colored_field(grid->N, grid->D, 4); // Draw every 4th grid point
            usleep(33000); // Display for 33ms so that terminal output is arond 30FPS
        } else {
            write_array_to_binary(grid->D, step, N, N, filename);
        }
        
    }

    // Free assigned memory
    grid_free(grid);
    deinit_FFT();
    return 0;
}

// Main function which runs 4 simulations each with variable viscosities. 
int main() {
    // Viscosities to be simulated
    float visc_arr[] = {0.1f, 0.01f};
    int N = sizeof(visc_arr) / sizeof(visc_arr[0]);
    char base_filename[100];

    for (int i = 0; i < N; i++) {
        // printf("Simulating with viscosity = %.4f\n", visc_arr[i]);

        // Initial simulation was designed for outputting binary files, so a naming
        // scheme is included here but it is unnecessary for terminal window rendering.
        sprintf(base_filename, "viscosity_investigation_visc_%.4f", visc_arr[i]);
        simulation(visc_arr[i], base_filename); 
    }

    return 0;
}
