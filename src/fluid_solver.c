#include <fftw3.h>
#include <math.h>
#include <string.h>

// Global variables for FFT plans
static fftwf_plan plan_rc, plan_cr;

/**
 * @brief Initializes real-to-complex (r2c) and complex-to-real (c2r) FFTs
 *
 * This functions creates FFTW plans for performing forward (r2c) and inverse (c2r)
 * 2D Fast Fourier Transforms on square arrays of size n * n. The input and output
 * arrays are not specified at initialization (NULL), allowing the plan to be created 
 * without binding to a specific memory location. 
 *
 * The plans for the forward and backward FFTs are stored in the global static variables
 * plan_rc and plan_cr respectively. The `FFTW_ESTIMATE` flag denotes that a simple 
 * heuristic is used to pick a (probably sub-optimal) plan quickly. For details see
 * https://www.fftw.org/fftw3_doc/Planner-Flags.html, and https://www.fftw.org/fftw3_doc/Using-Plans.html.
 *
 * @param n The logical size of the square array (n * n) to be transformed.
 */
void init_FFT(int n)
{
  plan_rc = fftwf_plan_dft_r2c_2d(n, n, NULL, NULL, FFTW_ESTIMATE);
  plan_cr = fftwf_plan_dft_c2r_2d(n, n, NULL, NULL, FFTW_ESTIMATE);
}


/**
 * @brief De-initializes the real-to-complex (r2c) and complex-to-real (c2r) FFT plans
 *
 * This function deallocates the FFT plans stored in the global static variables 
 * plan_rc and plan_cr.
 */
void deinit_FFT()
{
  fftwf_destroy_plan(plan_rc);
  fftwf_destroy_plan(plan_cr);
}


/**
 * @brief Executes a forward or inverse FFT (dependent on 's') on a given array 'u' using FFTW
 *
 * This macro acts as a wrapper to execute pre-computed FFTW plans. When using
 * s == 1 a forward (real-to-complex) FFT is performed and otherwise the inverse
 * (complex-to-real) FFT is performmed. The transforms are executed on an in-place
 * float or fftwf_complex array pointer 'u'.
 *
 * The macro requires that the FFTW plans 'plan_rc' and 'plan_cr' have already been 
 * initialized using init_FFT(). The input array 'u' is re-interpreted as a
 * float or fftwf_complex array depending on the direction of the FFT.
 * 
 * @param s An integer flag: set to 1 for forward FFT (r2c) and anything else for inverse FFT (c2r).
 * @param u Pointer to the float or fftwf_complex array that is to be transformed. The memory for 
 *          the array needs to have been correctly assigned, i.e. with fftwf_malloc().
 *
 */
#define FFT(s, u)                                                   \
  if (s == 1)                                                       \
    fftwf_execute_dft_r2c(plan_rc, (float *)u, (fftwf_complex *)u); \
  else                                                              \
    fftwf_execute_dft_c2r(plan_cr, (fftwf_complex *)u, (float *)u)


/**
 * @brief Applies forces input via @p u0, @p v0 to velocity field @p u, @p v.
 *
 * This function is used to apply forces input via @p u0, @p v0 to velocity field @p u, @p v.
 * Forces are applied simply as: (u, v) = (u, v) + dt * (f_x, f_y).
 * 
 * @param n: Grid size (n * n)
 * @param u: A pointer to the float array storing the x-components of the velocity field
 * @param v: A pointer to the float array storing the y-components of the velocity field
 * @param u0: A pointer to a float array storing the x-component of the force field 
 *            (suggestively named since u0, v0 will be repurposed as velocity buffers 
 *            following Jos Stam's implementation)
 * @param v0: A pointer to a float array storing the y-component of the force field
 * @param dt: Timestep
 *
 */
 void apply_forces(int n, float *u, float *v, float *u0, float *v0, float dt) {
    for (int i = 0; i < n * n; i++) {
        u[i] += dt * u0[i];
        // u0[i] = u[i]; // Handle copying elsewhere. A function called "apply_forces" should not silently copy into u0, v0
        v[i] += dt * v0[i];
        // v0[i] = v[i];
    }
 }

/**
 * @brief Advects a scalar field by a vector field (u, v).
 *
 * This function is used to advect scalar field values in @p IN by a velocity field
 * defined by @p u, @p v. Results are stored to @p OUT. Advection is done by backtracing
 * the velocity field a timestep @p dt, and assigning @p OUT the bilinearly interpolated
 * value of the @p IN field at the backtraced location.
 * 
 * @param n: Grid size (n * n)
 * @param u: A pointer to the float array storing the x-components of the velocity field
 * @param v: A pointer to the float array storing the y-components of the velocity field
 * @param OUT: A pointer to a float array storing the output advected scalar field
 * @param IN: A pointer to a float array storing the input scalar field to be advected
 * @param dt: Timestep
 *
 */
void advect_scalar_field(int n, float *u, float *v, float *OUT, float *IN, float dt) {
    
    float x, y, x0, y0, s, t;
    int i, j, i0, j0, i1, j1;

    for (x = 0.5 / n, i = 0; i < n; i++, x += 1.0 / n) {

        for (y = 0.5 / n, j = 0; j < n; j++, y += 1.0 / n) {

            // Backtrace the i,j'th voxel back a timestep dt
            // and determine the indices i0, j0 of the departing voxel
            x0 = n * (x - dt * u[i + n * j]) - 0.5; // This gives a floating point value in [0, N]
            y0 = n * (y - dt * v[i + n * j]) - 0.5;
            i0 = floor(x0);
            s = x0 - i0;
            i0 = (n + (i0 % n)) % n; // mod n because we have periodic B.C. and this might wrap around the edge
            i1 = (i0 + 1) % n;

            j0 = floor(y0);
            t = y0 - j0;
            j0 = (n + (j0 % n)) % n;
            j1 = (j0 + 1) % n; 

            // Linearly interpolate the scalar field from the backtraced 
            // points (i0, j0), (i0, j1) and (i1, j0), (i1, j1)
            OUT[i + n * j] = (1 - s) * ((1 - t) * IN[i0 + n * j0] + t * IN[i0 + n * j1]) +
                           s * ((1 - t) * IN[i1 + n * j0] + t * IN[i1 + n * j1]);
        }
    }

 }


 /**
 * @brief Applies the diffusion and projection operators in Fourier Space (FS) on arrays FSu0, FSv0
 *
 * This function is used to apply the diffusion and projection operators while in Fourier Space.
 * The function acts on the Fourier Transformed arrays @p FSu0, @p FSv0 by applying
 * a low-pass filter exp(-k^2 * @p dt * @p visc ) in order to treat diffusion, followed by
 * a projection onto the divergence free part of the velocity fields in the Fourier space. See
 * report or https://pages.cs.wisc.edu/~chaol/data/cs777/stam-stable_fluids.pdf for mathematical details.
 * 
 * @param n: Logical grid size (n * n)
 * @param FS_u0: A pointer to the fftwf_complex array storing the x-components of the Fourier transformed velocity field
 * @param FS_v0: A pointer to the fftwf_complex array storing the y-components of the Fourier transformed velocity field
 * @param visc: Kinematic viscosity
 * @param dt: Timestep
 *
 */
 void FS_diffuse_and_project(int n, float *FS_u0, float *FS_v0, float visc, float dt) {
    float U[2], V[2];

    // Take steps of length 2 since we have complex numbers in FS, and they are stored as
    // consecutive elements (real1, imaginary1, real2, imaginary2, ...)
    for (int i = 0; i <= n; i += 2) {

        float x = 0.5 * i;

        for (int j = 0; j < n; j++) {

            int l = j;
            if (j > n / 2) {
                l = j - n;
            }

            float r = x * x + l * l; // wave number squared
            if (r == 0.0) continue;
            
            // viscosity term
            float f = exp(-r * dt * visc);

            // Fourier transformed velocities
            U[0] = FS_u0[i + (n + 2) * j];
            V[0] = FS_v0[i + (n + 2) * j];
            U[1] = FS_u0[i + 1 + (n + 2) * j];
            V[1] = FS_v0[i + 1 + (n + 2) * j];

            // Projections and applying viscosity multiplication term (f)
            // In Fourier space the projection reads 
            // \hat{P}\hat{w}(k) = \hat{w}(k) - (1/k^2) * (k \cdot \hat{w}(k))k
            // Writing what this becomes for each component of w(k) = (u(k), v(k)) 
            // you obtain the equations below:
            FS_u0[i + (n + 2) * j] = f * ((1 - x * x / r) * U[0] - x * l / r * V[0]);
            FS_u0[i + 1 + (n + 2) * j] = f * ((1 - x * x / r) * U[1] - x * l / r * V[1]);

            FS_v0[i + (n + 2) * j] = f * (-l * x / r * U[0] + (1 - l * l / r) * V[0]);
            FS_v0[i + 1 + (n + 2) * j] = f * (-l * x / r * U[1] + (1 - l * l / r) * V[1]);
        }
    }
 }
 
 /**
 * @brief A function for evolving the velocity field @p u, @p v through a timestep @p dt with viscosity @p visc
 *
 * This function is used to evolve velocity field @p u, @p v through a timestep @p dt through a fluid with viscosity
 * @p visc. We assume an incompressible fluid. The mathematical details for this solution method can be checked from 
 * the report or https://pages.cs.wisc.edu/~chaol/data/cs777/stam-stable_fluids.pdf. We evolve over the timestep in 3 main stages:
 *  (1) External forces are applied to the velocity fields 
 *  (2) Self-advection is applied to the velocity fields,
 *  (3) A FFT is performed, to apply diffusion and a projection onto the divergence free solution in Fourier space. 
 * Then an inverse FFT and normalization is done to obtain the final solution. Steps (1), (2) and (3) are handled by
 * the functions apply_forces(...) (1), advect_scalar_field(...) (2) applied to each velocity field component (scalar) 
 * alone, and FS_diffuse_and_project(...) (3) with a call to the FFT(...) macro before and after the diffusion&projection, 
 * so that we enter and leave the Fourier space as needed.
 * 
 * @param n: Grid size (n * n)
 * @param u: A pointer to the float array storing the x-components of the velocity field
 * @param v: A pointer to the float array storing the y-components of the velocity field
 * @param u0: A pointer to a float array storing the x-component of the force field 
 *          (u0, v0 are repurposed as temporary buffer arrays to store velocities in the solution)
 * @param v0: A pointer to a float array storing the y-component of the force field
 * @param visc: A floating point number defining the kinematic viscosity of the fluid
 * @param dt: Timestep
 *
 */
void evolve_velocity_field(int n, float *u, float *v, float *u0, float *v0, float visc, float dt) {
    apply_forces(n, u, v, u0, v0, dt);

    // Repurpose force input arrays u0, v0 as velocity buffers
    memcpy(u0, u, sizeof(float) * n * (n + 2));
    memcpy(v0, v, sizeof(float) * n * (n + 2));

    // The velocity field x,y-components are scalar fields themselves, can advect them as such
    advect_scalar_field(n, u0, v0, u, u0, dt); 
    advect_scalar_field(n, u0, v0, v, v0, dt); 

    // Copy u, v into u0, v0, but with padding. Padding is necessary for FFT
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            u0[i + (n + 2) * j] = u[i + n * j];
            v0[i + (n + 2) * j] = v[i + n * j];
        }
    }

    // Fourier transforms for diffusion and projection step
    FFT(1, u0);
    FFT(1, v0);

    FS_diffuse_and_project(n, u0, v0, visc, dt); // FS_... in function name denotes that it gets called in Fourier Space

    FFT(-1, u0);
    FFT(-1, v0);

    // Normalize and copy buffers u0, v0 back to u, v
    float CNORM = 1.0 / (n * n); 
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            u[i + n * j] = CNORM * u0[i + (n + 2) * j]; 
            v[i + n * j] = CNORM * v0[i + (n + 2) * j];
        }
    }
}

 /**
 * @brief A function for evolving the dye field @p D, through a timestep @p dt by the velocity field @p u, @p v
 *
 * This function is used to evolve the dye field @p D by a timestep @p dt by the velocity field @p u, @p v. 
 * The function is essentially a wrapper for the advect_scalar_field(...) function, but the function also
 * does a memcpy of the given dye field @p D to the buffer @p D0, allowing for the buffer @p D0 to be passed in as 
 * an empty array. The buffer is used as the IN input scalar field to the advect_scalar_field(...) function,
 * and the output is stored to @p D = OUT.
 * 
 * @param n: Grid size (n * n)
 * @param u: A pointer to the float array storing the x-components of the velocity field
 * @param v: A pointer to the float array storing the y-components of the velocity field
 * @param D: A pointer to a float array storing the scalar dye field
 * @param D0: A pointer to a float array to be used as a tmp array for @p D. Can be supplied empty, but not unitialized
 * @param dt: Timestep
 *
 */
void evolve_dye_field(int n, float *u, float *v, float *D, float *D0, float dt) {
    memcpy(D0, D, sizeof(float) * n * n); // Populate D0 with current D values
    advect_scalar_field(n, u, v, D, D0, dt);
}
