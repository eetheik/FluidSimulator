// src/stable_solve.h
#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

typedef struct fftw_plan_struct *fftwf_plan;

void init_FFT(int n);
void deinit_FFT();
void apply_forces(int n, float *u, float *v, float *u0, float *v0, float dt);
void advect_scalar_field(int n, float *u, float *v, float *DST, float *SRC, float dt);
void FS_diffuse_and_project(int n, float *FSu0, float *FSv0, float visc, float dt);
void evolve_velocity_field(int n, float *u, float *v, float *u0, float *v0, float visc, float dt);
void evolve_dye_field(int n, float *u, float *v, float *D, float *D0, float dt);

#endif