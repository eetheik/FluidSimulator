# FluidSimulator
A 2-dimensional incompressible fluid with periodic boundary conditions is animated using an adaption of Jos Stam "stable-solve" algorithm with visual output to the termial window. The algorithm is detailed in https://pages.cs.wisc.edu/~chaol/data/cs777/stam-stable_fluids.pdf. The algorithm solves the incompressible Navier-Stokes equations

$$ \frac{\partial \vec{u}}{\partial t} + (\vec{u}\cdot\nabla)\vec{u} = -\frac{1}{\rho}\nabla p + \nu \nabla ^2\vec{u} + \vec{f}, \quad \text{and} \quad \nabla \cdot \vec{u} = 0,$$

with $\vec{u}$ a fluid velocity field, $t$ time, $\rho$ the density, $p$ the pressure, $\nu$ the kinematic viscosity and $\vec{f}$ any external force such as gravity. We take advantage of the fact that any vector field $\vec{w}$ can be decomposed into a divergence-free (i.e. $\nabla \cdot \vec{u} = 0$) part and the gradient of a scalar field $q$, as

$$ \vec{w} = \vec{u} + \nabla q. $$

One can then devise a projection operator $P$, which acts as $\vec{u}=P\vec{w} = \vec{w}-\nabla q$, (and consequently $P\vec{u} = \vec{u}, P\nabla p = 0$) so that the incompressible Navier-Stokes equations become

$$ \frac{\partial \vec{u}}{\partial t} = P\left(-(\vec{u} \cdot \nabla)\vec{u} + \nu\nabla^2 \vec{u} + \vec{f}\right). $$

This code implements periodic boundary conditions which enables the easy handling of the diffusion $(\nu \nabla^2 \vec{u})$ and projection $(P(\cdots))$ terms in the Fourier domain. Self-advection $(\vec{u}\cdot\nabla\vec{u})$ is handled in an unconditionally stable manner by matching the velocity of a fluid parcel to its biliniearly interpolated velocity a timestep $dt$ in the past, allowing the simulator to produce realistic animations with arbitrarily large timesteps (see paper). Forces are handled by a very simple Euler scheme. The simulation results are made tangible to a human viewer by adapting the self-advection part of the stable-solve algorithm to advect a scalar dye field using the velocity field. This is exactly analogous to observing food coloring get transported in a fluid. The scalar dye field is rendered to the terminal using ANSI color codes. 

The asymptotic time complexity for the velocity evolution is $O(N \log N)$ due to the FFT, which is inferior to the theoretically optimal $O(N)$ where each fluid parcel is consulted once, which is achievable using multigrid-solvers. 

# Examples
Two example simulations created with the FluidSimulator are shown below.

https://github.com/user-attachments/assets/81cbf7d2-4ea9-4d54-9d5a-fd7b3f9b5c43

https://github.com/user-attachments/assets/701a9dae-bbbe-44c5-ad5a-78bd5316f08b

# Requirements
- C, gcc (I used ggc13)
	- fftw3f (for Fast Fourier Transforms, note that the code uses floats) + standard C libraries
- Make utility
- Zooming out in your terminal window

# Creating a fluid simulation
Fluid simulations can be created in the `/run` directory. Two example simulations:
`kelvin_helmholtz_instability.c` and `viscosity_investigation.c` are provided. The
general outline for the simulation main function is as follows:
- Create the simulation grid using *grid_create* from `src/grid.h`
- Initialize the FFT using *init_FFT* from `src/fluid_solver.h`
- (Optional) Initialize a binary output file using *init_binary_output_file* from `src/io.h`
- Define the main simulation loop:
	- Apply forces
	- Evolve velocity fields using *evolve_velocity_field* from `src/fluid_solver.h`
	- Evolve dye field(s) using *evolve_dye_field* from `src/fluid_solver.h`
	- Render to terminal or write density field to binary using *render_to_terminal* or *write_array_to_binary* from `src/io.h`

# Running a fluid animation
In order to compile the fluid simulation file a Makefile is provied.
It is used as:

`make SIM=<filename>.cc`,

which compiles to a `./simulate` executable. The user should then zoom out in their
terminal window and run the executable

`./simulate`,

to render the simulation to the terminal.

# Code expansion ideas
- Command line inputs during a simulation, e.g. hotkeys for applying certain types of forces or adding dyes
- 3-dimensional fluid simulations (would not be possible in terminal anymore, or maybe ASCII characters instead of colored spaces, like the rotating donut)
