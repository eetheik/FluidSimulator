# FluidSimulator
A fluid simulator which outputs fluid animations to the terminal. 

A 2-dimensional incompressible fluid with periodic boundary conditions is animated
using an adaption of Jos Stam "stable-sovle" algorithm with visual output to the 
termial window. The algorithm 

# Requirements
- C, gcc (I used ggc13)
	- clib
	- fftw3 (for Fast Fourier Transforms)
- Make utility
- Zooming out in your terminal window

# Creating a fluid simulation
Fluid simulations can be created in the "/run" directory. Two example simulations:
"kelvin_helmholtz_instability.c" and "viscosity_investigation.c" are provided. The
general outline for the simulation main function is as follows:
	1) Create the simulation grid using grid_create(N) from src/grid.h
	2) Initialize the FFT using init_FFT(N) from src/fluid_solver.h
	(Optional) Initialize a binary output file using init_binary_output_file(...) from src/io.h
	3) Define the main simulation loop:
		- Apply forces
		- Evolve velocity fields using evolve_velocity_field(...) from src/fluid_solver.h
		- Evolve dye field(s) using evolve_dye_field(...) from src/fluid_solver.h
		- Render to terminal or write density field to binary using render_to_terminal(...) or write_array_to_binary from src/io.h

# Running a fluid animation
In order to compile the fluid simulation file a Makefile is provied.
It is used as:

"make SIM=<filename>.cc",

which compiles to a "./simulate" executable. The user should then zoom out in their
terminal window and run

"./simulate",

to render the simulation to the terminal.

# Code expansion ideas
- Command line inputs during a simulation, e.g. hotkeys for applying certain types of forces or adding dyes
- 3-dimensional fluid simulations (would not be possible in terminal anymore, or maybe ASCII characters instead of spaces, like the rotating donut)
