Dependencies:
* nvcc - NVIDIA CUDA Compiler
* g++  - any compatible version (4.5 and below)
* Cusp - Available at http://code.google.com/p/cusp-library/

To compile, set the environment variable CUSP_DIR to point to the directory 
with the Cusp library. For a bash shell, add the following line to ~/.bashrc or 
~/.bash_profile:

	export CUSP_DIR=/path/to/cusp_folder

And reload the file:

	source ~/.bashrc
	
Compile all files:
	
	make

Run the test:
	
	bin/cuIBM

Update: September 6th, 2012

ATTENTION: 
This is a work-in-progress code. Only the GPU-based routines work. The available
version solves flows over only one arbitrarily moving body.

The following time-stepping schemes have been tested for the plain 
Navier-Stokes and the Taira & Colonius solvers:

* Convection
	- EULER_EXPLICIT    : First-order Explicit Euler
	- ADAMS_BASHFORTH_2 : Second-order Adams-Bashforth

* Diffusion
	- EULER_EXPLICIT : First-order explicit Euler
	- EULER_IMPLICIT : First order implicit Euler
	- CRANK_NICOLSON : Crank-Nicolson

Examples available in the makefile are:

	- LidDrivenCavityRe100
	- LidDrivenCavityRe1000
	- cylinderRe40
	- cylinderRe550
	- cylinderRe3000
	- flappingRe75

Known issues:

* Linker step fails on Mac OS X
* Cannot choose the Krylov solver.
* TairaColoniusSolver does not calculate forces if multiple bodies are present.
* Forces are not calculated for the direct forcing method (FadlunEtAlSolver)
* FadlunEtAlSolver has not been fully validated
