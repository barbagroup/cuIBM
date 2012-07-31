Dependencies:
nvcc - NVIDIA CUDA Compiler
g++  - any compatible version (4.5 and below)
Cusp - Available at http://code.google.com/p/cusp-library/

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

Update: July 31st, 2012

ATTENTION: 
This is a work-in-progress code. Only the GPU-based routines work. The available 
version solves flows over arbitrary bodies, but options for moving bodies, 
Krylov solvers and preconditioners currently do not work. However, the maximum 
number of iterations and the relative tolerance for the linear solves can be 
set.

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

* Tested
	- LidDrivenCavityRe100
	- LidDrivenCavityRe1000
	- cylinderRe40
	- cylinderRe550

* Untested
	- cylinderRe3000

The bigger grids may take a long time to converge since the preconditioners 
have not yet been implemented.

Known issues:

* Linker step fails on Mac OS X
* Forces are not calculated for the direct forcing method (FadlunEtAlSolver)
* FadlunEtAlSolver has not been fully validated
