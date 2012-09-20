Updated: September 20th, 2012

Dependencies:

* nvcc - NVIDIA CUDA Compiler (tested with nvcc-4.2)
* g++  - any compatible version (nvcc-4.2 is compatible with g++4.6. nvcc-4.1 requires g++4.5 or below)
* Cusp - Available at http://code.google.com/p/cusp-library/ (tested with revision 558)

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

ATTENTION: 
This is a work-in-progress code. Only the GPU-based routines work. The available
version solves for flow over an arbitrarily moving body.

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

	LidDrivenCavityRe100
	LidDrivenCavityRe1000
	cylinderRe40
	cylinderRe550
	cylinderRe3000
	flappingRe75

Run
	
	make <examplename>

to run any of the examples. The biggest case requires a 2GB card.
	
As of now, post-processing requires MATLAB. Navigate into the folder
	
	scripts/matlab
	
and run the script with the following options:

	postprocess <foldername> <xmin> <xmax> <ymin> <ymax>

Proper python/gnuplot scripts for the validation cases will be included later.

Known issues:

* Linker step fails on Mac OS X
* Cannot choose the Krylov solver.
* TairaColoniusSolver and FadlunEtAlSolver fail if no body is present.
* TairaColoniusSolver does not calculate forces if multiple bodies are present.
* Forces are not calculated for the direct forcing method (FadlunEtAlSolver)
* FadlunEtAlSolver has not been fully validated
