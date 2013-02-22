# cuIBM - A GPU-based immersed boundary method code

*Updated: February 21st, 2013*

## Dependencies

* `nvcc` - NVIDIA CUDA Compiler (tested with `nvcc-4.2`)
* `g++`  - any compatible version (`nvcc-4.2` is compatible with `g++-4.6`. `nvcc-4.1` requires `g++-4.5` or below)
* `Cusp` - Available [here](http://code.google.com/p/cusp-library/) (tested with revision 558)

## Installation instructions

To compile, set the environment variable `CUSP_DIR` to point to the directory 
with the Cusp library. For a `bash` shell, add the following line to 
`~/.bashrc` or `~/.bash_profile`:

	export CUSP_DIR=/path/to/cusp_folder

And reload the file:

	source ~/.bashrc

Switch to the cuibm directory:

	cd $CUIBM_DIR

Compile all files:
	
	make

Run the test:
	
	bin/cuIBM

## Numerical schemes

### Temporal discretisation

The following schemes have been tested for the available solvers:

* Convection
	- `EULER_EXPLICIT`: First-order Explicit Euler
	- `ADAMS_BASHFORTH_2`: Second-order Adams-Bashforth

* Diffusion
	- `EULER_EXPLICIT`: First-order explicit Euler
	- `EULER_IMPLICIT`: First order implicit Euler
	- `CRANK_NICOLSON`: Crank-Nicolson

### Spatial discretisation 

The convection terms are calculated using a conservative symmetric 
finite-difference scheme, and the diffusion terms are calculated using a 
second-order central difference scheme.

## Examples

These are available in the default installation:

* `LidDrivenCavityRe100`: Flow in a lid-driven cavity with Reynolds number 
100.
* `LidDrivenCavityRe1000`: Flow in a lid-driven cavity with Reynolds number 
1000.
* `cylinderRe40`: Flow over a circular cylinder at Reynolds number 40. The 
flow eventually reaches a steady state.
* `cylinderRe100`: Flow over a circular cylinder at Reynolds number 100. The 
initial flow field has an asymmetric perturbation that triggers instability in 
the flow and vortex shedding is observed in the wake.
* `cylinderRe150`: Flow over a circular cylinder at Reynolds number 100. The 
initial flow field has an asymmetric perturbation that triggers instability in 
the flow and vortex shedding is observed in the wake.
* `cylinderRe550`: Initial flow over an impulsively started cylinder at 
Reynolds number 550.
* `cylinderRe3000`: Initial flow over an impulsively started cylinder at 
Reynolds number 3000.
* `flappingRe75`: Flow around a flapping airfoil.

### Run the tests
	
To run any of the examples:

	make <examplename>

The biggest case requires a graphics card with 2GB of memory.

## Post-processing	

Currently, post-processing requires MATLAB. Navigate into the folder
	
	cd $CUIBM_DIR/scripts/matlab
	
and run the script `postprocess.m` with the following options:

	postprocess <foldername> <xmin> <xmax> <ymin> <ymax>

`python`/`gnuplot` scripts for the validation cases will be included soon.

## Known issues

* Only GPU-based routines work.
* Linker step fails on Mac OS X
* Cannot choose the Krylov solver.
* TairaColoniusSolver and FadlunEtAlSolver fail if no body is present.
* Individual forces are not calculated if multiple bodies are present.
* FadlunEtAlSolver has not been validated for cases with moving bodies.
