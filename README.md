cuIBM - A GPU-based Immersed Boundary Method code
=================================================

*Updated: June 10th, 2013*

Currently, cuIBM runs only on Linux and has been tested on Ubuntu 11.10. It is not supported on Windows and Mac OS X.

Installation instructions
-------------------------

### Dependencies

Please ensure that the following dependencies are installed before compiling cuIBM:

* Mercurial revision control tool (`hg`)
* GNU C++ Compiler(`g++`)
* NVIDIA's CUDA Compiler (`nvcc`)
* CUSP Library (available [here](http://code.google.com/p/cusp-library/))

#### Mercurial (`hg`)

Install the package `mercurial`. On Ubuntu, this can be done via the Terminal using the following command:

    > sudo apt-get install mercurial

#### GNU C++ Compiler (`g++`)

Install `g++` using the following command:

    > sudo apt-get install g++

Check the version of G++ installed:

    > g++ --version

Other development and version control tools can be installed by following the instructions under Step 1 in the [CompilingEasyHowTo](https://help.ubuntu.com/community/CompilingEasyHowTo) page on the Ubuntu Community Help Wiki. Software developers will find it useful to install all of them.

#### NVIDIA's CUDA Compiler (`nvcc`)

[Download and install](https://developer.nvidia.com/cuda-downloads) the CUDA Toolkit.

Check the version of NVCC installed:

    > nvcc --version

cuIBM has been tested and compiles with NVCC versions 4.1, 4.2 and 5.0.

**IMPORTANT**: `nvcc-4.1` is compatible only with G++ version 4.5 (`g++-4.5`) or below. `nvcc-4.2` and `nvcc-5.0` are compatible with `g++-4.6` and below.

#### CUSP Library

CUSP is a library that provides routines to perform sparse linear algebra computations on Graphics Processing Units. It is written using the CUDA programming language and runs code on compatible NVIDIA GPUs. The version of CUSP that is currently used by cuIBM is hosted on [Google Code](http://code.google.com/p/cusp-library/).

The instructions here assume that the CUSP library is installed in a folder called `lib` in the home directory, but any other folder with write permissions can be used. Create a local copy of the CUSP library using the following commands:

    > mkdir -p $HOME/lib
    > cd $HOME/lib
    > hg clone https://code.google.com/p/cusp-library/

The above command will clone the latest revision of CUSP hosted on Google Code, i.e. revision 558. This is the version that cuIBM has been tested with and works.

The CUSP library is no longer being actively maintained on the Google Code website. It is currently hosted on [GitHub](https://github.com/cusplibrary/cusplibrary).

### Compiling cuIBM

cuIBM is currently hosted on [BitBucket](https://bitbucket.org/anushk/cuibm). Run the following commands to create a local copy in the folder `src` in the home directory (any other folder with appropriate read/write permissions can also be used):

    > mkdir -p $HOME/src
    > cd $HOME/src
    > hg clone https://bitbucket.org/anushk/cuibm

To compile, set the environment variable `CUSP_DIR` to point to the directory with the CUSP library. For a `bash` shell, add the following line to the file `~/.bashrc` or `~/.bash_profile`:

    export CUSP_DIR=/path/to/cusp/folder

which for the above case would be `$HOME/lib/cusp-library`

And reload the file:

    > source ~/.bashrc

Switch to the cuIBM directory:

    > cd /path/to/cuibm/folder

The path as per the above instructions is `$HOME/src/cuibm`.

Compile all the files:
	
    > make

Run the test:
	
    > bin/cuIBM
	
**IMPORTANT**: If your NVIDIA card supports only compute capability 1.3, then edit Line 13 of the Makefile in the cuIBM root directory before compiling: replace `compute_20` with `compute_13`.

Numerical schemes
-----------------

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

The convection terms are calculated using a conservative symmetric finite-difference scheme, and the diffusion terms are calculated using a second-order central difference scheme.

Examples
--------

The following are available in the default installation:

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

The biggest case (`cylinderRe3000`) requires a graphics card with 2GB of memory.

Post-processing
---------------

Currently, post-processing requires MATLAB. Navigate into the folder
	
	cd $CUIBM_DIR/scripts/matlab
	
and run the script `postprocess.m` with the following options:

	postprocess <foldername> <xmin> <xmax> <ymin> <ymax>

`python`/`gnuplot` scripts for the validation cases will be included soon.

Known issues
------------

* CPU routines do not work for FadlunEtAlSolver.
* SuLaiLinSolver does not work.
* Cannot specify which Krylov solver to use for solving the linear systems.
* TairaColoniusSolver and FadlunEtAlSolver fail if no body is present.
* Individual forces are not calculated if multiple bodies are present.
* FadlunEtAlSolver has not been tested for cases with multiple or moving bodies.

Contact
-------

Please e-mail [Anush Krishnan](mailto:anush@bu.edu) if you have any questions, suggestions or feedback.
