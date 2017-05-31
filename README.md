# cuIBM - A GPU-based Immersed Boundary Method code

cuIBM solves the 2D incompressible Navier-Stokes equations with an immersed-boundary method on a single CUDA-capable GPU device with the [CUSP](https://github.com/cusplibrary/cusplibrary) library.

Currently, cuIBM runs only on Unix-based systems (no support on Windows) and was last tested on Ubuntu 16.04.


## Installation

### Dependencies (last tested)

* GNU C++ Compiler(`g++-5.4`)
* NVIDIA's CUDA Compiler (`nvcc-8.0`)
* [CUSP](https://github.com/cusplibrary/cusplibrary) (`0.5.1`)

#### GNU C++ Compiler

Install `g++` using the following command:

    > sudo apt-get install g++

Check the version of G++ installed:

    > g++ --version

Other development and version control tools can be installed by following the 
instructions under Step 1 in the 
[CompilingEasyHowTo](https://help.ubuntu.com/community/CompilingEasyHowTo) page 
on the Ubuntu Community Help Wiki. Software developers will find it useful to 
install all of them.

#### NVIDIA's CUDA Compiler

[Download and install](https://developer.nvidia.com/cuda-downloads) the CUDA 
Toolkit.

Check the version of NVCC installed:

    > nvcc --version

#### CUSP Library

CUSP is a library that provides routines to perform sparse linear algebra 
computations on Graphics Processing Units. It is written using the CUDA 
programming language and runs code on compatible NVIDIA GPUs. 

CUSP is currently hosted on 
[GitHub](https://github.com/cusplibrary/cusplibrary). cuIBM has been tested 
and works with version 0.5.1, available for download 
[here](https://github.com/cusplibrary/cusplibrary/archive/v0.5.1.tar.gz).

    > mkdir -p $HOME/software/cusp/0.5.1
    > wget https://github.com/cusplibrary/cusplibrary/archive/v0.5.1.tar.gz
    > tar -xvf v0.5.1.tar.gz -C $HOME/software/cusp/0.5.1 --strip-components=1
    > rm -f v0.5.1.tar.gz


### Compiling cuIBM

Clone cuIBM:

    > mkdir -p $HOME/software
    > cd $HOME/software
    > git clone https://github.com/barbagroup/cuIBM.git

To compile cuIBM, set the environment variable `CUSP_DIR` to point to the directory with the CUSP library.

    > export CUSP_DIR=$HOME/software/cusp/0.5.1

We also recommend setting the environment variable `CUIBM_DIR` to point to the 
location of the cuIBM folder. While the code can be compiled and run without 
setting this variable, some of the validation scripts provided make use of it.

    > export CUIBM_DIR=$HOME/software/cuIBM

To compile cuIBM, navigate to the cuIBM directory:

    > cd $HOME/software/cuIBM

or

    > cd $CUIBM_DIR

then, run make:

    > make

Note: the current Makefile compiles cuIBM with the flag `-arch=sm_35` (compute capability of 3.5). You should adapt the Makefile if your NVidia card has a different compute capability.

Run an example (2D flow around cylinder at Re=40) to check that cuIBM has been built correctly:

    > bin/cuibm -directory examples/cylinder/Re40

Finally, you can add cuIBM to your path:

    > export PATH="$HOME/software/cuIBM/bin:$PATH"


## Examples

The following examples are available:

* `lidDrivenCavityRe100`: Flow in a lid-driven cavity with Reynolds number 100.
* `lidDrivenCavityRe1000`: Flow in a lid-driven cavity with Reynolds number 1000.
* `cylinderRe40`: Flow over a circular cylinder at Reynolds number 40. The 
flow eventually reaches a steady state.
* `cylinderRe75`: Flow over a circular cylinder at Reynolds number 75. The 
initial flow field has an asymmetric perturbation that triggers instability in 
the flow and vortex shedding is observed in the wake.
* `cylinderRe100`: Flow over a circular cylinder at Reynolds number 100. The 
initial flow field has an asymmetric perturbation that triggers instability in 
the flow and vortex shedding is observed in the wake.
* `cylinderRe150`: Flow over a circular cylinder at Reynolds number 150. The 
initial flow field has an asymmetric perturbation that triggers instability in 
the flow and vortex shedding is observed in the wake.
* `cylinderRe550`: Initial flow over an impulsively started cylinder at 
Reynolds number 550.
* `cylinderRe3000`: Initial flow over an impulsively started cylinder at 
Reynolds number 3000.
* `flappingRe75`: Flow around a flapping airfoil.
* `oscillatingCylinders`: Flow across two oscillating cylinders.
    
To run any of the examples list above:

    > make <examplename>

The biggest case (`cylinderRe3000`) requires a graphics card with 2GB of memory.

Post-processing scripts are available in th directory 
`$CUIBM_DIR/scripts/python`. The command-line argument `--help` will display 
the list of options for the script executed. 

For example, `plotVorticity.py` plots the contour of the vorticity field of the 
flow. To display its list of all the command-line options, run:

    > python $CUIBM_DIR/scripts/python/plotVorticity.py --help


## Known issues

* CPU routines do not work.
* DirectForcingSolver has not been tested for cases with multiple or moving 
bodies.

## Contact

Please e-mail [Anush Krishnan](mailto:k.anush@gmail.com), 
or [Olivier Mesnard](mailto:mesnardo@gwu.edu), if you have any 
questions, suggestions or feedback.

To report bugs, please use the GitHub issue tracking system. We are also open to pull-requests.
