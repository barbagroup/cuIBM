cuIBM - A GPU-based Immersed Boundary Method code
=================================================

**ATTENTION: The primary respository for cuIBM shall henceforth be hosted on 
[GitHub](https://github.com/barbagroup/cuIBM). However, we will continue to 
push changes to both Github and the 
[old Bitbucket repository](https://bitbucket.org/anushk/cuibm) and keep the 
code synchronised.**

Currently, cuIBM runs only on Unix-based systems and has been tested on 
Ubuntu 12.04. It is not supported on Windows.

Installation instructions
-------------------------

### Dependencies

Please ensure that the following dependencies are installed before compiling 
cuIBM:

* Git distributed version control system (`git`)
* GNU C++ Compiler(`g++`)
* NVIDIA's CUDA Compiler (`nvcc`)
* CUSP Library (available [here](https://github.com/cusplibrary/cusplibrary))

#### Git (`git`)

Install `git`. On Ubuntu, this can be done via the Terminal using the following 
command:

    > sudo apt-get install git-core

#### GNU C++ Compiler (`g++`)

Install `g++` using the following command:

    > sudo apt-get install g++

Check the version of G++ installed:

    > g++ --version

Other development and version control tools can be installed by following the 
instructions under Step 1 in the 
[CompilingEasyHowTo](https://help.ubuntu.com/community/CompilingEasyHowTo) page 
on the Ubuntu Community Help Wiki. Software developers will find it useful to 
install all of them.

#### NVIDIA's CUDA Compiler (`nvcc`)

[Download and install](https://developer.nvidia.com/cuda-downloads) the CUDA 
Toolkit.

Check the version of NVCC installed:

    > nvcc --version

cuIBM has been compiled and tested with NVCC versions 4.1 to 6.0.

**IMPORTANT**: `nvcc-4.1` is compatible only with G++ version 4.5 (`g++-4.5`) 
or below. `nvcc-4.2` and above are compatible with `g++-4.6` and below.

#### CUSP Library

CUSP is a library that provides routines to perform sparse linear algebra 
computations on Graphics Processing Units. It is written using the CUDA 
programming language and runs code on compatible NVIDIA GPUs. 

CUSP is currently hosted on 
[GitHub](https://github.com/cusplibrary/cusplibrary). cuIBM has been tested 
and works with version 0.4.0, available for download 
[here](https://github.com/cusplibrary/cusplibrary/archive/0.4.0.zip).

The instructions here assume that the CUSP library is to be installed in the 
folder `$HOME/lib`, but any other folder with write permissions can be used. 
Create a local copy of the CUSP library using the following commands:

    > mkdir -p $HOME/lib
    > cd $HOME/lib
    > wget https://github.com/cusplibrary/cusplibrary/archive/0.4.0.zip
    > unzip 0.4.0.zip

The folder `$HOME/lib/cusplibrary-0.4.0` is now created.

If you wish to use to latest version of CUSP, it can be cloned from the GitHub 
repository.

    > cd $HOME/lib
    > git clone https://github.com/cusplibrary/cusplibrary.git

which creates the folder `$HOME/lib/cusplibrary`.

### Compiling cuIBM

cuIBM is currently hosted on both [Github](https://github.com/barbagroup/cuIBM) 
and [BitBucket](https://bitbucket.org/anushk/cuibm). The Github repository is 
the official one, but both will reflect the latest changes.

Run the following commands to create a local copy of the repository in the 
folder `$HOME/src` (or any other folder with appropriate read/write/execute 
permissions):

    > mkdir -p $HOME/src
    > cd $HOME/src
    > git clone https://github.com/barbagroup/cuIBM.git

To compile, set the environment variable `CUSP_DIR` to point to the directory 
with the CUSP library. For a `bash` shell, add the following line to the file 
`~/.bashrc` or `~/.bash_profile`:

    export CUSP_DIR=/path/to/cusp/folder

which for the present case would be `$HOME/lib/cusplibrary-0.4.0`.

We also recommend setting the environment variable `CUIBM_DIR` to point to the 
location of the cuIBM folder. While the code can be compiled and run without 
setting this variable, some of the validation scripts provided make use of it.

    export CUIBM_DIR=/path/to/cuIBM/folder
    
which is `$HOME/src/cuIBM`, as per the above instructions.

Reload the file:

    > source ~/.bashrc

Switch to the cuIBM directory:

    > cd $HOME/src/cuIBM

or if you've set the environment variable,

    > cd $CUIBM_DIR

Compile all the files:
    
    > make

Run the test:
    
    > bin/cuIBM
    
**IMPORTANT**: If your NVIDIA card supports only compute capability 1.3, then 
edit Line 13 of the file `Makefile.inc` in the cuIBM root directory before 
compiling: replace `compute_20` with `compute_13`.

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

The convection terms are calculated using a conservative symmetric 
finite-difference scheme, and the diffusion terms are calculated using a 
second-order central difference scheme.

Examples
--------

The following are available in the default installation:

* `lidDrivenCavityRe100`: Flow in a lid-driven cavity with Reynolds number 
100.
* `lidDrivenCavityRe1000`: Flow in a lid-driven cavity with Reynolds number 
1000.
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

### Run the tests
    
To run any of the examples:

    > make <examplename>

The biggest case (`cylinderRe3000`) requires a graphics card with 2GB of memory.

Post-processing
---------------

The only currently available post-processing script is 
`$CUIBM_DIR/scripts/python/plotVorticity.py`. It plots the vorticity field of 
the flow at all the saved time steps. To display a list of all the command line 
options (which include the case folder, and the coordinates of the corners of 
the region of interest), run:

    > python $CUIBM_DIR/scripts/python/plotVorticity.py --help

To obtain the vorticity plots, navigate to a case folder (or specify it using 
the command line option) and run the script:

    > python $CUIBM_DIR/scripts/python/plotVorticity.py

Known issues
------------

* CPU routines do not work.
* Cannot specify which Krylov solver to use for solving the linear systems.
* TairaColoniusSolver and DirectForcingSolver fail if no body is present.
* DirectForcingSolver has not been tested for cases with multiple or moving 
bodies.

Contact
-------

Please e-mail [Anush Krishnan](mailto:k.anush@gmail.com) if you have any 
questions, suggestions or feedback.
