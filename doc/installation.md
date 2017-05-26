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


### Contact

Please e-mail [Anush Krishnan](mailto:k.anush@gmail.com), 
or [Olivier Mesnard](mailto:mesnardo@gwu.edu), if you have any 
questions, suggestions or feedback.

To report bugs, please use the GitHub issue tracking system. We are also open to pull-requests.
