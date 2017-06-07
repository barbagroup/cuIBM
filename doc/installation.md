### Dependencies (last tested)

* GNU C++ Compiler(`g++-5.4`)
* NVIDIA's CUDA Compiler (`nvcc-8.0`)
* [CUSP](https://github.com/cusplibrary/cusplibrary) (`0.5.1`)
* [Boost](https://www.boost.org) (`1.64.0`)

#### GNU C++ Compiler

Install `g++` using the following command:

    > sudo apt-get install g++

Check the version of g++ installed:

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

    > export CUSP_DIR=$HOME/software/cusp/0.5.1
    > mkdir -p $CUSP_DIR
    > wget https://github.com/cusplibrary/cusplibrary/archive/v0.5.1.tar.gz
    > tar -xvf v0.5.1.tar.gz -C $CUSP_DIR --strip-components=1
    > rm -f v0.5.1.tar.gz

#### Boost library

In cuIBM, we use the parser [YAML-CPP](https://github.com/jbeder/yaml-cpp) (version `0.5.1` bundled in cuIBM) that requires header files from the Boost library.

    > export BOOST_DIR=$HOME/software/boost/1.64.0
    > mkdir -p $BOOST_DIR
    > wget https://dl.bintray.com/boostorg/release/1.64.0/source/boost_1_64_0.tar.gz
    > tar -xvf boost_1_64_0.tar.gz -C $BOOST_DIR --strip-components=1
    > rm -f boost_1_64_0.tar.gz


### Compiling cuIBM

Clone cuIBM:

    > mkdir -p $HOME/software
    > cd $HOME/software
    > git clone https://github.com/barbagroup/cuIBM.git

To compile cuIBM, make sure you have set correctly the environment variables `CUSP_DIR` and `BOOST_DIR` to point to their respective folder:

    > export CUSP_DIR=$HOME/software/cusp/0.5.1
    > export BOOST_DIR=$HOME/software/boost/1.64.0

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
