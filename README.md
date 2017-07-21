# cuIBM - A GPU-based Immersed Boundary Method code

[![DOI_JOSS](http://joss.theoj.org/papers/10.21105/joss.00301/status.svg)](http://dx.doi.org/10.21105/joss.00301)
[![CITE_BIB](https://img.shields.io/badge/Cite%20cuIBM-bibtex-blue.svg)](http://www.doi2bib.org/#/doi/10.21105/joss.00301)

cuIBM solves the 2D incompressible Navier-Stokes equations with an immersed-boundary method on a single CUDA-capable GPU device with the [CUSP](https://github.com/cusplibrary/cusplibrary) library.

Currently, cuIBM runs only on Unix-based systems (no support on Windows) and was last tested on Ubuntu 16.04.

---

## Documentation

User's documentation is available on the [Wiki](https://github.com/barbagroup/cuIBM/wiki) pages of the cuIBM repository.

Doxygen API documentation is available [here](http://barbagroup.github.io/cuIBM).

---

## Installation

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

Users-documentation is available in the Wiki pages of the GitHub repository.


### Post-processing

To post-process the numerical solution from cuIBM, we provide Python scripts (present in each case in the folder `examples`).
They make use of the package [`snake`](https://github.com/mesnardo/snake) that is bundled in the `external` folder of cuIBM (version `0.3`).

To install `snake`:

    > cd $CUIBM_DIR/external/snake-0.3
    > python setup.py install

`snake` requires Python (2.7 or 3.5), Matplotlib, Scipy, and Pandas.

---

## Example: Flow over impulsively started cylinder (Re=550)

To run the example:

    > cd $CUIBM_DIR/examples/cylinder/Re550
    > cuibm

To plot the instantaneous drag coefficient and the vorticity field at saved time steps:

    > python scripts/plotDragCoefficient.py
    > python scripts/plotVorticity.py

Figures are saved in the folder `images` of the simulation directory.

We also provide a detailed description for several examples:
* [Lid-driven cavity flow](doc/lidDrivenCavity.md);
* [Flow around an impulsively started cylinder](doc/cylinder.md);
* [Flow generated by a flapping wing](doc/flapping.md);

with comparison to experimental and computational results from other studies.

---

## List of examples

* `lidDrivenCavityRe100`: lid-driven cavity flow at Reynolds number 100. (~5 seconds on a K40 GPU.)
* `lidDrivenCavityRe1000`: lid-driven cavity flow at Reynolds number 1000. (~81 seconds on a K40 GPU.)
* `lidDrivenCavityRe3200`: lid-driven cavity flow at Reynolds number 3200. (~5 minutes on a K40 GPU.)
* `lidDrivenCavityRe5000`: lid-driven cavity flow at Reynolds number 5000. (~11 minutes on a K40 GPU.)
* `cylinderRe40`: flow over a circular cylinder at Reynolds number 40. (~4 minutes on a K40 GPU.)
* `cylinderRe100`: flow over a circular cylinder at Reynolds number 100. The
initial flow field has an asymmetric perturbation that triggers instability in
the flow and vortex shedding is observed in the wake.
* `cylinderRe150`: flow over a circular cylinder at Reynolds number 150. The
initial flow field has an asymmetric perturbation that triggers instability in
the flow and vortex shedding is observed in the wake.
* `cylinderRe200`: flow over a circular cylinder at Reynolds number 150. The
initial flow field has an asymmetric perturbation that triggers instability in
the flow and vortex shedding is observed in the wake.
* `cylinderRe550`: initial flow over an impulsively started cylinder at
Reynolds number 550. (~4 minutes on a K40 GPU.)
* `cylinderRe3000`: initial flow over an impulsively started cylinder at
Reynolds number 3000. (~25 minutes on a K40 GPU.)
* `flappingRe75`: flow around a flapping foil at Reynolds number 75. (~10 hours on a K20 GPU.)
* `heavingRe500`: flow around heaving foil at Reynolds number 500.
* `oscillatingCylindersRe100`: flow across two oscillating cylinders at Reynolds number 100. (~25 minutes on a K20 GPU.)
* `snakeRe1000AoA30`: flow around the cross-section of a gliding snake forming a 30-degree angle of attack with the freestream at Reynolds number 1000. (~65 hours on a K20 GPU.)
* `snakeRe1000AoA35`: flow around the cross-section of a gliding snake forming a 35-degree angle of attack with the freestream at Reynolds number 1000. (~65 hours on a K20 GPU.)
* `snakeRe1000AoA30`: flow around the cross-section of a gliding snake forming a 30-degree angle of attack with the freestream at Reynolds number 2000. (~65 hours on a K20 GPU.)
* `snakeRe1000AoA35`: flow around the cross-section of a gliding snake forming a 35-degree angle of attack with the freestream at Reynolds number 2000. (~65 hours on a K20 GPU.)

To run any of the examples listed above:

    > cd $CUIBM_DIR
    > make <examplename>

or change directory to the corresponding example's folder and

    > cuibm

The biggest cases (for the gliding snake) requires a GPU device with at least 4GB of memory.

To test an example on a smaller number of time steps, you can change the value of the parameters `nt` in the input file `simParams.yaml` of the example directory. (Note that the flow solution is saved every `nsave` time steps which value can be changed in the same input file.)

---

## Papers published using cuIBM

* Layton, S. K., Krishnan, A., Barba, L. A. (2011). "cuIBM--a GPU-accelerated immersed boundary method," unpublished. Preprint on arXiv: https://arxiv.org/abs/1109.3524v1
* Krishnan, A., Barba, L. A. (2013). "Validation of the cuIBM code for Navier-Stokes equations with immersed boundary methods," Figshare: https://doi.org/10.6084/m9.figshare.92789.v3
* Krishnan, A., Socha, J. J., Vlachos, P. V., Barba, L. A. (2014). "Lift and wakes of flying snakes," _Physics of Fluids_, __26__:031901, [10.1063/1.4866444](https://dx.doi.org/10.1063/1.4866444). Preprint on arXiv: https://arxiv.org/abs/1309.2969
* Mesnard, O., Barba, L. A. (2017). "Reproducible and replicable computational fluid dynamics: it's harder than you think," _IEEE/AIP Computing in Science and Engineering_, in press. Preprint on arXiv: https://arxiv.org/abs/1605.04339

---

## Contact

Please e-mail [Anush Krishnan](mailto:k.anush@gmail.com)
or [Olivier Mesnard](mailto:mesnardo@gwu.edu), if you have any
questions, suggestions or feedback.

To report bugs, please use the GitHub issue tracking system. We are also open to pull-requests.
