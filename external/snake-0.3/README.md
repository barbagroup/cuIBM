# Snake

A collection of Python modules to post-process the numerical solution
of simulations run with one of the following software:
* [PetIBM](https://github.com/barbagroup/PetIBM)
* [cuIBM](https://github.com/barbaGroup/cuIBM)
* [IBAMR](https://github.com/IBAMR/IBAMR)
* [OpenFOAM](www.openfoam.com)


## Installation

### Dependencies (last tested)

* Python-2.7 or Python-3.5
* Numpy-1.11.1
* Scipy-0.17.1
* Matplotlib-1.5.3
* Pandas-0.19.2

To install the specific version of these packages we recommend using either [conda](https://conda.io/docs/get-started.html).

To create a new environment for using `snake` with Python-3.5 with conda you can do the following:

    > conda create -n snakepy35 python=3.5 numpy scipy matplotlib pandas
    > source activate snakepy35

or with Python-2.7:

    > conda create -n snakepy27 python=2.7 numpy scipy matplotlib pandas
    > source activate snakepy27

#### Optional
* PETSc-3.7.4 (for PetIBM post-processing; only the Python scripts are needed)
* VisIt-2.12.1 (IBAMR post-processing)
* OpenFOAM-2.3.1 and ThirdParty-2.3.1 (for OpenFOAM post-processing)


### Install `snake` using `git`

    > git clone https://github.com/mesnardo/snake.git
    > cd snake
    > python setup.py install
    > export SNAKE=$PWD

### Install `snake` using a specific release

    > VERSION="X.X.X"
    > wget https://github.com/mesnardo/snake/archive/v${VERSION}.tar.gz
    > mkdir snake-${VERSION}
    > tar -xzf v${VERSION}.tar.gz -C snake-${VERSION} --strip-components=1
    > rm -f v${VERSION}.tar.gz
    > cd snake-${VERSION}
    > python setup.py install
    > export SNAKE=$PWD

where `X.X.X` defines the version number to use.


## Notes

Some of the modules call the environment variable `$SNAKE` defined as the local directory of the `snake` repository.
The variable can be set by adding the following line to your `.bashrc` 
or `.bash_profile` file:

    > export SNAKE="/path/to/snake/directory"

The module for PetIBM calls the environment variable `$PETSC_DIR` defined as the local directory of PETSc.
The variable can be set by adding the following line to your `.bashrc` 
or `.bash_profile` file:

    > export PETSC_DIR="/path/to/petsc/directory"


Some examples to post-process the numerical solution from the four codes are provided in the `examples` folder.


Suggestions and bug-fix are welcome.
Contact: `mesnardo@gwu.edu`
