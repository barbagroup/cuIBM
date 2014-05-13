# this files contains all common options

# all suffixes that are subject to implicit rules (?) in the Makefiles
.SUFFIXES: .cu

# compiler
# nvcc - NVIDIA's CUDA compiler
CC = nvcc ${OSXOPTS}

# options for the compiler
# -arch=compute_20 compiles the program to use double precision on the GPU
# -O3 is an optimisation flag
OPTS = -arch=compute_20 -O3

# root directory of the project
# MAKEFILE_LIST is a special variable
# it lists all the makefiles in the current directory
# (Makefile and Makefile.inc in this case)
# lastword picks the last Makefile in that list
# dir returns the directory in which the last file is
# effectively giving us the current directory
PROJ_ROOT=$(dir $(lastword $(MAKEFILE_LIST)))

# include the programs in these folders while compiling
# src/include contains all the header files
# external/yaml-cpp/include contains all the headers for the YAML parser
INC = -I $(PROJ_ROOT)src/include -I $(PROJ_ROOT)external/yaml-cpp/include -I $(CUSP_DIR)

# external libraries used
LIBS = -L $(PROJ_ROOT)external/lib/ -lyaml-cpp

# used to generate the .a files (?)
AR = ar
ARFLAGS =

# (?)
RANLIB = ranlib

# command to delete files
RM = rm

# make command
MAKE = make

# directory which contains all the solver definitions
SOLVER_DIR = src/solvers/

# compile the .cu files to create the respective object files
# this is the suffix rule for converting .cu->.o
# suffix rules are defined as s1s2 when the rule is to convert from s1 to s2

# look up bash special parameters to find out what $? and $@ mean
# $^ lists all the dependencies (.cu files in this case)
# $? lists the dependencies that have changed since the last compile
# $@ lists the corresponding .o files (since the rule is from .cu to .o)
# otherwise, $@ is replaced by the name of the command
# -c creates an object file
# -o specifies the name of the file generated
.cu.o:
	$(CC) $(OPTS) $(INC) -c $? -o $@
