.SUFFIXES: .cu

# common definitions
CC = nvcc
OPTS = -arch=sm_20 -O3 #-Xcompiler -Wall
PROJ_ROOT=$(dir $(lastword $(MAKEFILE_LIST)))
INC = -I $(PROJ_ROOT)src/include -I $(PROJ_ROOT)external/yaml-cpp/include
LIBS = -L $(PROJ_ROOT)external/lib/ -lyaml-cpp
AR = ar
ARFLAGS =
RANLIB = ranlib
RM = rm
MAKE = make

SOLVER_DIR = src/solvers/

.cu.o:
	$(CC) $(OPTS) $(INC) -c $? -o $@
