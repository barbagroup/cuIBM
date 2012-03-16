.SUFFIXES: .cu .o

INCLUDE = src/include
#PRECISION = -DPRECISION=double
#OPTS = -fopenmp -Wall -pedantic -O3 -DNDEBUG
#LIBS =  -lgomp -lboost_program_options -lboost_filesystem -lboost_system
INC = -I ${INCLUDE}
#LIBDIR = -L ${BOOST_LIB}
#CUDA_INC = -I ${CUSP_DIR}
#CUDA_OPTS = -arch=sm_20 -O3 -use_fast_math -Xcompiler -fopenmp  #-g

NSSOLVERS = src/solvers/NavierStokes
NS = src/solvers/NavierStokes/NavierStokesSolver

OBJ += $(NSSOLVERS)/NavierStokesSolver.o\
       $(NSSOLVERS)/FadlunEtAlSolver.o

main: cuibm

cuibm: $(OBJ)
	nvcc $(INC) -c src/cuIBM.cu -o src/cuIBM.o
	nvcc $? src/cuIBM.o -o bin/cuIBM

clean:
	rm -rf $(NS)/*.o $(NSSOLVERS)/*.o src/*.o bin/cuIBM

.cu.o:
	nvcc $(INC) -c $< -o $@

.PHONY: cuibm
