# file: Makefile
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Compiles and links cuIBM code.


# compiler: nvcc is NVIDIA's CUDA compiler
CC = nvcc $(OSXOPTS)

# compiler options
# -arch=compute_20: compile to use double precision on GPU
# -O3: optimization flag
CCFLAGS = -arch=sm_35 -O3

# variables
RM = rm
MAKE = make

# root directory of the project
# return the absolute path of the directory  where is located the Makefile
# variable MAKEFILE_LIST lists all Makefiles in the working directory
# `lastword` picks the last element of the list
PROJ_ROOT = $(abspath $(dir $(lastword $(MAKEFILE_LIST))))

# source code directory
SRC_DIR = $(PROJ_ROOT)/src

# directory where object files are stored
BUILD_DIR = $(PROJ_ROOT)/build

# directory where binary executables are stored
BIN_DIR = $(PROJ_ROOT)/bin

# extension of source files
SRC_EXT = .cu

# cuIBM executable
TARGET = bin/cuIBM

# list all source files in the source directory
SRCS = $(shell find $(SRC_DIR) -type f -name *$(SRC_EXT))
# absolute path of all object files to be created
OBJS = $(patsubst $(SRC_DIR)/%, $(BUILD_DIR)/%, $(SRCS:$(SRC_EXT)=.o))

# include header files from cuIBM and CUSP library
INC = -I $(SRC_DIR) -I $(CUSP_DIR)

# path of the YAML static library
EXT_LIBS = $(PROJ_ROOT)/external/lib/libyaml-cpp.a
# include YAML header files
INC += -I $(PROJ_ROOT)/external/yaml-cpp/include


.PHONY: all

all: $(TARGET)

$(TARGET): $(OBJS) $(EXT_LIBS)
	@echo "\nLinking ..."
	@mkdir -p $(BIN_DIR)
	$(CC) $^ -o $@

$(EXT_LIBS):
	@echo "\nCreating static library $@ ..."
	cd external; $(MAKE) $(MFLAGS) all

$(BUILD_DIR)/%.o: $(SRC_DIR)/%$(SRC_EXT)
	@mkdir -p $(@D)
	$(CC) $(CCFLAGS) $(INC) -c $< -o $@

################################################################################

.PHONY: unittests unittest

UNITTESTS_DIR = $(PROJ_ROOT)/unitTests
TEST_SRCS = $(wildcard $(TEST_DIR)/*$(SRC_EXT))
TEST_OBJS = $(TEST_SRCS:$(SRC_EXT)=.o) $(filter-out $(BUILD_DIR)/cuIBM.o, $(OBJS))
TEST_BIN = $(BIN_DIR)/unitTests/$(lastword $(subst /, , $(TEST_DIR)))

unittests: convectionTerm_test diffusionTerm_test

convectionTerm_test:  export TEST_DIR = $(UNITTESTS_DIR)/convectionTerm
convectionTerm_test: 
	$(MAKE) unittest

diffusionTerm_test: export TEST_DIR = $(UNITTESTS_DIR)/diffusionTerm
diffusionTerm_test:
	$(MAKE) unittest

unittest: $(TEST_OBJS) $(EXT_LIBS)
	@echo "\nUnit-test: $(TEST_DIR) ..."
	@mkdir -p $(BIN_DIR)/unitTests
	$(CC) $^ -o $(TEST_BIN)

%.o: %$(SRC_EXT)
	$(CC) $(CCFLAGS) $(INC) -I $(TEST_DIR) -c $< -o $@

################################################################################

.PHONY: doc

DOC_DIR = $(PROJ_ROOT)/doc
DOXYGEN = doxygen

doc:
	@echo "\nGenerating Doxygen documentation ..."
	cd $(DOC_DIR); $(DOXYGEN) Doxyfile

################################################################################

.PHONY: clean cleanexternal cleanunittests cleandoc cleanall

clean:
	@echo "\nCleaning cuIBM ..."
	$(RM) -rf $(BUILD_DIR) $(BIN_DIR)

cleanexternal:
	@echo "\nCleaning external YAML ..."
	cd external; $(MAKE) $(MFLAGS) clean

cleanunittests:
	@echo "\nCleaning unitTests ..."
	find $(UNITTESTS_DIR) -type f -name *.o -delete
	$(RM) -rf $(BIN_DIR)/unitTests

cleandoc:
	@echo "\nCleaning documentation ..."
	find $(DOC_DIR) ! -name 'Doxyfile' -type f -delete
	find $(DOC_DIR)/* ! -name 'Doxyfile' -type d -delete

cleanall: clean cleanexternal cleanunittests cleandoc

################################################################################

# commands to run unit-tests

testConvection:
	bin/unitTests/convectionTerm -caseFolder cases/unitTests/convectionTerm/6
	bin/unitTests/convectionTerm -caseFolder cases/unitTests/convectionTerm/12
	bin/unitTests/convectionTerm -caseFolder cases/unitTests/convectionTerm/24

testDiffusion:
	bin/unitTests/diffusionTerm -caseFolder cases/unitTests/convectionTerm/6
	bin/unitTests/diffusionTerm -caseFolder cases/unitTests/convectionTerm/12
	bin/unitTests/diffusionTerm -caseFolder cases/unitTests/convectionTerm/24
	bin/unitTests/diffusionTerm -caseFolder cases/unitTests/convectionTerm/48

################################################################################

# commands to run cuIBM simulations

lidDrivenCavityRe100:
	bin/cuIBM -caseFolder cases/lidDrivenCavity/Re100

lidDrivenCavityRe1000:
	bin/cuIBM -caseFolder cases/lidDrivenCavity/Re1000

cylinder:
	bin/cuIBM -caseFolder cases/cylinder/test

cylinderRe40:
	bin/cuIBM -caseFolder cases/cylinder/Re40

cylinderRe550:
	bin/cuIBM -caseFolder cases/cylinder/Re550

cylinderRe3000:
	bin/cuIBM -caseFolder cases/cylinder/Re3000

cylinderRe75:
	bin/cuIBM -caseFolder cases/cylinder/Re75

cylinderRe100:
	bin/cuIBM -caseFolder cases/cylinder/Re100

cylinderRe150:
	bin/cuIBM -caseFolder cases/cylinder/Re150

cylinderDirectForcing:
	bin/cuIBM -caseFolder cases/cylinder/Re40 -ibmScheme DirectForcing

snakeRe1000AOA30:
	time bin/cuIBM -caseFolder cases/flyingSnake/Re1000_AoA30

snakeRe1000AOA35:
	time bin/cuIBM -caseFolder cases/flyingSnake/Re1000_AoA35

snakeRe2000AOA30:
	time bin/cuIBM -caseFolder cases/flyingSnake/Re2000_AoA30

snakeRe2000AOA35:
	time bin/cuIBM -caseFolder cases/flyingSnake/Re2000_AoA35

flappingRe75:
	time bin/cuIBM -caseFolder cases/flappingRe75

oscillatingCylinders:
	time bin/cuIBM -caseFolder cases/oscillatingCylinders
