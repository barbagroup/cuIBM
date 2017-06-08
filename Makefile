# file: Makefile
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Compiles and links cuIBM code.


# compiler: nvcc is NVIDIA's CUDA compiler
CC = nvcc $(OSXOPTS)

# compiler options
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
TARGET = bin/cuibm

# list all source files in the source directory
SRCS = $(shell find $(SRC_DIR) -type f -name *$(SRC_EXT))
# absolute path of all object files to be created
OBJS = $(patsubst $(SRC_DIR)/%, $(BUILD_DIR)/%, $(SRCS:$(SRC_EXT)=.o))

# include header files from cuIBM, CUSP library, and Boost library
INC = -I$(SRC_DIR) -I$(CUSP_DIR) -I$(BOOST_DIR)

# path of the YAML static library
EXT_LIBS = $(PROJ_ROOT)/external/lib/libyaml-cpp.a
# include YAML header files
INC += -I $(PROJ_ROOT)/external/yaml-cpp-0.5.3/include


.PHONY: all

all: $(TARGET)

$(TARGET): $(EXT_LIBS) $(OBJS)
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

.PHONY: tests testConvectionTerm testDiffusionTerm unittest

TESTS_DIR = $(PROJ_ROOT)/tests
TEST_SRCS = $(wildcard $(TEST_SRC_DIR)/*$(SRC_EXT))
TEST_OBJS = $(patsubst $(TEST_SRC_DIR)/%, $(TEST_BUILD_DIR)/%, $(TEST_SRCS:$(SRC_EXT)=.o))
TEST_OBJS += $(filter-out $(BUILD_DIR)/cuIBM.o, $(OBJS))
TEST_BIN = $(BIN_DIR)/tests/$(lastword $(subst /, , $(TEST_SRC_DIR)))

tests: testConvectionTerm testDiffusionTerm

testConvectionTerm: export TEST_SRC_DIR=$(TESTS_DIR)/convectionTerm
testConvectionTerm: export TEST_BUILD_DIR=$(BUILD_DIR)/tests/convectionTerm
testConvectionTerm:
	$(MAKE) unittest
	@echo "" > data.txt
	$(TEST_BIN) -directory tests/cases/6
	$(TEST_BIN) -directory tests/cases/12
	$(TEST_BIN) -directory tests/cases/24
	$(TEST_BIN) -directory tests/cases/48
	@python $(TESTS_DIR)/convergence.py --filepath data.txt
	@$(RM) data.txt

testDiffusionTerm: export TEST_SRC_DIR=$(TESTS_DIR)/diffusionTerm
testDiffusionTerm: export TEST_BUILD_DIR=$(BUILD_DIR)/tests/diffusionTerm
testDiffusionTerm:
	$(MAKE) unittest
	@echo "" > data.txt
	$(TEST_BIN) -directory tests/cases/6
	$(TEST_BIN) -directory tests/cases/12
	$(TEST_BIN) -directory tests/cases/24
	$(TEST_BIN) -directory tests/cases/48
	@python $(TESTS_DIR)/convergence.py --filepath data.txt
	@$(RM) data.txt

unittest: $(TEST_OBJS) $(EXT_LIBS)
	@echo "\nTest: $(TEST_DIR) ..."
	@mkdir -p $(BIN_DIR)/tests
	$(CC) $^ -o $(TEST_BIN)

$(TEST_BUILD_DIR)/%.o: $(TEST_SRC_DIR)/%$(SRC_EXT)
	@mkdir -p $(@D)
	$(CC) $(CCFLAGS) $(INC) -I $(TEST_SRC_DIR) -c $< -o $@

################################################################################

.PHONY: doc

DOC_DIR = $(PROJ_ROOT)/doc
DOXYGEN = doxygen

doc:
	@echo "\nGenerating Doxygen documentation ..."
	cd $(DOC_DIR); $(DOXYGEN) Doxyfile

################################################################################

.PHONY: clean cleanexternal cleantests cleanall

clean:
	@echo "\nCleaning cuIBM ..."
	$(RM) -rf $(BUILD_DIR) $(BIN_DIR)

cleanexternal:
	@echo "\nCleaning external YAML ..."
	cd external; $(MAKE) $(MFLAGS) clean

cleantests:
	@echo "\nCleaning tests ..."
	$(RM) -rf $(BIN_DIR)/tests
	$(RM) -rf $(BUILD_DIR)/tests

cleanall: clean cleanexternal cleantests

################################################################################

# commands to run cuIBM simulations

lidDrivenCavityRe100:
	bin/cuibm -directory examples/lidDrivenCavity/Re100

lidDrivenCavityRe1000:
	bin/cuibm -directory examples/lidDrivenCavity/Re1000

lidDrivenCavityRe3200:
	bin/cuibm -directory examples/lidDrivenCavity/Re3200

lidDrivenCavityRe5000:
	bin/cuibm -directory examples/lidDrivenCavity/Re5000

cylinderRe40:
	bin/cuibm -directory examples/cylinder/Re40

cylinderRe550:
	bin/cuibm -directory examples/cylinder/Re550

cylinderRe3000:
	bin/cuibm -directory examples/cylinder/Re3000

cylinderRe100:
	bin/cuibm -directory examples/cylinder/Re100

cylinderRe150:
	bin/cuibm -directory examples/cylinder/Re150

cylinderRe200:
	bin/cuibm -directory examples/cylinder/Re200

cylinderDirectForcing:
	bin/cuibm -directory examples/cylinder/Re40 -ibmScheme DirectForcing

snakeRe1000AOA30:
	bin/cuibm -directory examples/flyingSnake/Re1000AoA30

snakeRe1000AOA35:
	bin/cuibm -directory examples/flyingSnake/Re1000AoA35

snakeRe2000AOA30:
	bin/cuibm -directory examples/flyingSnake/Re2000AoA30

snakeRe2000AOA35:
	bin/cuibm -directory examples/flyingSnake/Re2000AoA35

flappingRe75:
	bin/cuibm -directory examples/flapping/Re75

heavingRe500:
	bin/cuibm -directory examples/heaving/Re500

oscillatingCylindersRe100:
	bin/cuibm -directory examples/oscillatingCylinders/Re100

convergence:
	cd examples/convergence/lidDrivenCavityRe100; $(MAKE) all
