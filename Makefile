# file: Makefile
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Compiles and links cuIBM code.


# compiler
CC = nvcc $(OSXOPTS)
CCFLAGS = -arch=compute_20 -O3

# variables
RM = rm
MAKE = make

# directories
PROJECT_ROOT = $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
SRC_DIR = $(PROJECT_ROOT)/src
BUILD_DIR = $(PROJECT_ROOT)/build
BIN_DIR = $(PROJECT_ROOT)/bin

SRC_EXT = .cu

# executable
TARGET = bin/cuIBM

SRCS = $(shell find $(SRC_DIR) -type f -name *$(SRC_EXT))
OBJS = $(patsubst $(SRC_DIR)/%, $(BUILD_DIR)/%, $(SRCS:$(SRC_EXT)=.o))

# header files from cuIBM and cusp
INC = -I $(SRC_DIR) -I $(CUSP_DIR)

# yaml
EXT_LIBS = $(PROJECT_ROOT)/external/lib/libyaml-cpp.a
INC += -I $(PROJECT_ROOT)/external/yaml-cpp/include


.PHONY: all

all: $(TARGET)

$(TARGET): $(OBJS) $(EXT_LIBS)
	@echo "\nLinking ..."
	@mkdir -p $(BIN_DIR)
	$(CC) $^ -o $@

external/lib/libyaml-cpp.a:
	@echo "\nCreating static library $@ ..."
	cd external; $(MAKE) $(MFLAGS) all

$(BUILD_DIR)/%.o: $(SRC_DIR)/%$(SRC_EXT)
	@mkdir -p $(@D)
	$(CC) $(CCFLAGS) $(INC) -c $< -o $@

################################################################################

.PHONY: unittests unittest

UNITTESTS_DIR = $(PROJECT_ROOT)/unitTests
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

.PHONY: clean cleanunittests cleanall

clean:
	@echo "\nCleaning ..."
	$(RM) -rf $(LIB_DIR) $(BUILD_DIR) $(BIN_DIR)

cleanunittests:
	@echo "\nCleaning unitTests ..."
	find $(UNITTESTS_DIR) -type f -name *.o -delete
	$(RM) -rf $(BIN_DIR)/unitTests

cleanall: clean cleanunittests
	cd external; $(MAKE) $(MFLAGS) clean

################################################################################

# unit-tests
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

# cases
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
