# file: Makefile
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Compiles cuIBM code.


CC = nvcc $(OSXOPTS)
CCFLAGS = -arch=compute_20 -O3

RM = rm
MAKE = make

PROJECT_ROOT = $(dir $(lastword $(MAKEFILE_LIST)))
SRC_DIR = $(PROJECT_ROOT)src
BUILD_DIR = $(PROJECT_ROOT)build
BIN_DIR = $(PROJECT_ROOT)bin

UNITTESTS_DIR = $(SRC_DIR)/unitTests

SRC_EXT = .cu

TARGET = bin/cuIBM
SRCS = $(shell find $(SRC_DIR) -type f -name *$(SRC_EXT) ! -path $(UNITTESTS_DIR)"*")
OBJS = $(patsubst $(SRC_DIR)/%, $(BUILD_DIR)/%, $(SRCS:$(SRC_EXT)=.o))
INC = -I $(SRC_DIR) -I $(CUSP_DIR)

EXT_LIBS = $(PROJECT_ROOT)external/lib/libyaml-cpp.a
INC += -I $(PROJECT_ROOT)external/yaml-cpp/include


.PHONY: all $(TARGET) clean cleanall

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

clean:
	@echo "\nCleaning ..."
	$(RM) -rf $(LIB_DIR) $(BUILD_DIR) $(BIN_DIR)

cleanall: clean
	cd external; $(MAKE) $(MFLAGS) clean


unitTests/convectionTerm: $(TESTLIBS) $(EXT_LIBS) src/helpers.o src/parameterDB.o src/unitTests/convectionTerm.o
	${CC} $^ -o bin/unitTests/convectionTerm

unitTests/diffusionTerm: $(TESTLIBS) $(EXT_LIBS) src/helpers.o src/parameterDB.o src/unitTests/diffusionTerm.o
	${CC} $^ -o bin/unitTests/diffusionTerm

testConvection:
	bin/unitTests/convectionTerm -caseFolder cases/unitTests/convectionTerm/6
	bin/unitTests/convectionTerm -caseFolder cases/unitTests/convectionTerm/12
	bin/unitTests/convectionTerm -caseFolder cases/unitTests/convectionTerm/24

testDiffusion:
	bin/unitTests/diffusionTerm -caseFolder cases/unitTests/convectionTerm/6
	bin/unitTests/diffusionTerm -caseFolder cases/unitTests/convectionTerm/12
	bin/unitTests/diffusionTerm -caseFolder cases/unitTests/convectionTerm/24
	bin/unitTests/diffusionTerm -caseFolder cases/unitTests/convectionTerm/48

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
