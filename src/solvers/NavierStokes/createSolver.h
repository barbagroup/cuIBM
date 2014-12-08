#pragma once

#include "NavierStokesSolver.h"
#include "DirectForcingSolver.h"
#include "FadlunEtAlSolver.h"
#include "TairaColoniusSolver.h"
#include "SLL0Solver.h"
#include "SLL1Solver.h"
#include "SLL2Solver.h"
#include "DiffusionSolver.h"
#include "DFModifiedSolver.h"
#include "DFImprovedSolver.h"

template <typename memoryType>
NavierStokesSolver<memoryType>* createSolver(parameterDB &paramDB, domain &domInfo);
