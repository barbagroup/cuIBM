#pragma once

#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <solvers/NavierStokes/DirectForcingSolver.h>
#include <solvers/NavierStokes/TairaColoniusSolver.h>
#include <solvers/NavierStokes/SLL0Solver.h>
#include <solvers/NavierStokes/SLL1Solver.h>
#include <solvers/NavierStokes/SLL2Solver.h>

template <typename memoryType>
NavierStokesSolver<memoryType>* createSolver(parameterDB &paramDB, domain &domInfo);
