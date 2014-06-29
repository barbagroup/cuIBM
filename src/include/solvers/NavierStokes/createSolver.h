/***************************************************************************//**
* \file createSolver.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Declaration of the function to create the Navier-Stokes solver
*/

#pragma once

#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <solvers/NavierStokes/DirectForcingSolver.h>
#include <solvers/NavierStokes/FadlunEtAlSolver.h>
#include <solvers/NavierStokes/TairaColoniusSolver.h>
#include <solvers/NavierStokes/SLL0Solver.h>
#include <solvers/NavierStokes/SLL1Solver.h>
#include <solvers/NavierStokes/SLL2Solver.h>

/***************************************************************************//**
* \fn
* \brief Create a Navier-Stokes solver
* \param paramDB database with all simulation parameters
* \param domInfo information related to the computational grid
* \return pointer to an object of the class NavierStokesSolver
*/
template <typename memoryType>
NavierStokesSolver<memoryType>* createSolver(parameterDB &paramDB, domain &domInfo);
