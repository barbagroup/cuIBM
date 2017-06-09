/**
 * \file DiffusionSolver.cu
 * \brief Implementation of the methods of the class \c DiffusionSolver.
 */


#include <sys/stat.h>

#include "DiffusionSolver.h"


/**
 * \brief Constructor. Copies the simulation parameters and the domain info.
 */
template <typename memoryType>
DiffusionSolver<memoryType>::DiffusionSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
} // DiffusionSolver


/**
 * \brief Sets the coefficients of the convection term to zero. This is done 
 *        because we are solving a pure diffusion problem.
 */
template <typename memoryType>
void DiffusionSolver<memoryType>::initialise()
{
	DirectForcingSolver<memoryType>::initialise();

	for(int i=0; i<NavierStokesSolver<memoryType>::intgSchm.subSteps; i++)
	{
		NavierStokesSolver<memoryType>::intgSchm.gamma[i] = 0.0;
		NavierStokesSolver<memoryType>::intgSchm.zeta[i] = 0.0;
	}
} // initialise


/**
 * \brief Does nothing.
 *
 * No Poisson equation needs to be solved for the unsteady diffusion problem.
 */
template <typename memoryType>
void DiffusionSolver<memoryType>::solvePoisson()
{
} // solvePoisson


/**
 * \brief The velocity at the next time step is the same as the intermediate velocity.
 */
template <typename memoryType>
void DiffusionSolver<memoryType>::projectionStep()
{
	NavierStokesSolver<memoryType>::logger.startTimer("projectionStep");

	NavierStokesSolver<memoryType>::q = NavierStokesSolver<memoryType>::qStar;

	NavierStokesSolver<memoryType>::logger.stopTimer("projectionStep");
} // projectionStep


template class DiffusionSolver<host_memory>;
template class DiffusionSolver<device_memory>;
