#include "DiffusionSolver.h"
#include <sys/stat.h>

template <typename memoryType>
void DiffusionSolver<memoryType>::initialise()
{
	DirectForcingSolver<memoryType>::initialise();

	for(int i=0; i<NavierStokesSolver<memoryType>::intgSchm.subSteps; i++)
	{
		NavierStokesSolver<memoryType>::intgSchm.gamma[i] = 0.0;
		NavierStokesSolver<memoryType>::intgSchm.zeta[i] = 0.0;
	}
}

template <typename memoryType>
void DiffusionSolver<memoryType>::solvePoisson()
{
}

template <typename memoryType>
void DiffusionSolver<memoryType>::projectionStep()
{
	NavierStokesSolver<memoryType>::logger.startTimer("projectionStep");

	NavierStokesSolver<memoryType>::q = NavierStokesSolver<memoryType>::qStar;

	NavierStokesSolver<memoryType>::logger.stopTimer("projectionStep");
}

template <typename memoryType>
DiffusionSolver<memoryType>::DiffusionSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

template class DiffusionSolver<host_memory>;
template class DiffusionSolver<device_memory>;
