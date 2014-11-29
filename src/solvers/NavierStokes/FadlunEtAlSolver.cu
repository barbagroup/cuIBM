#include "FadlunEtAlSolver.h"
#include <sys/stat.h>

template <typename memoryType>
FadlunEtAlSolver<memoryType>::FadlunEtAlSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::initialise()
{
	DirectForcingSolver<memoryType>::initialise();

	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
        ny = NavierStokesSolver<memoryType>::domInfo->ny;

	int numP  = nx*ny;

	pressure.resize(numP);
	cusp::blas::fill(pressure, 0.0);
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateQT()
{
	NavierStokesSolver<memoryType>::generateQT();
	cusp::transpose(NavierStokesSolver<memoryType>::Q, NavierStokesSolver<memoryType>::QT);
}

// rn = rn - zeta * Q.pressure
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::calculateExplicitLambdaTerms()
{
	// temp1 = Q.pressure
	cusp::multiply(NavierStokesSolver<memoryType>::Q, pressure, NavierStokesSolver<memoryType>::temp1);
	// rn = rn - temp1
	cusp::blas::axpy(NavierStokesSolver<memoryType>::temp1, NavierStokesSolver<memoryType>::rn, -1.0);
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::projectionStep()
{
	NavierStokesSolver<memoryType>::projectionStep();

	NavierStokesSolver<memoryType>::logger.startTimer("projectionStep");
	cusp::blas::axpy(NavierStokesSolver<memoryType>::lambda, pressure , 1.0);
	NavierStokesSolver<memoryType>::logger.stopTimer("projectionStep");
}

template class FadlunEtAlSolver<host_memory>;
template class FadlunEtAlSolver<device_memory>;
