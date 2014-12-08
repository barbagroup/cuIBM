#include "FadlunEtAlSolver.h"
#include <sys/stat.h>

template <typename memoryType>
FadlunEtAlSolver<memoryType>::FadlunEtAlSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateQT()
{
	NavierStokesSolver<memoryType>::generateQT();
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateC()
{
	NavierStokesSolver<memoryType>::generateC();
}

template class FadlunEtAlSolver<host_memory>;
template class FadlunEtAlSolver<device_memory>;

//----------------------------------------------------------------------------//

template <typename memoryType>
FEAModifiedSolver<memoryType>::FEAModifiedSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

template <typename memoryType>
void FEAModifiedSolver<memoryType>::generateQT()
{
	NavierStokesSolver<memoryType>::generateQT();
}

template <typename memoryType>
void FEAModifiedSolver<memoryType>::generateC()
{
	NavierStokesSolver<memoryType>::generateC();
}

template class FEAModifiedSolver<host_memory>;
template class FEAModifiedSolver<device_memory>;