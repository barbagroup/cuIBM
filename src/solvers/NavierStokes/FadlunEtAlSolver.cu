/***************************************************************************//**
* \file FadlunEtAlSolver.cu
* \author Krishnan, A. (anush@bu.edu)
* \brief Definition of the methods of the class \c FadlunEtAlSolver
*/

#include <solvers/NavierStokes/FadlunEtAlSolver.h>
#include <sys/stat.h>

/**
* \brief To be documented
*/
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateQT()
{
	NavierStokesSolver<memoryType>::generateQT();
	cusp::transpose(NavierStokesSolver<memoryType>::Q, NavierStokesSolver<memoryType>::QT);
}

/**
* \brief To be documented
*/
template <typename memoryType>
FadlunEtAlSolver<memoryType>::FadlunEtAlSolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

// specialization of the class
template class FadlunEtAlSolver<host_memory>;
template class FadlunEtAlSolver<device_memory>;
