/***************************************************************************//**
 * \file FadlunEtAlSolver.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c FadlunEtAlSolver.
 */


#include "FadlunEtAlSolver.h"
#include <sys/stat.h>


/**
 * \brief Generates the matrix \c QT for FadlunEtAlSolver.
 */
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::generateQT()
{
	NavierStokesSolver<memoryType>::generateQT();
	cusp::transpose(NavierStokesSolver<memoryType>::Q, NavierStokesSolver<memoryType>::QT);
}

/**
 * \brief Constructor. Copies the database and information about the computational grid.
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
