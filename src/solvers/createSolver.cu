/***************************************************************************//**
 * \file createSolver.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the function that creates the Navier-Stokes solver.
 */


#include "createSolver.h"


/**
 * \brief Creates the appropriate Navier-Stokes solver.
 *
 * The type of solver depends on the IBM used.
 * If there is no immersed body in the flow, 
 * an instance of the class \c NavierStokesSolver is created.
 * If \c TAIRA_COLONIUS is the IBM used, 
 * an instance of the class \c TairaColoniusSolver is created, 
 * that inherits from the class \c NavierStokesSolver.
 *
 * \param paramDB database that contains all the simulation parameters
 * \param domInfo information about the computational grid
 *
 * \return an instance of the appropriate class
 */
template <typename memoryType>
NavierStokesSolver<memoryType>* createSolver(parameterDB &paramDB, domain &domInfo)
{
	ibmScheme ibm = paramDB["simulation"]["ibmScheme"].get<ibmScheme>();
	NavierStokesSolver<memoryType> *solver = 0;
	switch(ibm)
	{
		case SAIKI_BIRINGEN:
			break;
		case TAIRA_COLONIUS:
			solver = new TairaColoniusSolver<memoryType>(&paramDB, &domInfo);
			break;
		case NAVIER_STOKES:
			solver = new NavierStokesSolver<memoryType>(&paramDB, &domInfo);
			break;
		case DIRECT_FORCING:
			solver = new DirectForcingSolver<memoryType>(&paramDB, &domInfo);
			break;
		case FADLUN_ET_AL:
			solver = new FadlunEtAlSolver<memoryType>(&paramDB, &domInfo);
			break;
		case DIFFUSION:
			solver = new DiffusionSolver<memoryType>(&paramDB, &domInfo);
			break;
		case DF_MODIFIED:
			solver = new DFModifiedSolver<memoryType>(&paramDB, &domInfo);
			break;
		case FEA_MODIFIED:
			solver = new FEAModifiedSolver<memoryType>(&paramDB, &domInfo);
			break;
		case DF_IMPROVED:
			solver = new DFImprovedSolver<memoryType>(&paramDB, &domInfo);
			break;
	}
	std::cout << "\nImmersed Boundary Method" << std::endl;
	std::cout <<   "------------------------" << std::endl;
	std::cout << solver->name() << '\n' << std::endl;
	return solver;
}

// specializations of the template function
template NavierStokesSolver<host_memory>* createSolver(parameterDB &paramDB, domain &domInfo);
template NavierStokesSolver<device_memory>* createSolver(parameterDB &paramDB, domain &domInfo);
