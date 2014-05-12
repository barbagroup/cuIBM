#include <solvers/NavierStokes/createSolver.h>

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
		case SLL0:
			solver = new SLL0Solver<memoryType>(&paramDB, &domInfo);
			break;
		case SLL1:
			solver = new SLL1Solver<memoryType>(&paramDB, &domInfo);
			break;
		case SLL2:
			solver = new SLL2Solver<memoryType>(&paramDB, &domInfo);
			break;
	}
	std::cout << "\nImmersed Boundary Method" << std::endl;
	std::cout <<   "------------------------" << std::endl;
	std::cout << solver->name() << '\n' << std::endl;
	return solver;
}

template NavierStokesSolver<host_memory>* createSolver(parameterDB &paramDB, domain &domInfo);
template NavierStokesSolver<device_memory>* createSolver(parameterDB &paramDB, domain &domInfo);
