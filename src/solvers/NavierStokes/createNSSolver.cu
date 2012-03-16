#include <solvers/NavierStokes/createNSSolver.h>

namespace allocator
{
	/**
	* \param a Description
	* \return
	*/
	NavierStokesSolver* createSolver(flowDescription &flow_desc, simulationParameters &sim_par, domain &dom_info)
	{
		NavierStokesSolver *solver;
		switch(sim_par.ibmSch)
		{
			case NAVIER_STOKES:
				solver = new NavierStokesSolver;
				break;
			case FADLUN_ET_AL:
				solver = new FadlunEtAlSolver;
				break;
		}
		solver->flowDesc = &flow_desc;
		solver->simPar = &sim_par;
		solver->domainInfo = &dom_info;
		return solver;
	}
}