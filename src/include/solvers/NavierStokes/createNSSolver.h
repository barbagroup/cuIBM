#pragma once

#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <solvers/NavierStokes/FadlunEtAlSolver.h>

namespace allocator
{
	/**
	* \param a Description
	* \return
	*/
	NavierStokesSolver* createSolver(flowDescription &flow_desc, simulationParameters &sim_par, domain &dom_info);
}