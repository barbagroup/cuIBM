#pragma once

#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <bodies.h>
/**
* Immersed boundary method described by Fadlun et al (2000)
*/
class FadlunEtAlSolver : public NavierStokesSolver
{
private:
	bodies B;
	void generateA();
	void updateA();
	void initialiseBodies();
	void updateBodies();
	void updateSolverState();
	
public:
	void initialise();
	//friend NavierStokesSolver* allocator::createSolver(flowDescription &flow_desc, simulationParameters &sim_par, domain &dom_info);
	//friend class NavierStokesSolver;
};