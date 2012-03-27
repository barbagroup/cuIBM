#pragma once

#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <bodies.h>
/**
* Immersed boundary method described by Fadlun et al (2000)
*/
template <typename memoryType>
class FadlunEtAlSolver : public NavierStokesSolver<memoryType>
{
private:
	bodies<memoryType> B;
	void generateA();
	void updateA();
	void initialiseBodies();
	void updateBodies();
	void updateSolverState();
	
public:
	void initialise();
};