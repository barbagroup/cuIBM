#pragma once

#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <bodies.h>
/**
* Immersed boundary method described by Fadlun et al (2000)
*/
template <typename Matrix, typename Vector>
class FadlunEtAlSolver : public NavierStokesSolver<Matrix, Vector>
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
};