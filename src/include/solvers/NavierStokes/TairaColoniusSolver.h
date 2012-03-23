#pragma once

#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <bodies.h>
/**
* Immersed boundary method described by Fadlun et al (2000)
*/
template <typename Matrix, typename Vector>
class TairaColoniusSolver : public NavierStokesSolver<Matrix, Vector>
{
private:
	bodies B;
	void initialiseArrays();
	void generateQT();
	void updateQT();
	void generateBC2();
	void initialiseBodies();
	void updateBodies();
	void updateSolverState();
	
public:
	void initialise();
	std::string name()
	{
		return "Taira & Colonius";
	}
};