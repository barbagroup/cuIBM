#pragma once

#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <bodies.h>
/**
* Immersed boundary method described by Taira and Colonius (2007)
*/
template <typename memoryType>
class TairaColoniusSolver : public NavierStokesSolver<memoryType>
{
private:
	bodies<memoryType> B;
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