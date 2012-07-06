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
	coo_matrix<int, real, memoryType> E, ET;
	
	array1d<real, memoryType> FxX, FxY, FxU;

	void generateQT(int *QTRows, int *QTCols, real *QTVals){}
	void generateQT();
	void updateQT();

	void generateBC2();

	void initialiseBodies();
	void updateBodies();
	void updateSolverState();

	void calculateForce1();
	void calculateForceTC();
	void calculateForce();
	void generateE();
	
public:
	void initialise();
	std::string name()
	{
		return "Taira & Colonius";
	}
};
