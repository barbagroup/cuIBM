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

//	array1d<int, memoryType>
//	     tags;
//	array1d<real, memoryType>
//	     coeffs;
	
	cusp::array1d<int, host_memory>
         tags;
	cusp::array1d<int, device_memory>
         tagsD;

	vecH coeffs;
	vecD coeffsD;

	
	void tagPoints();
//	void generateA();
//	void updateA();
	void generateL();
	void generateRN(int i);
	void updateRN();
	void generateBC1(int i);
	void updateBC1();
	void initialiseBodies();
	void updateBodies();
//	void updateSolverState(int i);
	
public:
	void initialise();
	std::string name()
	{
		return "Fadlun et al";
	}
};
