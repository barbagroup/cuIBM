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
/*
	// For the 1D interpolations
	cusp::array1d<int, host_memory>    tags;
	cusp::array1d<int, device_memory>  tagsD;
	vecH coeffs;
	vecD coeffsD;
*/	
	// For the 2D interpolations
	cusp::array1d<int, host_memory>    tagsX, tagsY;
	cusp::array1d<int, device_memory>  tagsXD, tagsYD;
	vecH coeffsX, coeffsY;
	vecD coeffsXD, coeffsYD;
	
	void tagPoints();
	void tagPoints(real *bx, real *by);

//	void generateA();
//	void updateA();

	void generateL();
	void generateQT(int *QTRows, int *QTCols, real *QTVals){}
	void generateQT();
	void updateQFadlun();
	
	void generateRN();
	void updateRN();
	void generateBC1();
	void updateBC1();

	void initialiseBodies();
	void updateBodies();
//	void updateSolverState();
	
public:
	void initialise();
	std::string name()
	{
		return "Fadlun et al";
	}
};
