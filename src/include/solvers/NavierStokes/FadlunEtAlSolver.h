/**
* @file  FadlunEtAlSolver.h
* @brief Modified direct forcing immersed boundary method.
*
* <b>Combined immersed-boundary finite-difference methods for three-dimensional complex flow simulations</b> \n
* Fadlun, E A, Verzicco, R, Orlandi, P, and Mohd-Yusof, J \n
* Journal of Computational Physics \n
* Volume 161 Number 1 \n
* 2000
*/

#pragma once

#include <solvers/NavierStokes/NSWithBody.h>
/**
* Immersed boundary method described by Fadlun et al (2000)
*/
template <typename memoryType>
class FadlunEtAlSolver : public NSWithBody<memoryType>
{
private:
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
	vecH coeffsX,  coeffsY,  uvX,  uvY;
	vecD coeffsXD, coeffsYD, uvXD, uvYD;
	
	void tagPoints();
	//void tagPoints(real *bx, real *by);
	void tagPoints(real *bx, real *by, real *uB, real *vB);

	void generateA(real alpha);
//	void updateA();

	void generateL();
	void generateQT(int *QTRows, int *QTCols, real *QTVals){} // check this!
	void generateQT();
	void updateQ();
	
	void assembleRHS1();
	void updateRHS1();

	void updateSolverState();
	void  calculateForceF();
	
public:
	FadlunEtAlSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	void initialise();
	std::string name()
	{
		return "Fadlun et al";
	}
};
