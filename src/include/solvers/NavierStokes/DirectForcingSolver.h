/**
* @file  DirectForcingSolver.h
* @brief Modified direct forcing immersed boundary method.
*/

#pragma once

#include <solvers/NavierStokes/NSWithBody.h>
/**
* A fully discrete formulation of the direct forcing method
* first proposed by Fadlun et al (2000) 
* with modifications by Kim et al (2001)
*/
template <typename memoryType>
class DirectForcingSolver : public NSWithBody<memoryType>
{
protected:
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
	
public:
	DirectForcingSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	void initialise();
	void writeData();
	std::string name()
	{
		return "Direct Forcing";
	}
};
