/***************************************************************************//**
* \file  DirectForcingSolver.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Modified direct forcing immersed boundary method.
*/

#pragma once

#include <solvers/NavierStokes/NSWithBody.h>

/********************//**
* \class DirectForcingSolver
* \brief A fully discrete formulation of the direct forcing method
*
* first proposed by Fadlun et al (2000) \n
* with modifications by Kim et al (2001)
*
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
	
	/**
	* \brief
	*/
	void tagPoints();
	
	//void tagPoints(real *bx, real *by);
	
	/**
	* \brief
	*/
	void tagPoints(real *bx, real *by, real *uB, real *vB);

	/**
	* \brief
	*/
	void generateA(real alpha);
	
	//void updateA();

	/**
	* \brief
	*/
	void generateL();
	
	/**
	* \brief
	*/
	void generateQT(int *QTRows, int *QTCols, real *QTVals){} // check this!
	
	/**
	* \brief
	*/
	void generateQT();
	
	/**
	* \brief
	*/
	void updateQ();
	
	/**
	* \brief
	*/
	void assembleRHS1();
	
	/**
	* \brief
	*/
	void updateRHS1();

	/**
	* \brief
	*/
	void updateSolverState();
	
public:
	/**
	* \brief Constructor of the class DirectForcingSolver
	*/
	DirectForcingSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	
	/**
	* \brief
	*/
	void initialise();
	
	/**
	* \brief
	*/
	void writeData();
	
	std::string name()
	{
		return "Direct Forcing";
	}
};
