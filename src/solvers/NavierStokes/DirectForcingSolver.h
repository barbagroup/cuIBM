/***************************************************************************//**
 * \file  DirectForcingSolver.h
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Modified direct forcing immersed boundary method.
 */


#pragma once

#include "NSWithBody.h"


/**
 * \class DirectForcingSolver
 * \brief A fully discrete formulation of the direct forcing method.
 *
 * First proposed by Fadlun et al (2000) \n
 * with modifications by Kim et al (2001).
 *
 */
template <typename memoryType>
class DirectForcingSolver : public NSWithBody<memoryType>
{
protected:
	// For the 1D interpolations
	cusp::array1d<int, host_memory>    tags;
	cusp::array1d<int, device_memory>  tagsD;
	vecH coeffs, uv;
	vecD coeffsD, uvD;
	
	/**
	 * \brief To be documented.
	 */
	void tagPoints();
	
	//void tagPoints(real *bx, real *by);
	
	/**
	 * \brief To be documented.
	 */
	void tagPoints(real *bx, real *by, real *uB, real *vB);

	virtual void generateA(real alpha);

	virtual void generateL();
	void generateQT(int *QTRows, int *QTCols, real *QTVals){} // check this!
	virtual void generateQT();
	void updateQ();
	
	virtual void assembleRHS1();
	void updateRHS1();

	virtual void updateSolverState();

	void writeMassFluxInfo();
	
public:
	/**
	 * \brief Constructor of the class \c DirectForcingSolver.
	 */
	DirectForcingSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	virtual void initialise();
	virtual void writeData();
	virtual std::string name()
	{
		return "Direct Forcing";
	}
};
