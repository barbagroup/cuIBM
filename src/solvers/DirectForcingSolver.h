/**
 * \file  DirectForcingSolver.h
 * \brief Declaration of the class \c DirectForcingSolver.
 */


#pragma once

#include "NSWithBody.h"


/**
 * \class DirectForcingSolver
 * \brief A fully discrete formulation of the direct forcing method.
 *
 * Based on the method first proposed by Fadlun et al (2000)
 * with modifications by Kim et al (2001).
 * It does not follow the same equations that they used, but use a 
 * fractional step method starting with the discretized equations,
 * based on the idea by Perot (1993).
 *
 * This method does not use an explicit pressure term in the step where
 * the intermediate velocity is calculated, and the pressure is directly 
 * obtained from the Poisson equation.
 */
template <typename memoryType>
class DirectForcingSolver : public NSWithBody<memoryType>
{
protected:
	// for the 1D interpolations
	cusp::array1d<int, host_memory>
		tags,   ///< vector used to differentiate forcing points from regular ones (host)
		tags2;  ///< vector used to differentiate forcing points from regular ones (host)
	cusp::array1d<int, device_memory>
		tagsD,   ///< vector used to differentiate forcing points from regular ones (device)
		tags2D;  ///< vector used to differentiate forcing points from regular ones (device)
	vecH coeffs,   ///< coefficients of interpolation (host)
	     coeffs2;  ///< other coefficients of interpolation; quadratic interpolation (host)
	vecD coeffsD,   ///< coefficients of interpolation (device)
	     coeffs2D;  ///< other coefficients of interpolation; quadratic interpolation (device)
	
	vecH uv;  ///< velocity field (host)
	vecD uvD; ///< velocity field (device)
	cusp::array1d<real, memoryType> pressure;  ///< pressure field
	
	// tag forcing points
	void tagPoints();
	void tagPoints(real *bx, real *by, real *uB, real *vB);

	// assemble the LHS matrix of the velocity system
	virtual void generateA(real alpha);

	// assemble the Laplacian operator
	virtual void generateL();
	
	// assemble the divergence operator
	void generateQT(int *QTRows, int *QTCols, real *QTVals){} // check this!
	virtual void generateQT();
	
	// update the gradient operator
	void updateQ();

	// compute the Poisson matrix
	virtual void generateC();
	
	// assemble the RHS of the velocity system
	virtual void assembleRHS1();
	// update the RHS of the velocity system
	void updateRHS1();

	// update the solver
	virtual void updateSolverState();

	// projection velocity onto divergence-free space
	virtual void projectionStep();

	// write info about mass flux
	void writeMassFluxInfo();
	
public:
	// constructor
	DirectForcingSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);

	// initialize the direct-forcing solver
	virtual void initialise();

	// write data
	virtual void writeData();
	
	/**
	 * \brief Returns the name of the solver as a string.
	 */
	virtual std::string name()
	{
		return "Direct Forcing";
	}

}; // DirectForcingSolver
