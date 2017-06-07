/**
 * \file convectionTerm.h
 * \brief Unit-test class to calculate the convection term.
 */


#pragma once

#include "solvers/NavierStokesSolver.h"
#include "utilities/bodies.h"


/**
 * \class convectionTerm
 * \brief Generic Navier-Stokes solver.
 */
template <typename memoryType>
class convectionTerm : public NavierStokesSolver<memoryType>
{
private:
	cusp::array1d<real, memoryType>
	    HExact;
	void initialiseFluxes();
	void initialiseFluxes(real *q);
	void initialiseExactSolution();
	void initialiseExactSolution(real *HExact);

public:
	void initialise();
	void writeData();
	void shutDown();

	convectionTerm(parameterDB &paramDB, domain &domInfo)
	{
		this->paramDB = &paramDB;
		this->domInfo = &domInfo;
	}
};
