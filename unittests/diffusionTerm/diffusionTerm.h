/**
 * \file diffusionTerm.h
 * \brief Unit-test class to calculate the diffusive term.
 */


#pragma once

#include "solvers/NavierStokesSolver.h"
#include "utilities/bodies.h"


/**
 * \class diffusionTerm
 * \brief Generic Navier-Stokes solver.
 */
template <typename memoryType>
class diffusionTerm : public NavierStokesSolver<memoryType>
{
private:
	cusp::array1d<real, memoryType>
	    rnExact;
	void initialiseFluxes();
	void initialiseFluxes(real *q);
	void initialiseExactSolution();
	void initialiseExactSolution(real *HExact, real nu, real dt);

public:
	void initialise();
	void writeData();
	void shutDown();

	diffusionTerm(parameterDB &paramDB, domain &domInfo)
	{
		this->paramDB = &paramDB;
		this->domInfo = &domInfo;
	}
};
