/**
* @file  convectionTerm.h
* @brief Unit test class to calculate the convection term.
*/

#pragma once

#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <bodies.h>

/**
* @brief Generic Navier-Stokes solver in the presence of immersed boundaries
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
