/**
* @file  FadlunEtAlSolver.h
* @brief Modified direct forcing immersed boundary method.
*/

#pragma once

#include "DirectForcingSolver.h"
/**
* The direct forcing method proposed by Fadlun et al (2000)
* with modifications by Kim et al (2001)
*/
template <typename memoryType>
class FadlunEtAlSolver : public DirectForcingSolver<memoryType>
{
private:
	cusp::array1d<real, memoryType> pressure;

	virtual void generateQT();
	virtual void calculateExplicitLambdaTerms();
	virtual void projectionStep();
	
public:
	FadlunEtAlSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	virtual void initialise();
	virtual std::string name()
	{
		return "Fadlun et al.";
	}
};
