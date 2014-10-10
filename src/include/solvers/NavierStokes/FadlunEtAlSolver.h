/**
* @file  FadlunEtAlSolver.h
* @brief Modified direct forcing immersed boundary method.
*/

#pragma once

#include <solvers/NavierStokes/DirectForcingSolver.h>
/**
* The direct forcing method proposed by Fadlun et al (2000)
* with modifications by Kim et al (2001)
*/
template <typename memoryType>
class FadlunEtAlSolver : public DirectForcingSolver<memoryType>
{
private:
	virtual void generateQT();
	
public:
	FadlunEtAlSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	virtual std::string name()
	{
		return "Fadlun et al.";
	}
};
