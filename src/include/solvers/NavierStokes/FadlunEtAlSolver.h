/***************************************************************************//**
* \file  FadlunEtAlSolver.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Modified direct forcing immersed boundary method
*/

#pragma once

#include <solvers/NavierStokes/DirectForcingSolver.h>

/********************//**
* \class FadlunEtAlSolver
* \brief Direct forcing method proposed by Fadlun et al (2000)
*
* with modifications by Kim et al (2001)
*
*/
template <typename memoryType>
class FadlunEtAlSolver : public DirectForcingSolver<memoryType>
{
private:
	/**
	* \brief
	*/
	void generateQT();
	
public:
	/**
	* \brief Constructor of the class FadlunEtAlSolver
	*/
	FadlunEtAlSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);
	
	std::string name()
	{
		return "Fadlun et al.";
	}
};
