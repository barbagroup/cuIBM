/**
* @file  TairaColoniusSolver.h
* @brief Solves the flow using the IB method described by Taira and Colonius (2007)
*
* <b>The immersed boundary method: a projection approach</b> \n
* Taira, K and Colonius, T \n
* Journal of Computational Physics \n
* Volume 225 Number 2 \n
* 2007
*/

#pragma once

#include <solvers/NavierStokes/NSWithBody.h>

/**
* @brief Immersed boundary method described by Taira and Colonius (2007)
*/
template <typename memoryType>
class TairaColoniusSolver : public NSWithBody<memoryType>
{
private:
	cusp::coo_matrix<int, real, memoryType> E, ET;
	
	void generateQT();
	void updateQT();
	void generateBC2();
	
	void updateSolverState();
	
	void calculateForce();
	void generateE();
	
public:
	void initialise();
	void writeData();
	std::string name()
	{
		return "Taira & Colonius";
	}
};
