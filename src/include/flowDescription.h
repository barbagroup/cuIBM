#pragma once

#include <types.h>
#include <body.h>
#include <utility>
/**
* Contains the description of the flow problem.
*/
class flowDescription
{
public:
	real nu;
	real initialU,
	     initialV;
	
	int numBodies;
	std::vector<body> B;
	std::pair <bcType, real> bcInfo[4][2];
	
	flowDescription(real NU=0.01, real IU=0.0, real IV=0.0, int NB=0)
	{
		nu = NU;
		initialU  = IU;
		initialV  = IV;
		numBodies = NB;
		bcInfo[XMINUS][0] = std::make_pair(DIRICHLET, 1.0);
		bcInfo[XMINUS][1] = std::make_pair(DIRICHLET, 0.0);
		bcInfo[XPLUS][0]  = std::make_pair(DIRICHLET, 0.0);
		bcInfo[XPLUS][1]  = std::make_pair(DIRICHLET, 0.0);
		bcInfo[YMINUS][0] = std::make_pair(DIRICHLET, 0.0);
		bcInfo[YMINUS][1] = std::make_pair(DIRICHLET, 0.0);
		bcInfo[YPLUS][0]  = std::make_pair(DIRICHLET, 0.0);
		bcInfo[YPLUS][1]  = std::make_pair(DIRICHLET, 0.0);
	}
};