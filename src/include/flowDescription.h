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
	
	std::pair <bcType, real> bcInfo[4][2];
	
	int numBodies;
	std::vector<body> B;
	
	flowDescription(real NU=0.01, real IU=0.0, real IV=0.0)
	{
		nu = NU;
		initialU  = IU;
		initialV  = IV;
		/*
		bcInfo[XMINUS][0] = std::make_pair(DIRICHLET, 0.0);
		bcInfo[XMINUS][1] = std::make_pair(DIRICHLET, 0.0);
		bcInfo[XPLUS][0]  = std::make_pair(DIRICHLET, 0.0);
		bcInfo[XPLUS][1]  = std::make_pair(DIRICHLET, 0.0);
		bcInfo[YMINUS][0] = std::make_pair(DIRICHLET, 0.0);
		bcInfo[YMINUS][1] = std::make_pair(DIRICHLET, 0.0);
		bcInfo[YPLUS][0]  = std::make_pair(DIRICHLET, 1.0);
		bcInfo[YPLUS][1]  = std::make_pair(DIRICHLET, 0.0);
		*/
		bcInfo[XMINUS][0] = std::make_pair(DIRICHLET, 1.0);
		bcInfo[XMINUS][1] = std::make_pair(DIRICHLET, 0.0);
		bcInfo[XPLUS][0]  = std::make_pair(CONVECTIVE, 1.0);
		bcInfo[XPLUS][1]  = std::make_pair(CONVECTIVE, 0.0);
		bcInfo[YMINUS][0] = std::make_pair(DIRICHLET, 1.0);
		bcInfo[YMINUS][1] = std::make_pair(DIRICHLET, 0.0);
		bcInfo[YPLUS][0]  = std::make_pair(DIRICHLET, 1.0);
		bcInfo[YPLUS][1]  = std::make_pair(DIRICHLET, 0.0);
		numBodies = 1;
		B.resize(numBodies);
		for(int k=0; k<numBodies; k++)
		{
			B[k].numPoints = 63;
			B[k].X0_x = 0.0;
			B[k].X0_y = 0.0;
			B[k].Theta0 = 0.0;
			B[k].X  = new real[B[k].numPoints];
			B[k].Y  = new real[B[k].numPoints];
			for(int i=0; i<B[k].numPoints; i++)
			{
				B[k].X[i] = 0.5*cos(i*2*M_PI/B[k].numPoints);
				B[k].Y[i] = 0.5*sin(i*2*M_PI/B[k].numPoints);
			}
		}
	}
};