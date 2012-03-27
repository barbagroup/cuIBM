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
	
	flowDescription()
	{
		nu = 0.01;
		initialU  = 0.0;
		initialV  = 0.0;
		
		bcInfo[XMINUS][0] = std::make_pair(DIRICHLET, 0.0);
		bcInfo[XMINUS][1] = std::make_pair(DIRICHLET, 0.0);
		bcInfo[XPLUS][0]  = std::make_pair(DIRICHLET, 0.0);
		bcInfo[XPLUS][1]  = std::make_pair(DIRICHLET, 0.0);
		bcInfo[YMINUS][0] = std::make_pair(DIRICHLET, 0.0);
		bcInfo[YMINUS][1] = std::make_pair(DIRICHLET, 0.0);
		bcInfo[YPLUS][0]  = std::make_pair(DIRICHLET, 1.0);
		bcInfo[YPLUS][1]  = std::make_pair(DIRICHLET, 0.0);
		
		numBodies = 0;
		
		/*
		nu = 0.025;
		initialU  = 1.0;
		initialV  = 0.0;
		
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
			B[k].X.resize(B[k].numPoints);
			B[k].Y.resize(B[k].numPoints);
			for(int i=0; i<B[k].numPoints; i++)
			{
				B[k].X[i] = 0.5*cos(i*2*M_PI/B[k].numPoints);
				B[k].Y[i] = 0.5*sin(i*2*M_PI/B[k].numPoints);
			}
			B[k].tFlag = false;
			B[k].rFlag = false;
			B[k].VelX = 0.0;
			B[k].VelY = 0.0;
			B[k].Omg  = 0.0;
			B[k].AmpX = 0.0;
			B[k].OmegaX = 0.0;
			B[k].PhiX = 0.0;
			B[k].AmpY = 0.0;
			B[k].OmegaY = 0.0;
			B[k].PhiY = 0.0;
			B[k].AmpT = 0.0;
			B[k].OmegaT = 0.0;
			B[k].PhiT = 0.0;
		}
		*/
	}
};