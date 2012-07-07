/**
*  Copyright (C) 2011 by Anush Krishnan, Simon Layton, Lorena Barba
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in
*  all copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*  THE SOFTWARE.
*/

/**
* \file
* \brief
*/
#pragma once

#include <types.h>

/**
* Class description
*/
class integrationScheme
{
public:
	int  subSteps;
	vecH gamma, // coefficent of the convection term in the current time step
	     zeta,  // coefficient of the convection  term in the previous time step 
	     alphaImplicit, // coefficient of the implicit diffusion term
	     alphaExplicit; // coefficient of the explicit diffusion term
	
	void initialise(timeScheme convScheme, timeScheme diffScheme)
	{
		std::cout << "Initialising integration scheme... ";
		switch(convScheme)
		{
			case EULER_EXPLICIT:
				subSteps = 1;
				break;
			case ADAMS_BASHFORTH_2:
				subSteps = 2;
				break;
			case RUNGE_KUTTA_3:
				subSteps = 3;
				break;
		}
		
		gamma.resize(subSteps);
		zeta.resize(subSteps);
		alphaExplicit.resize(subSteps);
		alphaImplicit.resize(subSteps);
		switch(convScheme)
		{
			case EULER_EXPLICIT:
				gamma[0] = 1;
				zeta[0]  = 0;
				break;
			case ADAMS_BASHFORTH_2:
				gamma[0] = 1.5;
				zeta[0]  = 0.5;
				break;
			case RUNGE_KUTTA_3:
				gamma[0] = 8.0/15;
				zeta[0]  = 0;
				gamma[1] = 5.0/12;
				zeta[1]  = -17.0/60;
				gamma[2] = 3.0/4;
				zeta[2]  = -5.0/12;
				break;
		}
		
		real aI, aE;
		switch(diffScheme)
		{
			case EULER_EXPLICIT:
				aI = 0.0;
				aE = 1.0;
				break;
			case EULER_IMPLICIT:
				aI = 1.0;
				aE = 0.0;
				break;
			case CRANK_NICOLSON:
				aI = 0.5;
				aE = 0.5;
		}
		for(int i=0; i<subSteps; i++)
		{
			alphaExplicit[i] = aE*(gamma[i]+zeta[i]);
			alphaImplicit[i] = aI*(gamma[i]+zeta[i]);
		}
		std::cout << "DONE! " << std::endl;
		for(int i=0; i<subSteps; i++)
		{
			std::cout << '[' <<  i << ']' << " " << gamma[i] << " " << zeta[i] << " " << alphaExplicit[i] << " " << alphaImplicit[i] << std::endl;
		}
	}
};
