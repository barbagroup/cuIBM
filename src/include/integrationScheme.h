/**
* @file integrationScheme.h
* @brief Specifies the time integration scheme used.
*/
#pragma once

#include <types.h>

/**
* @brief Specifies the time integration scheme used.
*/
class integrationScheme
{
public:
	int  subSteps;      ///< Number of substeps in each time step
	
	vecH gamma,         ///< Coefficent of the convection term in the current time step
	     zeta,          ///< Coefficient of the convection  term in the previous time step 
	     alphaImplicit, ///< Coefficient of the implicit diffusion term
	     alphaExplicit; ///< Coefficient of the explicit diffusion term
	
	/**
	* @brief Initialises the coefficients for the time-stepping schemes
	* @param convScheme Time-stepping scheme used for the convection term
	* @param diffScheme Time-stepping scheme used for the diffusion term
	*/
	void initialise(timeScheme convScheme, timeScheme diffScheme)
	{
		std::cout << "Initialising integration scheme... ";
		
		// Set the number of substeps required for the specified scheme
		switch(convScheme)
		{
			case EULER_EXPLICIT:
				subSteps = 1;
				break;
			case ADAMS_BASHFORTH_2:
				subSteps = 1;
				break;
			case RUNGE_KUTTA_3:
				subSteps = 3;
				break;
			default:
				break;
		}
		
		// Resize the arrays
		gamma.resize(subSteps);
		zeta.resize(subSteps);
		alphaExplicit.resize(subSteps);
		alphaImplicit.resize(subSteps);
		
		// Set the coefficients of the convectiont terms
		switch(convScheme)
		{
			case EULER_EXPLICIT:
				gamma[0] = 1;
				zeta[0]  = 0;
				break;
			case ADAMS_BASHFORTH_2:
				gamma[0] = 1.5;
				zeta[0]  = -0.5;
				break;
			case RUNGE_KUTTA_3:
				gamma[0] = 8.0/15;
				zeta[0]  = 0;
				gamma[1] = 5.0/12;
				zeta[1]  = -17.0/60;
				gamma[2] = 3.0/4;
				zeta[2]  = -5.0/12;
				break;
			default:
				break;
		}
		
		real aI, aE;
		// Set the coefficients of the diffusion terms
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
			default:
				break;
		}
		for(int i=0; i<subSteps; i++)
		{
			alphaExplicit[i] = aE*(gamma[i]+zeta[i]);
			alphaImplicit[i] = aI*(gamma[i]+zeta[i]);
		}	
		std::cout << "DONE! " << std::endl;
		
		// Print the coefficients
		for(int i=0; i<subSteps; i++)
		{
			std::cout << '[' <<  i << ']' << " " << gamma[i] << " " << zeta[i] << " " << alphaExplicit[i] << " " << alphaImplicit[i] << std::endl;
		}
	}
};
