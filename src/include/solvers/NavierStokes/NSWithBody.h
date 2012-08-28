/*
*  Copyright (C) 2012 by Anush Krishnan, Simon Layton, Lorena Barba
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
* @file  NSWithBody.h
* @brief Generic Navier-Stokes solver in the presence of immersed boundaries
*/

#pragma once

#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <bodies.h>

/**
* @brief Generic Navier-Stokes solver in the presence of immersed boundaries
*/
template <typename memoryType>
class NSWithBody : public NavierStokesSolver<memoryType>
{
protected:
	bodies<memoryType> B; ///< Bodies in the flow

#if 0
	void initialiseBodies();
	void updateBodies();
#endif

#if 1
	/// Initialise the bodies
	void initialiseBodies()
	{
		parameterDB &db = *NavierStokesSolver<memoryType>::paramDB;
		B.initialise(db, *NavierStokesSolver<memoryType>::domInfo);
	}
	
	/// Update the body information at each time step during motion
	void updateBodies()
	{
		parameterDB &db = *NavierStokesSolver<memoryType>::paramDB;
		real dt   = db["simulation"]["dt"].get<real>();
		real Time = dt*(NavierStokesSolver<memoryType>::timeStep+1);
		B.update(db, *NavierStokesSolver<memoryType>::domInfo, Time);
	};
#endif
};
