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
* @file  FadlunEtAlSolver.h
* @brief Modified direct forcing immersed boundary method.
*
* <b>Combined immersed-boundary finite-difference methods for three-dimensional complex flow simulations</b> \n
* Fadlun, E A, Verzicco, R, Orlandi, P, and Mohd-Yusof, J \n
* Journal of Computational Physics \n
* Volume 161 Number 1 \n
* 2000
*/

#pragma once

#include <solvers/NavierStokes/NSWithBody.h>
/**
* Immersed boundary method described by Fadlun et al (2000)
*/
template <typename memoryType>
class FadlunEtAlSolver : public NSWithBody<memoryType>
{
private:
/*
	// For the 1D interpolations
	cusp::array1d<int, host_memory>    tags;
	cusp::array1d<int, device_memory>  tagsD;
	vecH coeffs;
	vecD coeffsD;
*/	
	// For the 2D interpolations
	cusp::array1d<int, host_memory>    tagsX, tagsY;
	cusp::array1d<int, device_memory>  tagsXD, tagsYD;
	vecH coeffsX, coeffsY, uX, uY, vX, vY;
	vecD coeffsXD, coeffsYD, uXD, uYD, vXD, vYD;
	
	void tagPoints();
	void tagPoints(real *bx, real *by);

//	void generateA();
//	void updateA();

	void generateL();
	void generateQT(int *QTRows, int *QTCols, real *QTVals){} // check this!
	void generateQT();
	void updateQ();
	
	void assembleRHS1();
	void updateRHS1();

	void updateSolverState();
	void  calculateForceF();
	
public:
	void initialise();
	std::string name()
	{
		return "Fadlun et al";
	}
};
