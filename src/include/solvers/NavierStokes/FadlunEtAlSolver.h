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
* \file FadlunEtAlSolver.h
* \brief Direct forcing immersed boundary method.
*/

#pragma once

#include <solvers/NavierStokes/NavierStokesSolver.h>
#include <bodies.h>
/**
* Immersed boundary method described by Fadlun et al (2000)
*/
template <typename memoryType>
class FadlunEtAlSolver : public NavierStokesSolver<memoryType>
{
private:
	bodies<memoryType> B;
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
	vecH coeffsX, coeffsY;
	vecD coeffsXD, coeffsYD;
	
	void tagPoints();
	void tagPoints(real *bx, real *by);

//	void generateA();
//	void updateA();

	void generateL();
	void generateQT(int *QTRows, int *QTCols, real *QTVals){}
	void generateQT();
	void updateQFadlun();
	
	void generateRN();
	void updateRN();
	void generateBC1();
	void updateBC1();

	void initialiseBodies();
	void updateBodies();
//	void updateSolverState();
	
public:
	void initialise();
	std::string name()
	{
		return "Fadlun et al";
	}
};
