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

#include <solvers/NavierStokes/FadlunEtAlSolver.h>
#include <sys/stat.h>

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::initialise()
{
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
        ny = NavierStokesSolver<memoryType>::domInfo->ny;

	int numUV = (nx-1)*ny + nx*(ny-1);
	int numP  = nx*ny;
	
	NavierStokesSolver<memoryType>::initialiseCommon();
	
	NSWithBody<memoryType>::initialiseBodies();
	
	NavierStokesSolver<memoryType>::initialiseArrays(numUV, numP);
	
	NavierStokesSolver<memoryType>::logger.startTimer("allocateMemory");
/*	
	tags.resize(numUV);
	tagsD.resize(numUV);
	coeffs.resize(numUV);
	coeffsD.resize(numUV);
*/
	tagsX.resize(numUV);
	tagsXD.resize(numUV);
	tagsY.resize(numUV);
	tagsYD.resize(numUV);
	coeffsX.resize(numUV);
	coeffsXD.resize(numUV);
	coeffsY.resize(numUV);
	coeffsYD.resize(numUV);
	NavierStokesSolver<memoryType>::logger.startTimer("allocateMemory");
	
	tagPoints();
	std::cout << "Done tagging points!" << std::endl;
	
	NavierStokesSolver<memoryType>::assembleMatrices();
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::updateSolverState()
{	
/*	NavierStokesSolver<memoryType>::updateBoundaryConditions();
	if (B.bodiesMove)
	{
		updateBodies();
		//updateA();
	}
	
	if (NSWithBody<memoryType>::B.bodiesMove)
	{
		// update the locations of the body points
		NSWithBody<memoryType>::updateBodies();
		
		tagPoints();
		NavierStokesSolver<memoryType>::assembleMatrices();
		
		//NavierStokesSolver<memoryType>::logger.startTimer("preconditioner2");
		//NavierStokesSolver<memoryType>::PC2->update(NavierStokesSolver<memoryType>::C);
		//NavierStokesSolver<memoryType>::logger.stopTimer("preconditioner2");
	}*/
}

template <typename memoryType>
void FadlunEtAlSolver<memoryType>::assembleRHS1()
{
	NavierStokesSolver<memoryType>::assembleRHS1();
	
	NavierStokesSolver<memoryType>::logger.startTimer("updateRHS1");
	updateRHS1();
	NavierStokesSolver<memoryType>::logger.startTimer("updateRHS1");
}

#include "FadlunEtAl/tagPoints.inl"
#include "FadlunEtAl/generateL.inl"
#include "FadlunEtAl/updateRHS1.inl"
#include "FadlunEtAl/generateQT.inl"
#include "FadlunEtAl/calculateForce.inl"

template class FadlunEtAlSolver<host_memory>;
template class FadlunEtAlSolver<device_memory>;
