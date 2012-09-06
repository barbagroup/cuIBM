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

#include <solvers/NavierStokes/TairaColoniusSolver.h>
#include <sys/stat.h>

template <typename memoryType>
void TairaColoniusSolver<memoryType>::initialise()
{	
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
        ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	int numUV = (nx-1)*ny + nx*(ny-1);
	int numP  = nx*ny;
	
	NSWithBody<memoryType>::initialiseBodies();
	NavierStokesSolver<memoryType>::initialiseCommon();
	NavierStokesSolver<memoryType>::initialiseArrays(numUV, numP+2*NSWithBody<memoryType>::B.totalPoints);
	NavierStokesSolver<memoryType>::assembleMatrices();
	if(NSWithBody<memoryType>::B.totalPoints > 0)
	{
		generateE();
		FxX.resize(NSWithBody<memoryType>::B.numCellsX[0]);
		FxY.resize(NSWithBody<memoryType>::B.numCellsY[0]);
		FxU.resize( (NSWithBody<memoryType>::B.numCellsX[0]+1)*NSWithBody<memoryType>::B.numCellsY[0] );
	}
}

template <typename memoryType>
void TairaColoniusSolver<memoryType>::updateSolverState()
{
	NavierStokesSolver<memoryType>::logger.startTimer("updateSolverState");
	if (NSWithBody<memoryType>::B.bodiesMove)
	{
		
		NSWithBody<memoryType>::updateBodies();
		//updateQT();
		generateQT();
		generateE();
		NavierStokesSolver<memoryType>::generateC();
		NavierStokesSolver<memoryType>::PC2->update(NavierStokesSolver<memoryType>::C);
	}
	NavierStokesSolver<memoryType>::logger.startTimer("updateSolverState");
}

#include "TairaColonius/generateQT.inl"
#include "TairaColonius/generateBC2.inl"
#include "TairaColonius/calculateForce.inl"

template class TairaColoniusSolver<host_memory>;
template class TairaColoniusSolver<device_memory>;
