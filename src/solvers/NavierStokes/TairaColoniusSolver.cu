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

#include <solvers/NavierStokes/TairaColoniusSolver.h>
#include <sys/stat.h>

template <typename memoryType>
void TairaColoniusSolver<memoryType>::initialise()
{	
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
        ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	int numUV = (nx-1)*ny + nx*(ny-1);
	int numP  = nx*ny;
	
	initialiseBodies();
	NavierStokesSolver<memoryType>::initialiseCommon();
	NavierStokesSolver<memoryType>::initialiseArrays(numUV, numP+2*B.totalPoints);
	NavierStokesSolver<memoryType>::assembleMatrices();
	if(B.totalPoints > 0)
	{
		generateE();
		FxX.resize(B.numCellsX[0]);
		FxY.resize(B.numCellsY[0]);
		FxU.resize( (B.numCellsX[0]+1)*B.numCellsY[0] );
	}
}

template <typename memoryType>
void TairaColoniusSolver<memoryType>::updateBodies()
{
}

template <typename memoryType>
void TairaColoniusSolver<memoryType>::updateSolverState()
{
	NavierStokesSolver<memoryType>::updateBoundaryConditions();
	if (B.bodiesMove) {
		updateBodies();
		updateQT();
		NavierStokesSolver<memoryType>::generateC();
	}
}

template <typename memoryType>
void TairaColoniusSolver<memoryType>::calculateForceTC()
{
	int  nx = NavierStokesSolver<memoryType>::domInfo->nx,
	     ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	array1d<real, memoryType>
		f(2*B.totalPoints),
		F((nx-1)*ny + nx*(ny-1));
	
	real dx = NavierStokesSolver<memoryType>::domInfo->dx[ B.I[0] ],
	     dy = dx;
	
	thrust::copy(NavierStokesSolver<memoryType>::lambda.begin() + nx*ny, NavierStokesSolver<memoryType>::lambda.end(), f.begin());
	cusp::multiply(ET, f, F);
	NavierStokesSolver<memoryType>::forceX = (dx*dy)/dx*thrust::reduce( F.begin(), F.begin()+(nx-1)*ny );
	NavierStokesSolver<memoryType>::forceY = (dx*dy)/dy*thrust::reduce( F.begin()+(nx-1)*ny, F.end() );
}

template <typename memoryType>
void TairaColoniusSolver<memoryType>::calculateForce()
{
	calculateForceTC();
	//calculateForce1();
}

#include "TairaColonius/generateQT.inl"
#include "TairaColonius/generateBC2.inl"
#include "TairaColonius/initialiseBodies.inl"
#include "TairaColonius/calculateForce.inl"

template class TairaColoniusSolver<host_memory>;
template class TairaColoniusSolver<device_memory>;
