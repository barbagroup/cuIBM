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

#include <solvers/NavierStokes/kernels/calculateForce.h>

template <typename memoryType>
void TairaColoniusSolver<memoryType>::calculateForceTC()
{
	int  nx = NavierStokesSolver<memoryType>::domInfo->nx,
	     ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	cusp::array1d<real, memoryType>
		f(2*NSWithBody<memoryType>::B.totalPoints),
		F((nx-1)*ny + nx*(ny-1));
	
	real dx = NavierStokesSolver<memoryType>::domInfo->dx[ NSWithBody<memoryType>::B.I[0] ],
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
	//NSWithBody<memoryType>::calculateForceCV();
}
