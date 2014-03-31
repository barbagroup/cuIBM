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
void TairaColoniusSolver<memoryType>::calculateForce()
{
	typedef typename cusp::array1d<real, memoryType>::iterator ValueIterator;
	typedef typename cusp::array1d_view<ValueIterator>         View;
	
	int     nx = NavierStokesSolver<memoryType>::domInfo->nx,
	        ny = NavierStokesSolver<memoryType>::domInfo->ny,
	        numBodies = NSWithBody<memoryType>::B.numBodies,
	        totalPoints = NSWithBody<memoryType>::B.totalPoints,
	        ETRows = (nx-1)*ny + nx*(ny-1);
	real    dx, dy;
	View    f, fView;
	
	cusp::array1d<real, memoryType> F(ETRows), fTemp(2*totalPoints);
	
	dx = NavierStokesSolver<memoryType>::domInfo->dx[ NSWithBody<memoryType>::B.I[0] ],
	dy = dx;
	
	f = View(NavierStokesSolver<memoryType>::lambda.begin() + nx*ny, NavierStokesSolver<memoryType>::lambda.end());
	
	// loop through bodies
	for(int l=0; l < numBodies; l++)
	{
		dx = NavierStokesSolver<memoryType>::domInfo->dx[ NSWithBody<memoryType>::B.I[l] ],
		dy = dx;
		
		// x-component of the force
		fTemp = f;
		fView = View(fTemp.begin(), fTemp.begin() + NSWithBody<memoryType>::B.offsets[l]);
		cusp::blas::fill(fView, 0.0);
		if(l < numBodies-1)
		{
			fView = View(fTemp.begin() + NSWithBody<memoryType>::B.offsets[l+1], fTemp.end());
			cusp::blas::fill(fView, 0.0);
		}
		cusp::multiply(ET, fTemp, F);
		NSWithBody<memoryType>::B.forceX[l] = (dx*dy)/dx*thrust::reduce(F.begin(), F.end());
		
		// y-component of the force
		fTemp = f;
		fView = View(fTemp.begin(), fTemp.begin() + totalPoints + NSWithBody<memoryType>::B.offsets[l]);
		cusp::blas::fill(fView, 0.0);
		if(l < numBodies-1)
		{
			fView = View(fTemp.begin() + totalPoints + NSWithBody<memoryType>::B.offsets[l+1], fTemp.end());
			cusp::blas::fill(fView, 0.0);
		}
		cusp::multiply(ET, fTemp, F);
		NSWithBody<memoryType>::B.forceY[l] = (dx*dy)/dy*thrust::reduce(F.begin(), F.end());
	}
}
