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
