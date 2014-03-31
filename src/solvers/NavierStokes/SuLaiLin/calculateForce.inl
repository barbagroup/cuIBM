#include <solvers/NavierStokes/kernels/calculateForce.h>

template <typename memoryType>
void SuLaiLinSolver<memoryType>::calculateForce()
{
	int  nx = NavierStokesSolver<memoryType>::domInfo->nx,
	     ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	cusp::array1d<real, memoryType>
		tempF((nx-1)*ny + nx*(ny-1));
	
	real dx = NavierStokesSolver<memoryType>::domInfo->dx[ NSWithBody<memoryType>::B.I[0] ],
	     dy = dx;

	cusp::multiply(ET, f, tempF);
	NSWithBody<memoryType>::forceX = (dx*dy)/dx*thrust::reduce( tempF.begin(), tempF.begin()+(nx-1)*ny );
	NSWithBody<memoryType>::forceY = (dx*dy)/dy*thrust::reduce( tempF.begin()+(nx-1)*ny, tempF.end() );
}
