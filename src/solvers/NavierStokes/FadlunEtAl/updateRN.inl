#include <solvers/NavierStokes/kernels/generateRN.h>

template <>
void FadlunEtAlSolver<device_memory>::updateRN()
{
	const int blocksize = 256;
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numUV = (nx-1)*ny + nx*(ny-1);
	
	real *rn_r   = thrust::raw_pointer_cast(&(rn[0]));
//	int  *tags_r = thrust::raw_pointer_cast(&(tagsD[0]));
	int  *tagsX_r = thrust::raw_pointer_cast(&(tagsXD[0])),\
	     *tagsY_r = thrust::raw_pointer_cast(&(tagsYD[0]));
	
	dim3 dimGrid( int((numUV-0.5)/blocksize) + 1, 1);
	dim3 dimBlock(blocksize, 1);
	
//	kernels::updateRN <<<dimGrid, dimBlock>>> (rn_r, numUV, tags_r);
	kernels::updateRN <<<dimGrid, dimBlock>>> (rn_r, numUV, tagsX_r, tagsY_r);
}

template <>
void FadlunEtAlSolver<host_memory>::updateRN()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numUV = (nx-1)*ny + nx*(ny-1);
	
	for(int i=0; i<numUV; i++)
	{
//		if(tags[i]!=-1)
		if(tagsX[i]!=-1 || tagsY[i]!=-1)
			rn[i] = 0.0;
	}
}
