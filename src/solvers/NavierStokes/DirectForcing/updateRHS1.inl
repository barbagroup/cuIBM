/***************************************************************************//**
 * \file updateRHS1.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Updates the RHS.
 */


#include <solvers/NavierStokes/kernels/updateRHS1.h>

/**
 * \brief To be documented.
 */
template <>
void DirectForcingSolver<device_memory>::updateRHS1()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	int  numUV = (nx-1)*ny + nx*(ny-1);
	
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	
	real *rhs1_r = thrust::raw_pointer_cast(&(rhs1[0]));
	int  *tags_r = thrust::raw_pointer_cast(&(tagsD[0]));

	real *dx_r = thrust::raw_pointer_cast( &(domInfo->dxD[0]) ),
	     *dy_r = thrust::raw_pointer_cast( &(domInfo->dyD[0]) ),
	     *coeffs_r = thrust::raw_pointer_cast( &(coeffsD[0]) ),
	     *uv_r   = thrust::raw_pointer_cast( &(uvD[0]) );

	const int blocksize = 256;
	dim3 dimGrid( int( ((nx-1)*ny+nx*(ny-1)-0.5)/blocksize ) + 1, 1);
	dim3 dimBlock(blocksize, 1);
	
	// 1-d interpolation
	kernels::updateRHS1 <<<dimGrid, dimBlock>>> (rhs1_r, numUV, tags_r);
	
	kernels::updateRHS1X <<<dimGrid, dimBlock>>>(rhs1_r, nx, ny, dt, dx_r, tags_r, coeffs_r, uv_r);
	kernels::updateRHS1Y <<<dimGrid, dimBlock>>>(rhs1_r, nx, ny, dt, dy_r, tags_r, coeffs_r, uv_r);
}

/**
 * \brief To be documented.
 */
template <>
void DirectForcingSolver<host_memory>::updateRHS1()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny,
	     I = 0;
	
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			I = j*(nx-1)+i;

			rhs1[I] = (tags[I]==-1)*rhs1[I]
			          + ((tags[I]!=-1)*(1.0-coeffs[I])*uv[I]) * 0.5*(domInfo->dx[i+1]+domInfo->dx[i])/dt;
		}
	}
	
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			I = numU + j*nx + i;
		
			rhs1[I] = (tags[I]==-1)*rhs1[I]
			          + ((tags[I]!=-1)*((1.0-coeffs[I])*uv[I])) * 0.5*(domInfo->dy[j+1]+domInfo->dy[j])/dt;
		}
	}
}
