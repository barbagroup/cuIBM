#include <solvers/NavierStokes/kernels/calculateForce.h>

#define BSZ 16

template <>
void FadlunEtAlSolver<host_memory>::calculateForceF()
{
}

template <>
void FadlunEtAlSolver<device_memory>::calculateForceF()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	     
	real alpha = intgSchm.alphaImplicit[subStep],
	     nu = (*paramDB)["flow"]["nu"].get<real>();
	     
	cusp::array1d<real, device_memory>
	     f((nx-1)*ny + nx*(ny-1), 0.0),
	     temp((nx-1)*ny + nx*(ny-1));
	     
	// raw pointers for cup arrays
	real *f_r  = thrust::raw_pointer_cast(&f[0]),
	     *q_r  = thrust::raw_pointer_cast(&q[0]),
	     *rn_r = thrust::raw_pointer_cast(&rn[0]),
	     *dxD = thrust::raw_pointer_cast(&(domInfo->dxD[0])),
	     *dyD = thrust::raw_pointer_cast(&(domInfo->dyD[0]));
	
	int  *tagsX_r = thrust::raw_pointer_cast(&(tagsXD[0])),\
	     *tagsY_r = thrust::raw_pointer_cast(&(tagsYD[0]));
	
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	
	real h = 1;//domInfo->dx[ B.I[0] ];
	     
	dim3 dimGridx( int( (nx-1-0.5)/(BSZ-2) ) + 1, int( (ny-0.5)/(BSZ-2) ) + 1 ),
	     dimGridy( int( (nx-0.5)/(BSZ-2) ) + 1, int( (ny-1-0.5)/(BSZ-2) ) + 1 );
	dim3 dimBlock(BSZ, BSZ);
	
	// call the kernel
	kernels::forceX <<<dimGridx, dimBlock>>> (f_r, q_r, rn_r, tagsX_r, tagsY_r, nx, ny, dxD, dyD, dt, alpha, nu);
	
	//cusp::multiply(Q, lambda, temp);
	//cusp::blas::axpy(temp, f, 1.0);
	
	forceX = (h*h)*thrust::reduce( f.begin(), f.begin()+(nx-1)*ny );
	forceY = (h*h)*thrust::reduce( f.begin()+(nx-1)*ny, f.end() );
}
