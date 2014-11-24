#include "convectionTerm.h"

//##############################################################################
//                              INITIALISE
//##############################################################################

template <typename memoryType>
void convectionTerm<memoryType>::initialise()
{
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
	    ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	int numUV = (nx-1)*ny + nx*(ny-1);
	int numP  = nx*ny;
	
	NavierStokesSolver<memoryType>::initialiseCommon();
	NavierStokesSolver<memoryType>::initialiseArrays(numUV, numP);
	HExact.resize(numUV);
	cusp::blas::fill(HExact, 0.0);
	initialiseExactSolution();
	cusp::blas::scal(HExact, -1.0/NavierStokesSolver<memoryType>::intgSchm.gamma[0]);
}

/**
* \brief Sets the initial value of all the fluxes in the flow field
*/
template <>
void convectionTerm<host_memory>::initialiseFluxes()
{
	real *q_r = thrust::raw_pointer_cast(&(q[0]));
	initialiseFluxes(q_r);
	qStar = q;
}

template<>
void convectionTerm<device_memory>::initialiseFluxes()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	vecH qHost((nx-1)*ny+nx*(ny-1));
	
	// creating raw pointers
	real *qHost_r = thrust::raw_pointer_cast(&(qHost[0]));
	initialiseFluxes(qHost_r);
	q = qHost;
	qStar=q;
}

template <typename memoryType>
void convectionTerm<memoryType>::initialiseFluxes(real *q)
{
	real *xu = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->xu[0])),
	     *yu = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->yu[0])),
	     *xv = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->xv[0])),
	     *yv = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->yv[0])),
	     *dx = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->dx[0])),
	     *dy = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->dy[0]));
	int  nx = NavierStokesSolver<memoryType>::domInfo->nx,
	     ny = NavierStokesSolver<memoryType>::domInfo->ny,
	     numU  = (nx-1)*ny;
	real pi = 4*atan(1.0);
	
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			q[j*(nx-1) + i] = (sin(pi*xu[i])*sin(pi*yu[j]) + 1) * dy[j];
		}
	}
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			q[j*nx + i + numU] = sin(pi*xv[i])*sin(pi*yv[j]) * dx[i];
		}
	}
}

/**
* \brief Sets the exact solution of the convection term
*/
template <>
void convectionTerm<host_memory>::initialiseExactSolution()
{
	real *HExact_r = thrust::raw_pointer_cast(&(HExact[0]));
	initialiseExactSolution(HExact_r);
}

template<>
void convectionTerm<device_memory>::initialiseExactSolution()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	vecH HExactHost((nx-1)*ny+nx*(ny-1));
	cusp::blas::fill(HExactHost, 0.0);
	
	// creating raw pointers
	real *HExactHost_r = thrust::raw_pointer_cast(&(HExactHost[0]));
	initialiseExactSolution(HExactHost_r);
	HExact = HExactHost;
}

template <typename memoryType>
void convectionTerm<memoryType>::initialiseExactSolution(real *HExact)
{
	real *xu = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->xu[0])),
	     *yu = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->yu[0])),
	     *xv = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->xv[0])),
	     *yv = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->yv[0]));
	int  nx = NavierStokesSolver<memoryType>::domInfo->nx,
	     ny = NavierStokesSolver<memoryType>::domInfo->ny,
	     numU  = (nx-1)*ny;
	real pi = 4*atan(1.0);
	
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			HExact[j*(nx-1) + i] =   2*pi*cos(pi*xu[i])*sin(pi*yu[j])*(sin(pi*xu[i])*sin(pi*yu[j])+1)
			                           + pi*sin(pi*xu[i])*cos(pi*yu[j])*(2*sin(pi*xu[i])*sin(pi*yu[j])+1);
		}
	}
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			HExact[j*nx + i + numU] =   pi*sin(pi*xv[i])*sin(pi*xv[i])*sin(2*pi*yv[j])
			                              + pi*cos(pi*xv[i])*sin(pi*yv[j])*(2*sin(pi*xv[i])*sin(pi*yv[j])+1);
		}
	}
}

template <typename memoryType>
void convectionTerm<memoryType>::writeData()
{
	int  nx = NavierStokesSolver<memoryType>::domInfo->nx,
	     ny = NavierStokesSolver<memoryType>::domInfo->ny;
	cusp::array1d<real, memoryType> errorx((nx-1)*ny), errory(nx*(ny-1));

	typedef typename cusp::array1d<real, memoryType>::iterator ValueIterator;
	typedef typename cusp::array1d_view<ValueIterator>         View;

	View    Hx, HExactx, Hy, HExacty;

	Hx = View(NavierStokesSolver<memoryType>::H.begin(), NavierStokesSolver<memoryType>::H.begin()+(nx-1)*ny);
	Hy = View(NavierStokesSolver<memoryType>::H.begin()+(nx-1)*ny, NavierStokesSolver<memoryType>::H.end());

	HExactx = View(HExact.begin(), HExact.begin()+(nx-1)*ny);
	HExacty = View(HExact.begin()+(nx-1)*ny, HExact.end());

	cusp::blas::axpby(Hx, HExactx, errorx, -1.0, 1.0);
	cusp::blas::axpby(Hy, HExacty, errory, -1.0, 1.0);

	real errorxNorm  = cusp::blas::nrm2(errorx),
	     erroryNorm  = cusp::blas::nrm2(errorx),
	     HExactxNorm = cusp::blas::nrm2(HExactx),
	     HExactyNorm = cusp::blas::nrm2(HExacty);
	std::cout << std::setw(11) << nx << " x " << ny << " Relative L-2 Error in Hx: " << errorxNorm/HExactxNorm << std::endl;
	std::cout << std::setw(11) << nx << " x " << ny << " Relative L-2 Error in Hy: " << erroryNorm/HExactyNorm << std::endl;
}

template <typename memoryType>
void convectionTerm<memoryType>::shutDown()
{
	NavierStokesSolver<memoryType>::iterationsFile.close();
}

template class convectionTerm<host_memory>;
template class convectionTerm<device_memory>;
