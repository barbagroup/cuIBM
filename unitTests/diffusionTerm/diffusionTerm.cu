#include "diffusionTerm.h"

//##############################################################################
//                              INITIALISE
//##############################################################################

template <typename memoryType>
void diffusionTerm<memoryType>::initialise()
{
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
	    ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	int numUV = (nx-1)*ny + nx*(ny-1);
	int numP  = nx*ny;
	
	NavierStokesSolver<memoryType>::initialiseCommon();
	NavierStokesSolver<memoryType>::initialiseArrays(numUV, numP);
	NavierStokesSolver<memoryType>::intgSchm.gamma[0] = 0.0;
	NavierStokesSolver<memoryType>::intgSchm.zeta[0]  = 0.0;
	NavierStokesSolver<memoryType>::generateRN();
	rnExact.resize(numUV);
	cusp::blas::fill(rnExact, 0.0);
	initialiseExactSolution();
}

/**
* \brief Sets the initial value of all the fluxes in the flow field
*/
template <>
void diffusionTerm<host_memory>::initialiseFluxes()
{
	real *q_r = thrust::raw_pointer_cast(&(q[0]));
	initialiseFluxes(q_r);
	qStar = q;
}

template<>
void diffusionTerm<device_memory>::initialiseFluxes()
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
void diffusionTerm<memoryType>::initialiseFluxes(real *q)
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
void diffusionTerm<host_memory>::initialiseExactSolution()
{
	real *rnExact_r = thrust::raw_pointer_cast(&(rnExact[0]));
	real nu = (*paramDB)["flow"]["nu"].get<real>();
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	initialiseExactSolution(rnExact_r, nu, dt);
}

template<>
void diffusionTerm<device_memory>::initialiseExactSolution()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	vecH rnExactHost((nx-1)*ny+nx*(ny-1));
	cusp::blas::fill(rnExactHost, 0.0);
	
	// creating raw pointers
	real *rnExactHost_r = thrust::raw_pointer_cast(&(rnExactHost[0]));
	real nu = (*paramDB)["flow"]["nu"].get<real>();
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	initialiseExactSolution(rnExactHost_r, nu, dt);
	rnExact = rnExactHost;
}

template <typename memoryType>
void diffusionTerm<memoryType>::initialiseExactSolution(real *rnExact, real nu, real dt)
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
	real pi = 4*atan(1.0), u=0.0, v=0.0;
	//parameterDB *pDB = NavierStokesSolver<memoryType>::paramDB;
	real alpha = NavierStokesSolver<memoryType>::intgSchm.alphaExplicit[0];
	std::cout << "dt: " << dt << '\n';
	std::cout << "alpha: " << alpha << '\n';
	std::cout << "nu: " << nu << '\n';
	std::cout << NavierStokesSolver<memoryType>::intgSchm.gamma[0] << " " << NavierStokesSolver<memoryType>::intgSchm.zeta[0] << std::endl;
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			u = sin(pi*xu[i])*sin(pi*yu[j]) + 1;
			rnExact[j*(nx-1) + i]    = (u/dt - alpha*nu*2.0*pi*pi*sin(pi*xu[i])*sin(pi*yu[j]))*0.5*(dx[i]+dx[i+1]);
		}
	}
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			v = sin(pi*xv[i])*sin(pi*yv[j]);
			rnExact[j*nx + i + numU] = (v/dt - alpha*nu*2.0*pi*pi*sin(pi*xv[i])*sin(pi*yv[j]))*0.5*(dy[j]+dy[j+1]);
		}
	}
}

template <typename memoryType>
void diffusionTerm<memoryType>::writeData()
{
	int  nx = NavierStokesSolver<memoryType>::domInfo->nx,
	     ny = NavierStokesSolver<memoryType>::domInfo->ny;
	cusp::array1d<real, memoryType> errorx((nx-1)*ny), errory(nx*(ny-1));

	typedef typename cusp::array1d<real, memoryType>::iterator ValueIterator;
	typedef typename cusp::array1d_view<ValueIterator>         View;

	View    rnx, rnExactx, rny, rnExacty;

	rnx = View(NavierStokesSolver<memoryType>::rn.begin(), NavierStokesSolver<memoryType>::rn.begin()+(nx-1)*ny);
	rny = View(NavierStokesSolver<memoryType>::rn.begin()+(nx-1)*ny, NavierStokesSolver<memoryType>::rn.end());

	rnExactx = View(rnExact.begin(), rnExact.begin()+(nx-1)*ny);
	rnExacty = View(rnExact.begin()+(nx-1)*ny, rnExact.end());

	cusp::blas::axpby(rnx, rnExactx, errorx, -1.0, 1.0);
	cusp::blas::axpby(rny, rnExacty, errory, -1.0, 1.0);

	real errorxNorm  = cusp::blas::nrm2(errorx),
	     erroryNorm  = cusp::blas::nrm2(errorx),
	     rnExactxNorm = cusp::blas::nrm2(rnExactx),
	     rnExactyNorm = cusp::blas::nrm2(rnExacty);
	std::cout << std::setw(11) << nx << " x " << ny << " Relative L-2 Error in rnx: " << errorxNorm/rnExactxNorm << std::endl;
	std::cout << std::setw(11) << nx << " x " << ny << " Relative L-2 Error in rny: " << erroryNorm/rnExactyNorm << std::endl;
}

template <typename memoryType>
void diffusionTerm<memoryType>::shutDown()
{
	NavierStokesSolver<memoryType>::iterationsFile.close();
}

template class diffusionTerm<host_memory>;
template class diffusionTerm<device_memory>;
