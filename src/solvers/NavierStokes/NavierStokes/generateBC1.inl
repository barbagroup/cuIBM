#include <solvers/NavierStokes/kernels/generateBC1.h>

template <>
void NavierStokesSolver<device_memory>::generateBC1()
{
	logger.startTimer("generateBC1");
	
	real alpha = intgSchm.alphaImplicit[subStep];
	
	// raw pointers from cusp arrays
	real *bc1_r = thrust::raw_pointer_cast(&bc1[0]),
	     *q_r   = thrust::raw_pointer_cast(&q[0]),
	     *dx  = thrust::raw_pointer_cast(&(domInfo->dx[0])),
	     *dy  = thrust::raw_pointer_cast(&(domInfo->dy[0])),
	     *dxD = thrust::raw_pointer_cast(&(domInfo->dxD[0])),
	     *dyD = thrust::raw_pointer_cast(&(domInfo->dyD[0]));
	
	real *xminus = thrust::raw_pointer_cast(&(bc[XMINUS][0])),
	     *xplus  = thrust::raw_pointer_cast(&(bc[XPLUS][0])),
	     *yminus = thrust::raw_pointer_cast(&(bc[YMINUS][0])),
	     *yplus  = thrust::raw_pointer_cast(&(bc[YPLUS][0]));

	real dx0, dx1, dy0, dy1,
	     nu = (*paramDB)["flow"]["nu"].get<real>(),
         dt = (*paramDB)["simulation"]["dt"].get<real>();

	boundaryCondition **bcInfo = (*paramDB)["flow"]["boundaryConditions"].get<boundaryCondition **>();
	
	const int blocksize = 256;
	
	real C, beta;
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	
	// zero the bc1 vector
	cusp::blas::fill(bc1, 0.0);
	
	dim3 dimGridx( int((nx - 0.5)/blocksize) + 1, 1),
	     dimGridy( int((ny - 0.5)/blocksize) + 1, 1);
			
	dim3 dimBlock(blocksize, 1);
	
	/// bottom
	if(bcInfo[YMINUS][0].type == DIRICHLET)
	{
		/// u
		dy0	= 0.5*dy[0];
		dy1	= 0.5*(dy[0] + dy[1]);
		/// multiply by 0.5 for the Crank-Nicolson scheme and 2.0 for the non-uniform central difference
		C	= alpha * 2.0 * nu / (dy0 * (dy0+dy1));
		kernels::bc1DirichletU <<<dimGridx, dimBlock>>> (bc1_r, nx-1, nx, 0, 1, dxD, C, yminus);
		
		/// v
		C	= alpha * 2.0 * nu / (dy[0] * (dy[0]+dy[1]));
		kernels::bc1DirichletV <<<dimGridx, dimBlock>>> (bc1_r, nx, nx, numU, 0, 1, dyD, C, yminus, nx-1);
	}
	//else if(F.nbc.bottom_type==BC_CONVECTIVE)
	
	/// top
	if(bcInfo[YPLUS][0].type == DIRICHLET)
	{
		/// u
		dy0	= 0.5*(dy[ny-2] + dy[ny-1]);
		dy1	= 0.5*(dy[ny-1]);
		C	= alpha * 2.0 * nu / (dy1 * (dy0+dy1));
		kernels::bc1DirichletU <<<dimGridx, dimBlock>>> (bc1_r, nx-1, nx, (nx-1)*(ny-1), 1, dxD, C, yplus);
		
		/// v
		C	= alpha * 2.0 * nu / (dy[ny-1] * (dy[ny-1]+dy[ny-2]));
		kernels::bc1DirichletV <<<dimGridx, dimBlock>>> (bc1_r, nx, nx, numU, nx*(ny-2), 1, dyD, C, yplus, nx-1);
	}
	else if(bcInfo[YPLUS][0].type == CONVECTIVE)
	{
		/// u
		beta = bcInfo[YPLUS][1].value * dt / (0.5 * dy[ny-1]);
		dy0	= 0.5*(dy[ny-2] + dy[ny-1]);
		dy1	= 0.5*(dy[ny-1]);
		C	= alpha * 2.0 * nu / (dy1 * (dy0+dy1));
		kernels::bc1ConvectiveU <<<dimGridx, dimBlock>>> (bc1_r, nx-1, nx, (nx-1)*(ny-1), 1, dxD, dyD, C, yplus, q_r, beta);
		
		/// v
		beta = bcInfo[YPLUS][1].value * dt / dy[ny-1];
		C	= alpha * 2.0 * nu / (dy[ny-1] * (dy[ny-1]+dy[ny-2]));
		kernels::bc1ConvectiveV <<<dimGridx, dimBlock>>> (bc1_r, nx, nx, numU, nx*(ny-2), 1, dxD, dyD, C, yplus, nx-1, q_r, beta);
	}
	else if(bcInfo[YPLUS][0].type == SPECIAL)
	{
		/// u
		dy0	= 0.5*(dy[ny-2] + dy[ny-1]);
		dy1	= 0.5*(dy[ny-1]);
		C	= alpha * 2.0 * nu / (dy1 * (dy0+dy1));
		kernels::bc1SpecialU <<<dimGridx, dimBlock>>> (bc1_r, nx-1, nx, (nx-1)*(ny-1), 1, dxD, C, yplus, (timeStep+1)*dt);
		
		/// v
		C	= alpha * 2.0 * nu / (dy[ny-1] * (dy[ny-1]+dy[ny-2]));
		kernels::bc1DirichletV <<<dimGridx, dimBlock>>> (bc1_r, nx, nx, numU, nx*(ny-2), 1, dyD, C, yplus, nx-1);
	}
	
	/// left
	if(bcInfo[XMINUS][0].type == DIRICHLET)
	{
		/// u
		C = alpha * 2.0 * nu / ( dx[0] * (dx[0]+dx[1]) );
		kernels::bc1DirichletU <<<dimGridy, dimBlock>>> (bc1_r, ny, nx, 0, (nx-1), dxD, C, xminus);
		
		/// v
		dx0	= 0.5*dx[0];
		dx1	= 0.5*(dx[0] + dx[1]);
		C = alpha * 2.0 * nu / (dx0 * (dx0+dx1));
		kernels::bc1DirichletV <<<dimGridy, dimBlock>>> (bc1_r, ny-1, nx, numU, 0, nx, dyD, C, xminus, ny);
	}
	//else if(F.nbc.left_type==BC_CONVECTIVE)

	/// right
	if(bcInfo[XPLUS][0].type == DIRICHLET)
	{
		/// u
		C = alpha * 2.0 * nu / ( dx[nx-1] * (dx[nx-1]+dx[nx-2]) );
		kernels::bc1DirichletU <<<dimGridy, dimBlock>>> (bc1_r, ny, nx, nx-2, (nx-1), dxD, C, xplus);
		
		/// v
		dx0	= 0.5*(dx[nx-1] + dx[nx-2]);
		dx1	= 0.5*(dx[nx-1]);
		C	= alpha * 2.0 * nu / (dx1 * (dx0+dx1));
		kernels::bc1DirichletV <<<dimGridy, dimBlock>>> (bc1_r, ny-1, nx, numU, nx-1, nx, dyD, C, xplus, ny);
	}
	else if(bcInfo[XPLUS][0].type == CONVECTIVE)
	{
		/// u
		beta = bcInfo[XPLUS][0].value * dt / dx[nx-1];
		C = alpha * 2.0 * nu / ( dx[nx-1] * (dx[nx-1]+dx[nx-2]) );
		kernels::bc1ConvectiveU <<<dimGridy, dimBlock>>> (bc1_r, ny, nx, nx-2, nx-1, dxD, dyD, C, xplus, q_r, beta);
		
		/// v
		beta = bcInfo[XPLUS][0].value * dt / (0.5*dx[nx-1]);
		dx0	= 0.5*(dx[nx-1] + dx[nx-2]);
		dx1	= 0.5*(dx[nx-1]);
		C	= alpha * 2.0 * nu / (dx1 * (dx0+dx1));
		kernels::bc1ConvectiveV <<<dimGridy, dimBlock>>> (bc1_r, ny-1, nx, numU, nx-1, nx, dxD, dyD, C, xplus, ny, q_r, beta);
	}
	
	logger.stopTimer("generateBC1");
}

template <>
void NavierStokesSolver<host_memory>::generateBC1()
{
	logger.startTimer("generateBC1");
	
	real alpha = intgSchm.alphaImplicit[subStep];
	
	// raw pointers from cusp arrays
	real *dx  = thrust::raw_pointer_cast(&(domInfo->dx[0])),
	     *dy  = thrust::raw_pointer_cast(&(domInfo->dy[0]));
	
	real dx0, dx1, dy0, dy1,
	     nu = (*paramDB)["flow"]["nu"].get<real>(),
         dt = (*paramDB)["simulation"]["dt"].get<real>();

	boundaryCondition **bcInfo = (*paramDB)["flow"]["boundaryConditions"].get<boundaryCondition **>();
	
	real C,     // coefficient that the boundary velocity is multiplied with for the diffusion term
	     beta;  // coefficient that is used to convect the velocity to the boundary 
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	
	// zero the bc1 vector
	cusp::blas::fill(bc1, 0.0);
		
	/// bottom
	if(bcInfo[YMINUS][0].type == DIRICHLET)
	{
		// u
		dy0	= 0.5*dy[0];
		dy1	= 0.5*(dy[0] + dy[1]);
		C	= alpha * 2.0 * nu / (dy0 * (dy0+dy1));  // multiply by 0.5 for the Crank-Nicolson scheme and 2.0 for the non-uniform central difference
		for(int I=0; I<nx-1; I++)
		{
			bc1[I] += bc[YMINUS][I] * C * 0.5*(dx[I] + dx[I+1]);
		}
		
		/// v
		C	= alpha * 2.0 * nu / (dy[0] * (dy[0]+dy[1]));		
		for(int I=0; I<nx; I++)
		{	
			bc1[numU + I] += bc[YMINUS][I+(nx-1)] * C * 0.5*(dy[0] + dy[1]);
		}

	}
	//else if(F.nbc.bottom_type==BC_CONVECTIVE)
	
	/// top
	if(bcInfo[YPLUS][0].type == DIRICHLET)
	{
		// u
		dy0	= 0.5*(dy[ny-2] + dy[ny-1]);
		dy1	= 0.5*(dy[ny-1]);
		C	= alpha * 2.0 * nu / (dy1 * (dy0+dy1));		
		for(int I=0; I<nx-1; I++)
		{
			bc1[I+(nx-1)*(ny-1)] += bc[YPLUS][I] * C * 0.5*(dx[I] + dx[I+1]);
		}
		
		// v
		C	= alpha * 2.0 * nu / (dy[ny-1] * (dy[ny-1]+dy[ny-2]));
		for(int I=0; I<nx; I++)
		{	
			bc1[numU + I + nx*(ny-2)] += bc[YPLUS][I+(nx-1)] * C * 0.5*(dy[ny-1] + dy[ny-2]);
		}
	}
	else if(bcInfo[YPLUS][0].type == CONVECTIVE)
	{
		/// u
		beta = bcInfo[YPLUS][1].value * dt / (0.5 * dy[ny-1]);  // 0.5 because the velocity is being convected only half a cell
		dy0	= 0.5*(dy[ny-2] + dy[ny-1]);
		dy1	= 0.5*(dy[ny-1]);
		C	= alpha * 2.0 * nu / (dy1 * (dy0+dy1));
		for(int I=0; I<nx-1; I++)
		{
			bc[YPLUS][I] = (1.0-beta)*bc[YPLUS][I] + beta*q[I + (nx-1)*(ny-1)]/dy[ny-1];
			bc1[I + (nx-1)*(ny-1)] += bc[YPLUS][I] * C * 0.5*(dx[I] + dx[I+1]);
		}
		
		/// v
		beta = bcInfo[YPLUS][1].value * dt / dy[ny-1];
		C	= alpha * 2.0 * nu / (dy[ny-1] * (dy[ny-1]+dy[ny-2]));
		for(int I=0; I<nx; I++)
		{	
			bc[YPLUS][I+(nx-1)] = (1.0-beta)*bc[YPLUS][I+(nx-1)] + beta*q[numU + I + nx*(ny-2)]/dx[I];
			bc1[numU + I + nx*(ny-2)] += bc[YPLUS][I+(nx-1)] * C * 0.5*(dy[ny-2] + dy[ny-1]);
		}
	}
	
	/// left
	if(bcInfo[XMINUS][0].type == DIRICHLET)
	{
		/// u
		C = alpha * 2.0 * nu / ( dx[0] * (dx[0]+dx[1]) );
		for(int I=0; I<ny; I++)
		{
			bc1[I*(nx-1)] += bc[XMINUS][I] * C * 0.5*(dx[0] + dx[1]);
		}
		
		/// v
		dx0	= 0.5*dx[0];
		dx1	= 0.5*(dx[0] + dx[1]);
		C = alpha * 2.0 * nu / (dx0 * (dx0+dx1));
		for(int I=0; I<ny-1; I++)
		{
			bc1[numU + I*nx] += bc[XMINUS][I+ny] * C * 0.5*(dy[I] + dy[I+1]);
		}
	}
	//else if(F.nbc.left_type==BC_CONVECTIVE)

	/// right
	if(bcInfo[XPLUS][0].type == DIRICHLET)
	{
		/// u
		C = alpha * 2.0 * nu / ( dx[nx-1] * (dx[nx-1]+dx[nx-2]) );
		for(int I=0; I<ny; I++)
		{
			bc1[(nx-2) + I*(nx-1)] += bc[XPLUS][I] * C * 0.5*(dx[nx-2] + dx[nx-1]);
		}
		
		/// v
		dx0	= 0.5*(dx[nx-1] + dx[nx-2]);
		dx1	= 0.5*(dx[nx-1]);
		C	= alpha * 2.0 * nu / (dx1 * (dx0+dx1));
		for(int I=0; I<ny-1; I++)
		{
			bc1[numU + I*nx + (nx-1)] += bc[XPLUS][I+ny] * C * 0.5*(dy[I] + dy[I+1]);
		}
	}
	else if(bcInfo[XPLUS][0].type == CONVECTIVE)
	{
		/// u
		beta = bcInfo[XPLUS][0].value * dt / dx[nx-1];
		C = alpha * 2.0 * nu / ( dx[nx-1] * (dx[nx-1]+dx[nx-2]) );
		for(int I=0; I<ny; I++)
		{
			bc[XPLUS][I] = (1.0-beta)*bc[XPLUS][I] + beta*q[(nx-2) + I*(nx-1)]/dy[I];
			bc1[(nx-2) + I*(nx-1)] += bc[XPLUS][I] * C * 0.5*(dx[nx-2] + dx[nx-1]);
		}
		
		/// v
		beta = bcInfo[XPLUS][0].value * dt / (0.5*dx[nx-1]);
		dx0	= 0.5*(dx[nx-1] + dx[nx-2]);
		dx1	= 0.5*(dx[nx-1]);
		C	= alpha * 2.0 * nu / (dx1 * (dx0+dx1));
		for(int I=0; I<ny-1; I++)
		{
			bc[XPLUS][I+ny] = (1.0-beta)*bc[XPLUS][I+ny] + beta*q[numU + (nx-1) + I*nx]/dx[nx-1];
			bc1[numU + (nx-1) + I*nx] += bc[XPLUS][I+ny] * C * 0.5*(dy[I] + dy[I+1]);
		}
	}
	
	logger.stopTimer("generateBC1");
}
