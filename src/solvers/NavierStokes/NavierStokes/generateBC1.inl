template <>
void NavierStokesSolver<cooH, vecH>::generateBC1()
{
}

__global__
void bc1_dirichlet_u(real *bc1, int N, int nx, int offset, int stride, real *dx, real C, real *bc)
{
	int idx	= threadIdx.x + blockIdx.x*blockDim.x;
	
	if ( idx >= N ) return;	/// check if idx too high
	
	int	I	= (offset + idx*stride),
		i	= I % (nx-1);
	bc1[I] += bc[idx] * C * 0.5*(dx[i] + dx[i+1]);
}

__global__
void bc1_dirichlet_v(real *bc1, int N, int nx, int numU, int offset, int stride, real *dy, real C, real *bc, int numUbc)
{
	int idx	= threadIdx.x + blockIdx.x*blockDim.x;
		
	if ( idx >= N ) return;	/// check if idx too high
	
	int	I	= (offset + idx*stride),
		j	= I / nx;
	bc1[numU + I] += bc[idx+numUbc] * C * 0.5*(dy[j] + dy[j+1]);
}

__global__
void bc1_convective_u(real *bc1, int N, int nx, int offset, int stride, real *dx, real *dy, real C, real *bc, real *q, real alpha)
{
	int idx	= threadIdx.x + blockIdx.x*blockDim.x;
	
	if ( idx >= N ) return;	/// check if idx too high
	
	int	I = (offset + idx*stride),
		i = I % (nx-1),
		j = I / (nx-1);
	
	bc[idx] = (1.0-alpha)*bc[idx] + alpha*q[I]/dy[j];
	
	bc1[I] += bc[idx] * C * 0.5*(dx[i] + dx[i+1]);
}

__global__
void bc1_convective_v(real *bc1, int N, int nx, int numU, int offset, int stride, real *dx, real *dy, real C, real *bc, int numUbc, real *q, real alpha)
{
	int idx	= threadIdx.x + blockIdx.x*blockDim.x;
		
	if ( idx >= N ) return;	/// check if idx too high
	
	int	I = (offset + idx*stride),
		i = I % nx,
		j = I / nx;
	
	bc[idx+numUbc] = (1.0-alpha)*bc[idx+numUbc] + alpha*q[numU + I]/dx[i];
	
	bc1[numU + I] += bc[idx+numUbc] * C * 0.5*(dy[j] + dy[j+1]);
}

template <>
void NavierStokesSolver<cooD, vecD>::generateBC1()
{
	// raw pointers from cusp arrays
	real *bc1_r	= thrust::raw_pointer_cast(&bc1[0]),
	     *q_r   = thrust::raw_pointer_cast(&q[0]),
	     *dx    = thrust::raw_pointer_cast(&(domInfo->dx[0])),
	     *dy    = thrust::raw_pointer_cast(&(domInfo->dy[0])),
	     *dxD_r = thrust::raw_pointer_cast(&(domInfo->dxD[0])),
	     *dyD_r = thrust::raw_pointer_cast(&(domInfo->dyD[0]));
	
	//real *xminus = thrust::raw_pointer_cast(&(bcNP1[XMINUS][0])),
	//     *xplus  = thrust::raw_pointer_cast(&(bcNP1[XPLUS][0])),
	//     *yminus = thrust::raw_pointer_cast(&(bcNP1[YMINUS][0])),
	//     *yplus  = thrust::raw_pointer_cast(&(bcNP1[YPLUS][0]));
	
	real *xminus = thrust::raw_pointer_cast(&(bcN[XMINUS][0])),
	     *xplus  = thrust::raw_pointer_cast(&(bcN[XPLUS][0])),
	     *yminus = thrust::raw_pointer_cast(&(bcN[YMINUS][0])),
	     *yplus  = thrust::raw_pointer_cast(&(bcN[YPLUS][0]));

	real dx0, dx1, dy0, dy1,
	     nu = flowDesc->nu;
	
	const int blocksize = 256;
	
	real C, alpha, omega = 1.0;
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny,
	     numUV = numU + nx*(ny-1);
	
	dim3 dimGrid( int((numUV-0.5)/blocksize) + 1, 1);
	dim3 dimBlock(blocksize, 1);
	
	// call the kernels
	// zero the bc1 vector
	cusp::blas::fill(bc1, 0.0);
	
	dim3	dimGridx( int((nx - 0.5)/blocksize) + 1, 1),
			dimGridy( int((ny - 0.5)/blocksize) + 1, 1);
	
	/// bottom
	if(flowDesc->bcInfo[YMINUS][0].first == DIRICHLET)
	{
		/// u
		dy0	= 0.5*dy[0];
		dy1	= 0.5*(dy[0] + dy[1]);
		/// multiply by 0.5 for the Crank-Nicolson scheme and 2.0 for the non-uniform central difference
		C	= omega * 2.0 * nu / (dy0 * (dy0+dy1));
		bc1_dirichlet_u <<<dimGridx, dimBlock>>> (bc1_r, nx-1, nx, 0, 1, dxD_r, C, yminus);
		
		/// v
		C	= omega * 2.0 * nu / (dy[0] * (dy[0]+dy[1]));
		bc1_dirichlet_v <<<dimGridx, dimBlock>>> (bc1_r, nx, nx, numU, 0, 1, dyD_r, C, yminus, nx-1);
	}
	//else if(F.nbc.bottom_type==BC_CONVECTIVE)
	
	/// top
	if(flowDesc->bcInfo[YPLUS][0].first == DIRICHLET)
	{
		/// u
		dy0	= 0.5*(dy[ny-2] + dy[ny-1]);
		dy1	= 0.5*(dy[ny-1]);
		C	= omega * 2.0 * nu / (dy1 * (dy0+dy1));
		bc1_dirichlet_u <<<dimGridx, dimBlock>>> (bc1_r, nx-1, nx, (nx-1)*(ny-1), 1, dxD_r, C, yplus);
		
		/// v
		C	= omega * 2.0 * nu / (dy[ny-1] * (dy[ny-1]+dy[ny-2]));
		bc1_dirichlet_v <<<dimGridx, dimBlock>>> (bc1_r, nx, nx, numU, nx*(ny-2), 1, dyD_r, C, yplus, nx-1);
	}
	else if(flowDesc->bcInfo[YPLUS][0].first == CONVECTIVE)
	{
		/// u
		alpha = flowDesc->bcInfo[YPLUS][1].second * simPar->dt / (0.5 * dy[ny-1]);
		dy0	= 0.5*(dy[ny-2] + dy[ny-1]);
		dy1	= 0.5*(dy[ny-1]);
		C	= omega * 2.0 * nu / (dy1 * (dy0+dy1));
		bc1_convective_u <<<dimGridx, dimBlock>>> (bc1_r, nx-1, nx, (nx-1)*(ny-1), 1, dxD_r, dyD_r, C, yplus, q_r, alpha);
		
		/// v
		alpha = flowDesc->bcInfo[YPLUS][1].second * simPar->dt / dy[ny-1];
		C	= omega * 2.0 * nu / (dy[ny-1] * (dy[ny-1]+dy[ny-2]));
		bc1_convective_v <<<dimGridx, dimBlock>>> (bc1_r, nx, nx, numU, nx*(ny-2), 1, dxD_r, dyD_r, C, yplus, nx-1, q_r, alpha);
		
		//cudaMemcpy(F.nbc.top, F.nbc.top_d, (2*nx-1)*sizeof(real), cudaMemcpyDeviceToHost);
	}
	
	/// left
	if(flowDesc->bcInfo[XMINUS][0].first == DIRICHLET)
	{
		/// u
		C = omega * 2.0 * nu / ( dx[0] * (dx[0]+dx[1]) );
		bc1_dirichlet_u <<<dimGridy, dimBlock>>> (bc1_r, ny, nx, 0, (nx-1), dxD_r, C, xminus);
		
		/// v
		dx0	= 0.5*dx[0];
		dx1	= 0.5*(dx[0] + dx[1]);
		C = omega * 2.0 * nu / (dx0 * (dx0+dx1));
		bc1_dirichlet_v <<<dimGridy, dimBlock>>> (bc1_r, ny-1, nx, numU, 0, nx, dyD_r, C, xminus, ny);
	}
	//else if(F.nbc.left_type==BC_CONVECTIVE)

	/// right
	if(flowDesc->bcInfo[XPLUS][0].first == DIRICHLET)
	{
		/// u
		C = omega * 2.0 * nu / ( dx[nx-1] * (dx[nx-1]+dx[nx-2]) );
		bc1_dirichlet_u <<<dimGridy, dimBlock>>> (bc1_r, ny, nx, nx-2, (nx-1), dxD_r, C, xplus);
		
		/// v
		dx0	= 0.5*(dx[nx-1] + dx[nx-2]);
		dx1	= 0.5*(dx[nx-1]);
		C	= omega * 2.0 * nu / (dx1 * (dx0+dx1));
		bc1_dirichlet_v <<<dimGridy, dimBlock>>> (bc1_r, ny-1, nx, numU, nx-1, nx, dyD_r, C, xplus, ny);
	}
	else if(flowDesc->bcInfo[XPLUS][0].first == CONVECTIVE)
	{
		/// u
		alpha = flowDesc->bcInfo[XPLUS][0].second * simPar->dt / dx[nx-1];
		C = omega * 2.0 * nu / ( dx[nx-1] * (dx[nx-1]+dx[nx-2]) );
		bc1_convective_u <<<dimGridy, dimBlock>>> (bc1_r, ny, nx, nx-2, nx-1, dxD_r, dyD_r, C, xplus, q_r, alpha);
		
		/// v
		alpha = flowDesc->bcInfo[XPLUS][0].second * simPar->dt / (0.5*dx[nx-1]);
		dx0	= 0.5*(dx[nx-1] + dx[nx-2]);
		dx1	= 0.5*(dx[nx-1]);
		C	= omega * 2.0 * nu / (dx1 * (dx0+dx1));
		bc1_convective_v <<<dimGridy, dimBlock>>> (bc1_r, ny-1, nx, numU, nx-1, nx, dxD_r, dyD_r, C, xplus, ny, q_r, alpha);
		
		//cudaMemcpy(F.nbc.right, F.nbc.right_d, (2*ny-1)*sizeof(real), cudaMemcpyDeviceToHost);
	}
}