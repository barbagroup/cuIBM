template <>
void NavierStokesSolver<cooH, vecH>::generateM()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  I;
	real value;
	
	M.resize(numUV, numUV, numUV);
	Minv.resize(numUV, numUV, numUV);
	
	for (int j=0; j < ny; j++)
	{
		for (int i=0; i < nx-1; i++)
		{
			I = j*(nx-1) + i;
			value = 0.5*(domInfo->dx[i+1]+domInfo->dx[i]) / domInfo->dy[j] / simPar->dt;
			
			M.row_indices[I] = I;
			M.column_indices[I] = I;
			M.values[I] = value;
			
			Minv.row_indices[I] = I;
			Minv.column_indices[I] = I;
			Minv.values[I] = 1.0/value;
		}
	}
	for (int j=0; j < ny-1; j++)
	{
		for (int i=0; i < nx; i++)
		{
			I = j*nx + i + numU;
			
			value  = 0.5*(domInfo->dy[j+1]+domInfo->dy[j]) / domInfo->dx[i] / simPar->dt;
			
			M.row_indices[I] = I;
			M.column_indices[I] = I;
			M.values[I] = value;
			
			Minv.row_indices[I] = I;
			Minv.column_indices[I] = I;
			Minv.values[I] = 1.0/value;
		}
	}
}

template <>
void NavierStokesSolver<cooD, vecD>::generateM()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  I;
	real value;
	
	cooH MHost(numUV, numUV, numUV),
	     MinvHost(numUV, numUV, numUV);
	
	for (int j=0; j < ny; j++)
	{
		for (int i=0; i < nx-1; i++)
		{
			I = j*(nx-1) + i;
			value = 0.5*(domInfo->dx[i+1]+domInfo->dx[i]) / domInfo->dy[j] / simPar->dt;
			
			MHost.row_indices[I] = I;
			MHost.column_indices[I] = I;
			MHost.values[I] = value;
			
			MinvHost.row_indices[I] = I;
			MinvHost.column_indices[I] = I;
			MinvHost.values[I] = 1.0/value;
		}
	}
	for (int j=0; j < ny-1; j++)
	{
		for (int i=0; i < nx; i++)
		{
			I = j*nx + i + numU;
			
			value  = 0.5*(domInfo->dy[j+1]+domInfo->dy[j]) / domInfo->dx[i] / simPar->dt;
			
			MHost.row_indices[I] = I;
			MHost.column_indices[I] = I;
			MHost.values[I] = value;
			
			MinvHost.row_indices[I] = I;
			MinvHost.column_indices[I] = I;
			MinvHost.values[I] = 1.0/value;
		}
	}
	
	M = MHost;
	Minv = MinvHost;
}
/*
template<>
void NavierStokesSolver<cooD, vecD>::generateM()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  numP  = nx*ny;
	M.resize(numUV, numUV, numUV);
	Minv.resize(numUV, numUV, numUV);
	int  I;
	real value;
	
	const int blockSize = 256;
	dim3 dimGrid( int((nx*ny-0.5)/blocksize) + 1, 1);
	dim3 dimBlock(blocksize, 1);
	fillM_u <<<dimGrid, dimBlock>>> (M, Minv, nx, ny, domInfo->dxD, domInfo->dyD, simPar->dt);
	fillM_v <<<dimGrid, dimBlock>>> (M, Minv, nx, ny, domInfo->dxD, domInfo->dyD, simPar->dt);
	
}

__global__
void fillM_u(cooD &M, cooD &Minv, int nx, int ny, real *dx, real *dy, real dt)
{
	int I = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(I >= (nx-1)*ny) return;
	
	int  i = I % (nx-1);
	int  j = I / (nx-1);
	real value = 0.5*(dx[i]+dx[i+1])/dy[j]/dt;
	
	M.row_indices[I] = I;
	M.row_indices[I] = I;
	M.values[I] = value;
	
	Minv.row_indices[I] = I;
	Minv.row_indices[I] = I;
	Minv.values[I] = 1.0/value;
}

__global__
void fillM_v(real *dx, real *dy, real dt)
{
	int I = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(I >= nx*(ny-1)) return;
	
	int  numU = (nx-1)*ny;
	int  i = I % nx;
	int  j = I / nx;
	real value = 0.5*(dy[j]+dx[j+1])/dx[i]/dt;
	
	M.row_indices[I+numU] = I+numU;
	M.row_indices[I+numU] = I+numU;
	M.values[I+numU] = value;
	
	Minv.row_indices[I+numU] = I+numU;
	Minv.row_indices[I+numU] = I+numU;
	Minv.values[I+numU] = 1.0/value;
}*/