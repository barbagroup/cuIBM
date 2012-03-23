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