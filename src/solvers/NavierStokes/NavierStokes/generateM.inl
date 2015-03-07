/***************************************************************************//**
 * \file generateM.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods to generate the mass matrix and its inverse.
 */


#include <solvers/NavierStokes/kernels/generateM.h>


/**
 * \brief Generates the mass matrix and its inverse (on the host).
 *
 * The mass diagonal matrix is stored in a COO format. 
 * The time-increment value is included in the matrix.
 *
 */
template <>
void NavierStokesSolver<host_memory>::generateM()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  I;
	real value;
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	
	M.resize(numUV, numUV, numUV);
	Minv.resize(numUV, numUV, numUV);
	
	for (int j=0; j < ny; j++)
	{
		for (int i=0; i < nx-1; i++)
		{
			I = j*(nx-1) + i;
			value = 0.5*(domInfo->dx[i+1]+domInfo->dx[i]) / domInfo->dy[j] / dt;
			
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
			
			value  = 0.5*(domInfo->dy[j+1]+domInfo->dy[j]) / domInfo->dx[i] / dt;
			
			M.row_indices[I] = I;
			M.column_indices[I] = I;
			M.values[I] = value;
			
			Minv.row_indices[I] = I;
			Minv.column_indices[I] = I;
			Minv.values[I] = 1.0/value;
		}
	}
}

/**
 * \brief Generates the mass matrix and its inverse (on the device).
 *
 * The mass diagonal matrix is stored in a COO format. 
 * The time-increment value is included in the matrix.
 *
 */
template<>
void NavierStokesSolver<device_memory>::generateM()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;

	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	
	real *dxD = thrust::raw_pointer_cast(&(domInfo->dxD[0])),
	     *dyD = thrust::raw_pointer_cast(&(domInfo->dyD[0]));
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	M.resize(numUV, numUV, numUV);
	Minv.resize(numUV, numUV, numUV);
	
	int  *MRows = thrust::raw_pointer_cast(&(M.row_indices[0])),
	     *MCols = thrust::raw_pointer_cast(&(M.column_indices[0])),
	     *MinvRows = thrust::raw_pointer_cast(&(Minv.row_indices[0])),
	     *MinvCols = thrust::raw_pointer_cast(&(Minv.column_indices[0]));
	
	real *MVals = thrust::raw_pointer_cast(&(M.values[0])),
	     *MinvVals = thrust::raw_pointer_cast(&(Minv.values[0]));
	
	const int blockSize = 256;
	dim3 dimGrid( int((nx*ny-0.5)/blockSize) + 1, 1);
	dim3 dimBlock(blockSize, 1);
	kernels::fillM_u <<<dimGrid, dimBlock>>> (MRows, MCols, MVals, MinvRows, MinvCols, MinvVals, nx, ny, dxD, dyD, dt);
	kernels::fillM_v <<<dimGrid, dimBlock>>> (MRows, MCols, MVals, MinvRows, MinvCols, MinvVals, nx, ny, dxD, dyD, dt);
}
