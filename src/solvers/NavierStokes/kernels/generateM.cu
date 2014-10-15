#include "generateM.h"

namespace kernels
{

__global__
void fillM_u(int *MRows, int *MCols, real *MVals, int *MinvRows, int *MinvCols, real *MinvVals, int nx, int ny, real *dx, real *dy, real dt)
{
	int I = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(I >= (nx-1)*ny) return;
	
	int  i = I % (nx-1);
	int  j = I / (nx-1);
	real value = 0.5*(dx[i]+dx[i+1])/dy[j]/dt;
	
	MRows[I] = I;
	MCols[I] = I;
	MVals[I] = value;
	
	MinvRows[I] = I;
	MinvCols[I] = I;
	MinvVals[I] = 1.0/value;
}

__global__
void fillM_v(int *MRows, int *MCols, real *MVals, int *MinvRows, int *MinvCols, real *MinvVals, int nx, int ny, real *dx, real *dy, real dt)
{
	int I = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(I >= nx*(ny-1)) return;
	
	int  numU = (nx-1)*ny;
	int  i = I % nx;
	int  j = I / nx;
	real value = 0.5*(dy[j]+dy[j+1])/dx[i]/dt;
	
	MRows[I+numU] = I+numU;
	MCols[I+numU] = I+numU;
	MVals[I+numU] = value;
	
	MinvRows[I+numU] = I+numU;
	MinvCols[I+numU] = I+numU;
	MinvVals[I+numU] = 1.0/value;
}

} // end of namespace kernels
