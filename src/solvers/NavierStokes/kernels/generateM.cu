/***************************************************************************//**
* \file generateM.cu
* \author Krishnan, A. (anush@bu.edu)
* \biref Definition of the kernels required to generate matrix M
*/

#include <solvers/NavierStokes/kernels/generateM.h>

/********************//**
* \namespace kernels
* \brief Contain all custom-written CUDA kernels
*/
namespace kernels
{

/**
* \brief Assemble the part of M that corresponds to the x-velocity values
*
* \param MRows raw pointer to the array storing the row indices of the elements of M
* \param MCols raw pointer to the array storing the column indices of the elements of M
* \param MVals raw pointer to the array storing the values of the elements of M
* \param MinvRows raw pointer to the array storing the row indices of the elements of M-inverse
* \param MinvCols raw pointer to the array storing the column indices of the elements of M-inverse
* \param MinvVals raw pointer to the array storing the values of the elements of M
* \param nx number of cells in the x-direction
* \param ny number of cells in the y-direction
* \param dx raw pointer to the array storing the cell widths in the x-direction
* \param dy raw pointer to the array storing the cell widths in the y-direction
* \param dt time increment
*/
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

/**
* \brief Assemble the part of M that corresponds to the y-velocity values
*
* \param MRows raw pointer to the array storing the row indices of the elements of M
* \param MCols raw pointer to the array storing the column indices of the elements of M
* \param MVals raw pointer to the array storing the values of the elements of M
* \param MinvRows raw pointer to the array storing the row indices of the elements of M-inverse
* \param MinvCols raw pointer to the array storing the column indices of the elements of M-inverse
* \param MinvVals raw pointer to the array storing the values of the elements of M
* \param nx number of cells in the x-direction
* \param ny number of cells in the y-direction
* \param dx raw pointer to the array storing the cell widths in the x-direction
* \param dy raw pointer to the array storing the cell widths in the y-direction
* \param dt time increment
*/
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
