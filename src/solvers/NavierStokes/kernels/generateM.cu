/***************************************************************************//**
 * \file generateM.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the kernels to generate the mass matrix and its inverse.
 */


#include "generateM.h"


/**
 * \namespace kernels
 * \brief Contains all custom-written CUDA kernels.
 */
namespace kernels
{

/**
 * \brief Computes an element of the mass matrix and its inverse for a x-velocity node.
 *
 * \param MRows row index of elements of the mass matrix
 * \param MCols column index of elements of the mass matrix
 * \param MVals value of elements of the mass matrix
 * \param MinvRows row index of elements of the mass matrix inverse
 * \param MinvCols column index of elements of the mass matrix inverse
 * \param MinvVals value of elements of the mass matrix inverse
 * \param nx number of cells in the x-direction
 * \param ny number of cells in the y-direction
 * \param dx cell-widths in the x-direction
 * \param dy cell-widths in the y-direction
 * \param dt time-increment
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
 * \brief Computes an element of the mass matrix and its inverse for a y-velocity node.
 *
 * \param MRows row index of elements of the mass matrix
 * \param MCols column index of elements of the mass matrix
 * \param MVals value of elements of the mass matrix
 * \param MinvRows row index of elements of the mass matrix inverse
 * \param MinvCols column index of elements of the mass matrix inverse
 * \param MinvVals value of elements of the mass matrix inverse
 * \param nx number of cells in the x-direction
 * \param ny number of cells in the y-direction
 * \param dx cell-widths in the x-direction
 * \param dy cell-widths in the y-direction
 * \param dt time-increment
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
