/***************************************************************************//**
* \file generateA.cu
* \author Krishnan, A. (anush@bu.edu)
* \brief Contain definition of the kernels required to generate the matrix A
*/

#include <solvers/NavierStokes/kernels/generateA.h>

/********************//**
* namespace kernels
* \brief Contain all the custom-written CUDA kernels
*/
namespace kernels
{

/**
* \brief To be documented
*/
__global__
void generateA(int *ARows, int *ACols, real *AVals, real *MVals, int *LRows, int *LCols, real *LVals, int ASize, real alpha)
{
	for (int I=threadIdx.x + blockIdx.x*blockDim.x; I<ASize; I += blockDim.x*gridDim.x)
	{
		ARows[I] = LRows[I];
		ACols[I] = LCols[I];
		AVals[I] = -alpha*LVals[I] + (LRows[I]==LCols[I])*MVals[LRows[I]];
	}
}

/**
* \brief To be documented
*/
__global__
void generateADirectForcing(int *ARows, int *ACols, real *AVals, real *MVals, int *LRows, int *LCols, real *LVals, int ASize, real alpha, int *tagsX, int *tagsY)
{
	for(int I=threadIdx.x + blockIdx.x*blockDim.x; I<ASize; I += blockDim.x*gridDim.x)
	{
		ARows[I] = LRows[I];
		ACols[I] = LCols[I];
		AVals[I] =   (tagsX[LRows[I]] == -1 && tagsY[LRows[I]] == -1)*(-alpha*LVals[I]) // if the current location is untagged, add -alpha*L
		           + (tagsX[LRows[I]] != -1 || tagsY[LRows[I]] != -1)*(-LVals[I]) // if the current location is tagged, add -L
		           + (LRows[I]==LCols[I])*MVals[LRows[I]]; // if it is a diagonal, add M
	}
}
	
} // end of namespace kernels
