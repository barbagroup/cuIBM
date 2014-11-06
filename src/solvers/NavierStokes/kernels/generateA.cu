/***************************************************************************//**
 * \file generateA.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the kernels required to generate the matrix
 *        resulting from the implicit terms in the momentum equation.
 */


#include "generateA.h"


/**
 * namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */
namespace kernels
{

/**
 * \brief Generates a block of the matrix resulting from implicit terms in the momentum equation.
 *
 * It assembles the matrix \c A as a combination 
 * of the Laplacian matrix \c L and the mass matrix \c M.
 *
 * \param ARows rows of the COO matrix \c A
 * \param ACols columns of the COO matrix \c A
 * \param AVals values of the COO matrix \c A
 * \param MVals values of the COO matrix \c M
 * \param LRows rows of the COO matrix \c L
 * \param LCols columns of the COO matrix \c L
 * \param LVals values of the COO matrix \c A
 * \param ASize number of entries of the COO matrix \c A
 * \param alpha implicit coefficient of the diffusive scheme
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
 * \brief Generates a block of the matrix resulting from implicit terms in the momentum equation
 *        for the direct forcing method.
 *
 * It assembles the matrix \c A as a combination 
 * of the Laplacian matrix \c L and the mass matrix \c M.
 *
 * \param ARows rows of the COO matrix \c A
 * \param ACols columns of the COO matrix \c A
 * \param AVals values of the COO matrix \c A
 * \param MVals values of the COO matrix \c M
 * \param LRows rows of the COO matrix \c L
 * \param LCols columns of the COO matrix \c L
 * \param LVals values of the COO matrix \c A
 * \param ASize number of entries of the COO matrix \c A
 * \param alpha implicit coefficient of the diffusive scheme
 * \param tagsX tag to check if the node is next to an immersed boundary
 * \param tagsY tag to check if the node is next to an immersed boundary
 */
__global__
void generateADirectForcing(int *ARows, int *ACols, real *AVals, real *MVals, int *LRows, int *LCols, real *LVals, int ASize, real alpha, int *tags)
{
	for(int I=threadIdx.x + blockIdx.x*blockDim.x; I<ASize; I += blockDim.x*gridDim.x)
	{
		ARows[I] = LRows[I];
		ACols[I] = LCols[I];
		AVals[I] =   (tags[LRows[I]] == -1)*(-alpha*LVals[I]) // if the current location is untagged, add -alpha*L
		           + (tags[LRows[I]] != -1)*(-LVals[I]) // if the current location is tagged, add -L
		           + (LRows[I]==LCols[I])*MVals[LRows[I]]; // if it is a diagonal, add M
	}
}
	
} // end of namespace kernels
