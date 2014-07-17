/***************************************************************************//**
* \file generateVelB.cu
* \author Krishnan, A. (anush@bu.edu)
* \brief Definition of the CUDA kernels required to generate vector \c velB
*/

#include <solvers/NavierStokes/kernels/generateVelB.h>

/**
* \namespace kernels
* \brief Contain all custom-written CUDA kernels
*/
namespace kernels
{

/**
* \brief CUDA kernel to store element of x- and y- body velocity arrays
* 		into one single array
*
* \param velB vector that contains both x- and y- velocity components
* \param uB x-velocity of body points (all bodies included)
* \param vB y-velocity of body points (all bodies included)
* \param totalPoints number of body points (all bodies included)
*/
__global__
void fill_velB(real *velB, real *uB, real *vB, int totalPoints)
{
	int k = threadIdx.x + blockIdx.x*blockDim.x;
	if(k<totalPoints)
	{
		velB[k] = uB[k];
		velB[k + totalPoints] = vB[k];
	}
}

} // end of namespace kernels
