/***************************************************************************//**
 * \file generateVelB.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the kernels to generate body-velocities.
 */


#include "generateVelB.h"


/**
 * \namespace kernels
 * \brief Contains all custom-written CUDA kernels.
 */
namespace kernels
{

/**
 * \brief Stores an element of the u- and v- body-velocities into one single array.
 *
 * \param velB vector that contains both u- and v- velocities
 * \param uB u-velocity of body points (all bodies included)
 * \param vB v-velocity of body points (all bodies included)
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
