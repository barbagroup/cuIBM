/***************************************************************************//**
* \file generateVelB.cu
* \author Krishnan, A. (anush@bu.edu)
* \brief Definition of the kernels required to generate vector VelB
*/

#include <solvers/NavierStokes/kernels/generateVelB.h>

/********************//**
* \namespace kernels
* \brief Contain all custom-written CUDA kernels
*/
namespace kernels
{

/**
* \brief To be documented
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
