#include "generateVelB.h"
namespace kernels
{

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

}
