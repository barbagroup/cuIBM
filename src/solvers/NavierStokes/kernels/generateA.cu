#include <solvers/NavierStokes/kernels/generateA.h>

namespace kernels
{

__global__
void generateA(int *ARows, int *ACols, real *AVals, real *MVals, int *LRows, int *LCols, real *LVals, int ASize, real alpha)
{
	int  I = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(I >= ASize) return;

	ARows[I] = LRows[I];
	ACols[I] = LCols[I];
	AVals[I] = -alpha*LVals[I] + (LRows[I]==LCols[I])*MVals[LRows[I]];
}
	
} // end of namespace kernels