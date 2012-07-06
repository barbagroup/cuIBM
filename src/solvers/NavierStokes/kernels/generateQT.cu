#include <solvers/NavierStokes/kernels/generateQT.h>

namespace kernels
{

__global__
void updateQFadlun(int *QRows, int *QCols, real *QVals, int QSize, int *tags)
{
	int  I = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(I >= QSize) return;
	
	QVals[I] *= ( tags[QRows[I]] == -1 );
}

__global__
void updateQFadlun(int *QRows, int *QCols, real *QVals, int QSize, int *tagsX, int *tagsY)
{
	int  I = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(I >= QSize) return;
	
	QVals[I] *= ( tagsX[QRows[I]] == -1 && tagsY[QRows[I]] == -1 );
}

} // end of namespace kernels
