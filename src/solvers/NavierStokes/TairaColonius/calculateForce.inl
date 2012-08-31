/*
*  Copyright (C) 2012 by Anush Krishnan, Simon Layton, Lorena Barba
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in
*  all copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*  THE SOFTWARE.
*/

#include <solvers/NavierStokes/kernels/calculateForce.h>

template <>
void TairaColoniusSolver<host_memory>::calculateForce1()
{
//	int  nx = domInfo->nx, \
	     ny = domInfo->ny;
}

template <>
void TairaColoniusSolver<device_memory>::calculateForce1()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	
	real *FxX_r = thrust::raw_pointer_cast(&FxX[0]),
	     *FxY_r = thrust::raw_pointer_cast(&FxY[0]),
	     *FxU_r = thrust::raw_pointer_cast(&FxU[0]),
		 *q_r   = thrust::raw_pointer_cast(&q[0]),
		 *qOld_r   = thrust::raw_pointer_cast(&qOld[0]),
		 *lambda_r = thrust::raw_pointer_cast(&lambda[0]),
		 *dxD = thrust::raw_pointer_cast(&(domInfo->dxD[0])),
		 *dyD = thrust::raw_pointer_cast(&(domInfo->dyD[0]));

	const int blockSize = 256;
	dim3 dimGrid( int((B.numCellsX[0]+B.numCellsY[0]-0.5)/blockSize)+1, 1 );
	dim3 dimBlock(blockSize, 1);

	kernels::dragLeftRight <<<dimGrid, dimBlock>>> (FxY_r, q_r, lambda_r, dxD, dyD, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);
	kernels::dragBottomTop <<<dimGrid, dimBlock>>> (FxX_r, q_r, dxD, dyD, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	dim3 dimGrid2( int( ( (B.numCellsX[0]+1)*B.numCellsY[0]-0.5 )/blockSize )+1, 1 );
	
	kernels::dragUnsteady <<<dimGrid2, dimBlock>>> (FxU_r, q_r, qOld_r, dxD, dyD, dt, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);
	force1 = thrust::reduce(FxX.begin(), FxX.end()) + thrust::reduce(FxY.begin(), FxY.end()) + thrust::reduce(FxU.begin(), FxU.end());
}
