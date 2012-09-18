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

#include <solvers/NavierStokes/kernels/generateQT.h>

namespace kernels
{

__global__
void updateQ(int *QRows, int *QCols, real *QVals, int QSize, int *tags)
{
	int I = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(I >= QSize) return;
	
	QVals[I] *= ( tags[QRows[I]] == -1 );
}

__global__
void updateQ(int *QRows, int *QCols, real *QVals, int QSize, int *tagsX, int *tagsY)
{
	int I = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(I >= QSize) return;
	
	QVals[I] *= ( tagsX[QRows[I]] == -1 && tagsY[QRows[I]] == -1 );
}

void generateQT(int *QTRows, int *QTCols, real *QTVals, int nx, int ny)
{
	int  numU = (nx-1)*ny;
	
	int Iu, Iv;
	int row = 0;
	int num_elements = 0;
	
	/// QT is an (np + 2*nb) x nuv matrix
	
	/// Generate the GT part
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx; i++)
		{
			Iu = j*(nx-1) + i;
			Iv = j*nx + i + numU;
			
			if(i>0)
			{
				QTRows[num_elements] = row;
				QTCols[num_elements] = Iu - 1;
				QTVals[num_elements] = 1;
				num_elements++;
			}
			if(i<nx-1)
			{
				QTRows[num_elements] = row;
				QTCols[num_elements] = Iu;
				QTVals[num_elements] = -1;
				num_elements++;
			}
			if(j>0)
			{
				QTRows[num_elements] = row;
				QTCols[num_elements] = Iv - nx;
				QTVals[num_elements] = 1;
				num_elements++;
			}
			if(j<ny-1)
			{
				QTRows[num_elements] = row;
				QTCols[num_elements] = Iv;
				QTVals[num_elements] = -1;
				num_elements++;
			}
			row++;
		}
	}
}



} // end of namespace kernels
