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

#define BLOCKSIZE 256

template <>
void TairaColoniusSolver<host_memory>::updateQT()
{
}

template <>
void TairaColoniusSolver<device_memory>::updateQT()
{
	logger.startTimer("updateQT");
	
	int  *QTRows = thrust::raw_pointer_cast(&(QT.row_indices[0])),
	     *QTCols = thrust::raw_pointer_cast(&(QT.column_indices[0]));
	real *QTVals = thrust::raw_pointer_cast(&(QT.values[0]));
	
	int  *ERows = thrust::raw_pointer_cast(&(E.row_indices[0])),
	     *ECols = thrust::raw_pointer_cast(&(E.column_indices[0]));
	real *EVals = thrust::raw_pointer_cast(&(E.values[0]));
	
	real *x  = thrust::raw_pointer_cast(&(domInfo->xD[0])),
	     *y  = thrust::raw_pointer_cast(&(domInfo->yD[0])),
	     *dx = thrust::raw_pointer_cast(&(domInfo->dxD[0]));
	
	real *xB = thrust::raw_pointer_cast(&(B.x[0])),
	     *yB = thrust::raw_pointer_cast(&(B.y[0]));
	
	int *I = thrust::raw_pointer_cast(&(B.I[0])),
	    *J = thrust::raw_pointer_cast(&(B.J[0]));
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	     
	dim3 dimGrid(int((B.totalPoints-0.5)/BLOCKSIZE)+1, 1);
	dim3 dimBlock(BLOCKSIZE, 1);

	kernels::updateQT <<<dimGrid, dimBlock>>>
	         (QTRows, QTCols, QTVals,
	          ERows,  ECols,  EVals,
	          nx, ny, x, y, dx,
              B.totalPoints, xB, yB, I, J);
    
    logger.stopTimer("updateQT");
    
	logger.startTimer("transposeQT");
	cusp::transpose(QT, Q);
	logger.stopTimer("transposeQT");
	
	logger.startTimer("transposeE");
	cusp::transpose(E, ET);
	logger.stopTimer("transposeE");
}

template <>
void TairaColoniusSolver<host_memory>::generateQT()
{
	logger.startTimer("generateQT");
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  numP  = numU + ny;
	
	/// QT is an (np + 2*nb) x nuv matrix
	QT.resize(numP + 2*B.totalPoints, numUV, 4*numP-2*(nx+ny) + 24*B.totalPoints);
	
	int  *QTRows = thrust::raw_pointer_cast(&(QT.row_indices[0])),
	     *QTCols = thrust::raw_pointer_cast(&(QT.column_indices[0]));
	real *QTVals = thrust::raw_pointer_cast(&(QT.values[0]));
	
	kernels::generateQT(QTRows, QTCols, QTVals, nx, ny);

	logger.stopTimer("generateQT");
	
	updateQT();	
}

template <>
void TairaColoniusSolver<device_memory>::generateQT()
{
	logger.startTimer("generateQT");
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  numP  = numU + ny;
	
	/// QT is an (np + 2*nb) x nuv matrix
	cooH QTHost(numP + 2*B.totalPoints, numUV, 4*numP-2*(nx+ny) + 24*B.totalPoints);
	
	int  *QTRows = thrust::raw_pointer_cast(&(QTHost.row_indices[0])),
	     *QTCols = thrust::raw_pointer_cast(&(QTHost.column_indices[0]));
	real *QTVals = thrust::raw_pointer_cast(&(QTHost.values[0]));
	
	kernels::generateQT(QTRows, QTCols, QTVals, nx, ny);

	QT = QTHost;
	
	logger.stopTimer("generateQT");
	
	updateQT();
}

#if 0
template <>
void TairaColoniusSolver<device_memory>::generateQT()
{
	logger.startTimer("generateQT");
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  numP  = numU + ny;
	
	/// QT is an (np + 2*nb) x nuv matrix
	cooH QTHost(numP + 2*B.totalPoints, numUV, 4*numP-2*(nx+ny) + 24*B.totalPoints);
	
	int  *QTRows = thrust::raw_pointer_cast(&(QTHost.row_indices[0])),
	     *QTCols = thrust::raw_pointer_cast(&(QTHost.column_indices[0]));
	real *QTVals = thrust::raw_pointer_cast(&(QTHost.values[0]));
	
	kernels::generateQT(QTRows, QTCols, QTVals, nx, ny);

	// Temparawari
	int row = numP;
	int num_elements = 4*numP-2*(nx+ny);
	
/*	
	/// Generate E
	int  I, J, col;
	real xB, yB, x, y, alpha, Dx;
	int cur = 0;
	for(int k=0; k<B.numBodies; k++)
	{
		Dx = B.h[k];
		alpha = Dx*Dx;
		for(int l=0; l<B.numPoints[k]; l++)
		{	
			xB = B.x[cur];
			yB = B.y[cur];
			I  = B.I[cur];
			J  = B.J[cur];
			for(int j=J-1; j<=J+1; j++)
			{
				for(int i=I-2; i<=I+1; i++)
				{
					col = j*(nx-1) + i;
					x   = domInfo->x[i+1];
					y   = 0.5*(domInfo->y[j] + domInfo->y[j+1]);
					QTHost.row_indices[num_elements] = row;
					QTHost.column_indices[num_elements] = col;
					QTHost.values[num_elements] = alpha*delta(x-xB, y-yB, Dx)/Dx; // Dx = Dy near the body
					num_elements++;
				}
			}
			cur++;
			row++;
		}
	}
	cur=0;
	for(int k=0; k<B.numBodies; k++)
	{
		Dx = B.h[k];
		alpha = Dx*Dx;
		for(int l=0; l<B.numPoints[k]; l++)
		{
			xB	= B.x[cur];
			yB	= B.y[cur];
			I	= B.I[cur];
			J	= B.J[cur];
			for(int j=J-2; j<=J+1; j++)
			{
				for(int i=I-1; i<=I+1; i++)
				{
					col	= j*nx + i + numU;
					x	= 0.5*(domInfo->x[i] + domInfo->x[i+1]);
					y	= domInfo->y[j+1];
					QTHost.row_indices[num_elements] = row;
					QTHost.column_indices[num_elements] = col;
					QTHost.values[num_elements] = alpha*delta(x-xB, y-yB, Dx)/Dx;
					num_elements++;
				}
			}
			cur++;
			row++;
		}
	}
*/	
	/// Generate the E part
/*	for(int k=0; k<B.totalPoints; k++)
	{
		Dx = NavierStokesSolver<Matrix, Vector>::domInfo->dx[B.I[k]];
		alpha = Dx*Dx;
		xB = B.x[k];
		yB = B.y[k];
		I  = B.I[k];
		J  = B.J[k];
		for(int j=J-1; j<=J+1; j++)
		{
			for(i=I-2; i<=I+1; i++)
			{
				col	= j*(nx-1) + i;
				x	= NavierStokesSolver<Matrix, Vector>::domInfo->x[i+1];
				y	= 0.5*(NavierStokesSolver<Matrix, Vector>::domInfo->y[j] + NavierStokesSolver<Matrix, Vector>::domInfo->y[j+1]);
				QTHost.row_indices[num_elements] = row;
				QTHost.column_indices[num_elements] = col;
				QTHost.values[num_elements] = alpha*delta(x-xB, y-yB, Dx)/Dx; // Dx = Dy near the body
				num_elements++;
			}
		}
		row++;
	}
	for(int k=0; k<B.totalPoints; k++)
	{
		Dx = NavierStokesSolver<Matrix, Vector>::domInfo->dx[B.I[k]];
		alpha = Dx*Dx;
		xB	= B.x[k];
		yB	= B.y[k];
		I	= B.I[k];
		J	= B.J[k];
		for(j=J-2; j<=J+1; j++)
		{
			for(i=I-1; i<=I+1; i++)
			{
				col	= j*nx + i + numU;
				x	= 0.5*(NavierStokesSolver<Matrix, Vector>::domInfo->x[i] + NavierStokesSolver<Matrix, Vector>::domInfo->x[i+1]);
				y	= NavierStokesSolver<Matrix, Vector>::domInfo->y[j+1];
				QTHost.row_indices[num_elements] = row;
				QTHost.column_indices[num_elements] = col;
				QTHost.values[num_elements] = alpha*delta(x-xB, y-yB, Dx)/Dx;
				num_elements++;
			}
		}
		row++;
	}
*/
	QT = QTHost;
	
	logger.stopTimer("generateQT");
	
	logger.startTimer("transposeQT");
	cusp::transpose(QT, Q);
	logger.stopTimer("transposeQT");
}
#endif
