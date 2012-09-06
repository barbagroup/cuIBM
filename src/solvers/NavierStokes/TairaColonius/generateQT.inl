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

real dh_roma(real x, real h)
{
	real r = fabs(x)/h;
	
	if(r>1.5)
		return 0.0;
	else if(r>0.5 && r<=1.5)
		return 1.0/(6*h)*( 5.0 - 3.0*r - sqrt(-3.0*(1-r)*(1-r) + 1.0) );
	else
		return 1.0/(3*h)*( 1.0 + sqrt(-3.0*r*r + 1.0) );
}

real delta(real x, real y, real h)
{
	return dh_roma(x, h) * dh_roma(y, h);
}

template <>
void TairaColoniusSolver<host_memory>::generateQT()
{
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
	
	// Temparawari
	int row = numP;
	int num_elements = 4*numP-2*(nx+ny);
	
	/// Generate E
	int  I, J, col;
	real xB, yB, x, y, alpha, Dx;
	int cur = 0;
	for(int k=0; k<B.numBodies; k++)
	{
		Dx = domInfo->dx[B.I[B.offsets[k]]];
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
					QT.row_indices[num_elements] = row;
					QT.column_indices[num_elements] = col;
					QT.values[num_elements] = alpha*delta(x-xB, y-yB, Dx)/Dx; // Dx = Dy near the body
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
		Dx = domInfo->dx[B.I[B.offsets[k]]];
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
					QT.row_indices[num_elements] = row;
					QT.column_indices[num_elements] = col;
					QT.values[num_elements] = alpha*delta(x-xB, y-yB, Dx)/Dx;
					num_elements++;
				}
			}
			cur++;
			row++;
		}
	}
	
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
				NavierStokesSolver<Matrix, Vector>::QT.row_indices[num_elements] = row;
				NavierStokesSolver<Matrix, Vector>::QT.column_indices[num_elements] = col;
				NavierStokesSolver<Matrix, Vector>::QT.values[num_elements] = alpha*delta(x-xB, y-yB, Dx)/Dx; // Dx = Dy near the body
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
				NavierStokesSolver<Matrix, Vector>::QT.row_indices[num_elements] = row;
				NavierStokesSolver<Matrix, Vector>::QT.column_indices[num_elements] = col;
				NavierStokesSolver<Matrix, Vector>::QT.values[num_elements] = alpha*delta(x-xB, y-yB, Dx)/Dx;
				num_elements++;
			}
		}
		row++;
	}
*/
	cusp::transpose(QT, Q);
}

template <>
void TairaColoniusSolver<device_memory>::generateQT()
{
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

/*
	int Iu, Iv;
	int row = 0;
	int num_elements = 0;
	
	/// Generate the GT part
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx; i++)
		{
			Iu = j*(nx-1) + i;
			Iv = j*nx + i + numU;
			
			if(i>0)
			{
				QTHost.row_indices[num_elements] = row;
				QTHost.column_indices[num_elements] = Iu - 1;
				QTHost.values[num_elements] = 1;
				num_elements++;
			}
			if(i<nx-1)
			{
				QTHost.row_indices[num_elements] = row;
				QTHost.column_indices[num_elements] = Iu;
				QTHost.values[num_elements] = -1;
				num_elements++;
			}
			if(j>0)
			{
				QTHost.row_indices[num_elements] = row;
				QTHost.column_indices[num_elements] = Iv - nx;
				QTHost.values[num_elements] = 1;
				num_elements++;
			}
			if(j<ny-1)
			{
				QTHost.row_indices[num_elements] = row;
				QTHost.column_indices[num_elements] = Iv;
				QTHost.values[num_elements] = -1;
				num_elements++;
			}
			row++;
		}
	}
*/
	// Temparawari
	int row = numP;
	int num_elements = 4*numP-2*(nx+ny);
	
	/// Generate E
	int  I, J, col;
	real xB, yB, x, y, alpha, Dx;
	int cur = 0;
	for(int k=0; k<B.numBodies; k++)
	{
		Dx = domInfo->dx[B.I[B.offsets[k]]];
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
		Dx = domInfo->dx[B.I[B.offsets[k]]];
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
	cusp::transpose(QT, Q);
	//cusp::print(QT);
	std::cout << "Generated QT!" << std::endl;
}

template <typename memoryType>
void TairaColoniusSolver<memoryType>::updateQT()
{	
}

template <>
void TairaColoniusSolver<host_memory>::generateE()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	
	/// QT is an (np + 2*nb) x nuv matrix
	E.resize(2*B.totalPoints, numUV, 24*B.totalPoints);
	
	int row = 0;
	int num_elements = 0;

	/// Generate E
	int  I, J, col;
	real xB, yB, x, y, alpha, Dx;
	int cur = 0;
	for(int k=0; k<B.numBodies; k++)
	{
		Dx = domInfo->dx[B.I[B.offsets[k]]];
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
					E.row_indices[num_elements] = row;
					E.column_indices[num_elements] = col;
					E.values[num_elements] = alpha*delta(x-xB, y-yB, Dx)/Dx; // Dx = Dy near the body
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
		Dx = domInfo->dx[B.I[B.offsets[k]]];
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
					E.row_indices[num_elements] = row;
					E.column_indices[num_elements] = col;
					E.values[num_elements] = alpha*delta(x-xB, y-yB, Dx)/Dx;
					num_elements++;
				}
			}
			cur++;
			row++;
		}
	}
	cusp::transpose(E, ET);
}

template <>
void TairaColoniusSolver<device_memory>::generateE()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	
	/// QT is an (np + 2*nb) x nuv matrix
	cooH EHost(2*B.totalPoints, numUV, 24*B.totalPoints);
	
	int row = 0;
	int num_elements = 0;

	/// Generate E
	int  I, J, col;
	real xB, yB, x, y, alpha, Dx;
	int cur = 0;
	for(int k=0; k<B.numBodies; k++)
	{
		Dx = domInfo->dx[B.I[B.offsets[k]]];
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
					EHost.row_indices[num_elements] = row;
					EHost.column_indices[num_elements] = col;
					EHost.values[num_elements] = alpha*delta(x-xB, y-yB, Dx)/Dx; // Dx = Dy near the body
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
		Dx = domInfo->dx[B.I[B.offsets[k]]];
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
					EHost.row_indices[num_elements] = row;
					EHost.column_indices[num_elements] = col;
					EHost.values[num_elements] = alpha*delta(x-xB, y-yB, Dx)/Dx;
					num_elements++;
				}
			}
			cur++;
			row++;
		}
	}
	E = EHost;
	cusp::transpose(E, ET);
}
