/**
*  Copyright (C) 2011 by Anush Krishnan, Simon Layton, Lorena Barba
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

#if 0
template <typename memoryType>
void NavierStokesSolver<memoryType>::generateQT(int *QTRows, int *QTCols, real *QTVals)
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;

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

template <>
void NavierStokesSolver<host_memory>::generateQT()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  numP  = numU + ny;
	
	QT.resize(numP, numUV, 4*numP-2*(nx+ny));
	
	int  *QTRows = thrust::raw_pointer_cast(&(QT.row_indices[0])),
	     *QTCols = thrust::raw_pointer_cast(&(QT.column_indices[0]));
	real *QTVals = thrust::raw_pointer_cast(&(QT.values[0]));
	
	generateQT(QTRows, QTCols, QTVals);
	cusp::transpose(QT, Q);
}

template <>
void NavierStokesSolver<device_memory>::generateQT()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  numP  = numU + ny;
	
	cooH QTHost(numP, numUV, 4*numP-2*(nx+ny));
	
	int  *QTRows = thrust::raw_pointer_cast(&(QTHost.row_indices[0])),
	     *QTCols = thrust::raw_pointer_cast(&(QTHost.column_indices[0]));
	real *QTVals = thrust::raw_pointer_cast(&(QTHost.values[0]));
	
	generateQT(QTRows, QTCols, QTVals);
	QT = QTHost;
	cusp::print(QT);
	cusp::transpose(QT, Q);
}
#endif

#if 1
template <>
void NavierStokesSolver<host_memory>::generateQT()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  numP  = numU + ny;
	QT.resize(numP, numUV, 4*numP-2*(nx+ny));
	
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
				QT.row_indices[num_elements] = row;
				QT.column_indices[num_elements] = Iu - 1;
				QT.values[num_elements] = 1;
				num_elements++;
			}
			if(i<nx-1)
			{
				QT.row_indices[num_elements] = row;
				QT.column_indices[num_elements] = Iu;
				QT.values[num_elements] = -1;
				num_elements++;
			}
			if(j>0)
			{
				QT.row_indices[num_elements] = row;
				QT.column_indices[num_elements] = Iv - nx;
				QT.values[num_elements] = 1;
				num_elements++;
			}
			if(j<ny-1)
			{
				QT.row_indices[num_elements] = row;
				QT.column_indices[num_elements] = Iv;
				QT.values[num_elements] = -1;
				num_elements++;
			}
			row++;
		}
	}
	
	cusp::transpose(QT, Q);
}

template <>
void NavierStokesSolver<device_memory>::generateQT()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  numP  = numU + ny;
	
	cooH QTHost(numP, numUV, 4*numP-2*(nx+ny));
	
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
	
	QT = QTHost;
	cusp::transpose(QT, Q);
}
#endif
