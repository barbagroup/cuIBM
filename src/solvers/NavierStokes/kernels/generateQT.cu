/***************************************************************************//**
* \file generateQT.cu
* \author Krishnan, A. (anush@bu.edu)
* \brief Definition of the kernels required to generate the matrix \c QT
*/

#include <solvers/NavierStokes/kernels/generateQT.h>

/**
* \brief Calculate the discrete Delta function from Rome et al. (1999)
*
* \param x x- or y- component of the vector defined between two points
* \param h grid-spacing
*
* \return the value of the discrete delta function 
*/
__device__ \
real dhRomaDeviceQT(real x, real h)
{
	real r = fabs(x)/h;
	
	if(r>1.5)
		return 0.0;
	else if(r>0.5 && r<=1.5)
		return 1.0/(6*h)*( 5.0 - 3.0*r - sqrt(-3.0*(1-r)*(1-r) + 1.0) );
	else
		return 1.0/(3*h)*( 1.0 + sqrt(-3.0*r*r + 1.0) );
}

/**
* \brief Calculate the discrete delta function in 2D
*
* \param x x-component of the vector defined between two points
* \param y y-component of the vector defined between two points
* \param h grid-spacing
*
* \return the value of the discrete Delta function in 2D
*/
__device__ \
real deltaDeviceQT(real x, real y, real h)
{
	return dhRomaDeviceQT(x, h) * dhRomaDeviceQT(y, h);
}

/**
* \namespace kernels
* \brief Contain all custom-written CUDA kernels
*/
namespace kernels
{

/**
* \brief To be documented
*/
__global__ \
void updateQ(int *QRows, int *QCols, real *QVals, int QSize, int *tags)
{
	int I = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(I >= QSize) return;
	
	QVals[I] *= ( tags[QRows[I]] == -1 );
}

/**
* \brief To be documented
*/
__global__ \
void updateQ(int *QRows, int *QCols, real *QVals, int QSize, int *tagsX, int *tagsY)
{
	int I = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(I >= QSize) return;
	
	QVals[I] *= ( tagsX[QRows[I]] == -1 && tagsY[QRows[I]] == -1 );
}

/**
* \brief Generate the matrix \c QT.
*
* It actually generates the GT part (discrete divergence operator).
*
* \param QTRows row index of elements of the COO-format divergence matrix \c QT
* \param QTCols column index of elements of the COO-format divergence matrix \c QT
* \param QTVals value of elements of the COO-format divergence matrix \c QT
* \param nx number of cells in the x-direction
* \param ny number of cells in the y-direction
*/
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

/**
* \brief Kernel to update element of the matrix \c QT on the device.
*
* It also updates element of the interpolation matrix \c E.
*
* \param QTRows row index of elements of the COO-format matrix \c QT
* \param QTCols column index of elements of the COO-format matrix \c QT
* \param QTVals value of elements of the COO-format matrix \c QT
* \param ERows row index of elements of the COO-format interpolation matrix \c E
* \param ECols column index of elements of the COO-format interpolation matrix \c E
* \param EVals value of elements of the COO-format interpolation matrix \c E
* \param nx number of cells in the x-direction
* \param ny number of cells in the y-direction
* \param x x-component of grid points
* \param y y-component of grid points
* \param dx cell widths in the x-direction
* \param totalPoints number of body points (all bodies included)
* \param xB x-component of body points (all bodies included)
* \param yB y-component of body points (all bodies included)
* \param I x-index of grid cells in which body points are located
* \param J y-index of grid cells in which body points are located
*/
__global__ \
void updateQT(int *QTRows, int *QTCols, real *QTVals,
              int *ERows,  int *ECols,  real *EVals,
              int nx, int ny, real *x, real *y, real *dx,
              int totalPoints, real *xB, real *yB, int *I, int *J)
{
	int bodyIdx = threadIdx.x + blockIdx.x*blockDim.x;
	
	if(bodyIdx >= totalPoints)
		return;
	
	int  Ib=I[bodyIdx],
	     Jb=J[bodyIdx],
	     QTIdx = 4*nx*ny - 2*(nx+ny) + bodyIdx*12,
	     EIdx  = bodyIdx*12,
	     i, j;

	real Dx = dx[Ib];
	
	// populate x-components
	for(j=Jb-1; j<=Jb+1; j++)
	{
		for(i=Ib-2; i<=Ib+1; i++)
		{
			QTRows[QTIdx] = bodyIdx + nx*ny;
			ERows[EIdx] = bodyIdx;
			
			QTCols[QTIdx] = j*(nx-1) + i;
			ECols[EIdx] = QTCols[QTIdx];
			
			QTVals[QTIdx] = Dx*deltaDeviceQT(x[i+1]-xB[bodyIdx], 0.5*(y[j]+y[j+1])-yB[bodyIdx], Dx);
			EVals[EIdx] = QTVals[QTIdx];
			
			QTIdx++;
			EIdx++;
		}
	}

	// populate y-components
	for(j=Jb-2; j<=Jb+1; j++)
	{
		for(i=Ib-1; i<=Ib+1; i++)
		{
			QTRows[QTIdx+12*totalPoints-12] = bodyIdx + nx*ny + totalPoints;
			ERows[EIdx+12*totalPoints-12] = bodyIdx + totalPoints;
			
			QTCols[QTIdx+12*totalPoints-12] = j*nx + i + (nx-1)*ny;
			ECols[EIdx+12*totalPoints-12] = QTCols[QTIdx+12*totalPoints-12];
			
			QTVals[QTIdx+12*totalPoints-12] = Dx*deltaDeviceQT(0.5*(x[i]+x[i+1])-xB[bodyIdx], y[j+1]-yB[bodyIdx], Dx);
			EVals[EIdx+12*totalPoints-12] = QTVals[QTIdx+12*totalPoints-12];
			
			QTIdx++;
			EIdx++;
		}
	}
}

/**
* \brief Update the matrix \c QT on the host.
*
* It also updates the interpolation matrix \c E.
*
* \param QTRows row index of elements of the COO-format matrix \c QT
* \param QTCols column index of elements of the COO-format matrix \c QT
* \param QTVals value of elements of the COO-format matrix \c QT
* \param ERows row index of elements of the COO-format interpolation matrix \c E
* \param ECols column index of elements of the COO-format interpolation matrix \c E
* \param EVals value of elements of the COO-format interpolation matrix \c E
* \param nx number of cells in the x-direction
* \param ny number of cells in the y-direction
* \param x x-component of grid points
* \param y y-component of grid points
* \param dx cell widths in the x-direction
* \param totalPoints number of body points (all bodies included)
* \param xB x-component of body points (all bodies included)
* \param yB y-component of body points (all bodies included)
* \param I x-index of grid cells in which body points are located
* \param J y-index of grid cells in which body points are located
*/
void updateQTHost(int *QTRows, int *QTCols, real *QTVals,
              int *ERows,  int *ECols,  real *EVals,
              int nx, int ny, real *x, real *y, real *dx,
              int totalPoints, real *xB, real *yB, int *I, int *J)
{
	for(int bodyIdx=0; bodyIdx<totalPoints; bodyIdx++)
	{
		int  Ib=I[bodyIdx],
		     Jb=J[bodyIdx],
		     QTIdx = 4*nx*ny - 2*(nx+ny) + bodyIdx*12,
		     EIdx  = bodyIdx*12,
		     i, j;

		real Dx = dx[Ib];
	
		// populate x-components
		for(j=Jb-1; j<=Jb+1; j++)
		{
			for(i=Ib-2; i<=Ib+1; i++)
			{
				QTRows[QTIdx] = bodyIdx + nx*ny;
				ERows[EIdx] = bodyIdx;
			
				QTCols[QTIdx] = j*(nx-1) + i;
				ECols[EIdx] = QTCols[QTIdx];
			
				QTVals[QTIdx] = Dx*delta(x[i+1]-xB[bodyIdx], 0.5*(y[j]+y[j+1])-yB[bodyIdx], Dx);
				EVals[EIdx] = QTVals[QTIdx];
			
				QTIdx++;
				EIdx++;
			}
		}

		// populate y-components
		for(j=Jb-2; j<=Jb+1; j++)
		{
			for(i=Ib-1; i<=Ib+1; i++)
			{
				QTRows[QTIdx+12*totalPoints-12] = bodyIdx + nx*ny + totalPoints;
				ERows[EIdx+12*totalPoints-12] = bodyIdx + totalPoints;
			
				QTCols[QTIdx+12*totalPoints-12] = j*nx + i + (nx-1)*ny;
				ECols[EIdx+12*totalPoints-12] = QTCols[QTIdx+12*totalPoints-12];
			
				QTVals[QTIdx+12*totalPoints-12] = Dx*delta(0.5*(x[i]+x[i+1])-xB[bodyIdx], y[j+1]-yB[bodyIdx], Dx);
				EVals[EIdx+12*totalPoints-12] = QTVals[QTIdx+12*totalPoints-12];
			
				QTIdx++;
				EIdx++;
			}
		}
	}
}

} // end of namespace kernels
