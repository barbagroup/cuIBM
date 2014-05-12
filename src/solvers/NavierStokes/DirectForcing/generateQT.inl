#include <solvers/NavierStokes/kernels/generateQT.h>

template <>
void DirectForcingSolver<host_memory>::updateQ()
{
}

template <>
void DirectForcingSolver<device_memory>::updateQ()
{
	const int blocksize = 256;
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  QSize = 4*nx*ny-2*(nx+ny);
	
	int  *QRows = thrust::raw_pointer_cast(&(Q.row_indices[0])),
	     *QCols = thrust::raw_pointer_cast(&(Q.column_indices[0]));
//	int  *tags_r = thrust::raw_pointer_cast(&(tagsD[0]));
	int  *tagsX_r = thrust::raw_pointer_cast(&(tagsXD[0])),\
	     *tagsY_r = thrust::raw_pointer_cast(&(tagsYD[0]));

	real *QVals = thrust::raw_pointer_cast(&(Q.values[0]));
	
	dim3 dimGrid( int((QSize-0.5)/blocksize) + 1, 1);
	dim3 dimBlock(blocksize, 1);
	
//	kernels::updateQ <<<dimGrid, dimBlock>>> (QRows, QCols, QVals, QSize, tags_r);
	kernels::updateQ <<<dimGrid, dimBlock>>> (QRows, QCols, QVals, QSize, tagsX_r, tagsY_r);
}

template <typename memoryType>
void DirectForcingSolver<memoryType>::generateQT()
{
	NavierStokesSolver<memoryType>::generateQT();
	updateQ();
	cusp::transpose(NavierStokesSolver<memoryType>::Q, NavierStokesSolver<memoryType>::QT);
}

/*
template <>
void DirectForcingSolver<host_memory>::generateQT()
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
				QT.values[num_elements] = (tags[Iu-1] != Iu); //1;
				num_elements++;
			}
			if(i<nx-1)
			{
				QT.row_indices[num_elements] = row;
				QT.column_indices[num_elements] = Iu;
				QT.values[num_elements] = -(tags[Iu] != Iu-1 || Iu == 0); //-1;
				// the Iu==0 condition is there to take care of the special case
				num_elements++;
			}
			if(j>0)
			{
				QT.row_indices[num_elements] = row;
				QT.column_indices[num_elements] = Iv - nx;
				QT.values[num_elements] = (tags[Iv-nx] != Iv);//1;
				num_elements++;
			}
			if(j<ny-1)
			{
				QT.row_indices[num_elements] = row;
				QT.column_indices[num_elements] = Iv;
				QT.values[num_elements] = -(tags[Iv] != Iv-nx);//-1;
				num_elements++;
			}
			row++;
		}
	}
	
	cusp::transpose(QT, Q);
	
	updateQDirectForcing();
}

template <>
void DirectForcingSolver<device_memory>::generateQT()
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
				QTHost.values[num_elements] = (tags[Iu-1] != Iu); //1;
				num_elements++;
			}
			if(i<nx-1)
			{
				QTHost.row_indices[num_elements] = row;
				QTHost.column_indices[num_elements] = Iu;
				QTHost.values[num_elements] = -(tags[Iu] != Iu-1 || Iu == 0); //-1;
				// the Iu==0 condition is there to take care of the special case
				num_elements++;
			}
			if(j>0)
			{
				QTHost.row_indices[num_elements] = row;
				QTHost.column_indices[num_elements] = Iv - nx;
				QTHost.values[num_elements] = (tags[Iv-nx] != Iv);//1;
				num_elements++;
			}
			if(j<ny-1)
			{
				QTHost.row_indices[num_elements] = row;
				QTHost.column_indices[num_elements] = Iv;
				QTHost.values[num_elements] = -(tags[Iv] != Iv-nx);//-1;
				num_elements++;
			}
			row++;
		}
	}
	
	QT = QTHost;
//	cusp::print(QT);
	cusp::transpose(QT, Q);
	
	updateQDirectForcing();
	
	cusp::transpose(Q, QT);
//	cusp::print(Q);
}*/


