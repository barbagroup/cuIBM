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