template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::generateM()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  I;
	real value;
	
	M.resize(numUV, numUV, numUV);
	Minv.resize(numUV, numUV, numUV);
	
	for (int j=0; j < ny; j++)
	{
		for (int i=0; i < nx-1; i++)
		{
			I = j*(nx-1) + i;
			value = 0.5*(domInfo->dx[i+1]+domInfo->dx[i]) / domInfo->dy[j] / simPar->dt;
			
			M.row_indices[I] = I;
			M.column_indices[I] = I;
			M.values[I] = value;
			
			Minv.row_indices[I] = I;
			Minv.column_indices[I] = I;
			Minv.values[I] = 1.0/value;
		}
	}
	for (int j=0; j < ny-1; j++)
	{
		for (int i=0; i < nx; i++)
		{
			I = j*nx + i + numU;
			
			value  = 0.5*(domInfo->dy[j+1]+domInfo->dy[j]) / domInfo->dx[i] / simPar->dt;
			
			M.row_indices[I] = I;
			M.column_indices[I] = I;
			M.values[I] = value;
			
			Minv.row_indices[I] = I;
			Minv.column_indices[I] = I;
			Minv.values[I] = 1.0/value;
		}
	}
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::generateL()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	real Cx0 = 0.0, Cx1 = 0.0, Cy0 = 0.0, Cy1 = 0.0,
	     scale = 0.0,
	     dx0 = 1.0, dx1 = 1.0, dy0 = 1.0, dy1 = 1.0;
	int  numU = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  I,
	     num_elements = 0;
			
	L.resize(numUV, numUV, 5*numUV-4*(nx+ny)+4);
	
	///x-component
	for (int j=0; j < ny; j++)
	{
		if(j == 0)
		{
			dy0 = 0.5*domInfo->dy[0];
			dy1 = 0.5*(domInfo->dy[0]+domInfo->dy[1]);
		}
		else if(j == ny-1)
		{
			dy0 = 0.5*(domInfo->dy[ny-2]+domInfo->dy[ny-1]);
			dy1 = 0.5*domInfo->dy[ny-1];
		}
		else
		{
			dy0 = 0.5*(domInfo->dy[j]+domInfo->dy[j-1]);
			dy1 = 0.5*(domInfo->dy[j]+domInfo->dy[j+1]);
		}
		Cy0 = 2.0 * flowDesc->nu / ( dy1*(dy1+dy0) );
		Cy1 = 2.0 * flowDesc->nu / ( dy0*(dy1+dy0) );
		
		for (int i=0; i < nx-1; i++)
		{
			I = j*(nx-1) + i;										///< calculate the row of the matrix
			
			Cx0 = 2.0 * flowDesc->nu / ( domInfo->dx[i+1]*(domInfo->dx[i+1]+domInfo->dx[i]) );
			Cx1 = 2.0 * flowDesc->nu / ( domInfo->dx[i]*(domInfo->dx[i+1]+domInfo->dx[i]) );
			
			scale = 0.5*(domInfo->dx[i+1]+domInfo->dx[i]);					///< scaling factor (to obtain the normalised matrix)
			
			// south
			if(j>0)													///< no south coefficient for the bottom row
			{
				L.row_indices[num_elements] = I;
				L.column_indices[num_elements] = I - (nx-1);
				L.values[num_elements] = Cy1 * scale / domInfo->dy[j-1];
				num_elements++;
			}
			// west
			if(i>0)													///< no west coefficient for the leftmost column			
			{
				L.row_indices[num_elements] = I;
				L.column_indices[num_elements] = I - 1;
				L.values[num_elements] = Cx1 * scale / domInfo->dy[j];
				num_elements++;
			}
			// centre
			L.row_indices[num_elements] = I;
			L.column_indices[num_elements] = I;
			L.values[num_elements] = (-Cy0-Cy1-Cx0-Cx1) * scale / domInfo->dy[j];
			num_elements++;
			// east
			if(i<nx-2)											// no east coefficient for the rightmost column
			{
				L.row_indices[num_elements] = I;
				L.column_indices[num_elements] = I + 1;
				L.values[num_elements] = Cx0 * scale / domInfo->dy[j];
				num_elements++;
			}
			// north
			if(j<ny-1)
			{
				L.row_indices[num_elements] = I;
				L.column_indices[num_elements] = I + (nx-1);
				L.values[num_elements] = Cy0 * scale / domInfo->dy[j+1];
				num_elements++;
			}
		}
	}
	
	///y-component
	for (int j=0; j < ny-1; j++)
	{
		Cy0 = 2.0 * flowDesc->nu / ( domInfo->dy[j+1]*(domInfo->dy[j+1]+domInfo->dy[j]) );
		Cy1 = 2.0 * flowDesc->nu / ( domInfo->dy[j]*(domInfo->dy[j+1]+domInfo->dy[j]) );
		
		for (int i=0; i < nx; i++)
		{
			I = j*nx + i + numU;		// calculate the row of the matrix
			
			if((i>0) && (i<nx-1))
			{
				dx0 = 0.5*(domInfo->dx[i]+domInfo->dx[i-1]);
				dx1 = 0.5*(domInfo->dx[i]+domInfo->dx[i+1]);
			}
			else if(i==0)
			{
				dx0 = 0.5*domInfo->dx[i];
				dx1 = 0.5*(domInfo->dx[i]+domInfo->dx[i+1]);
			}
			else if(i==nx-1)
			{
				dx0 = 0.5*(domInfo->dx[i]+domInfo->dx[i-1]);
				dx1 = 0.5*domInfo->dx[i];
			}
			Cx0 = 2.0 * flowDesc->nu /( dx1*(dx1+dx0) );
			Cx1 = 2.0 * flowDesc->nu /( dx0*(dx1+dx0) );
			
			scale = 0.5*(domInfo->dy[j+1]+domInfo->dy[j]);	// scaling factor (to obtain the normalised matrix)
			
			// south
			if(j>0)													// no south coefficient for the bottom row
			{
				L.row_indices[num_elements] = I;
				L.column_indices[num_elements] = I - nx;
				L.values[num_elements] = Cy1 * scale / domInfo->dx[i];
				num_elements++;
			}
			// west
			if(i>0)													// no west coefficient for the leftmost column			
			{
				L.row_indices[num_elements] = I;
				L.column_indices[num_elements] = I - 1;
				L.values[num_elements] = Cx1 * scale / domInfo->dx[i-1];
				num_elements++;
			}
			// centre
			L.row_indices[num_elements] = I;
			L.column_indices[num_elements] = I;
			L.values[num_elements] = (-Cy0-Cy1-Cx0-Cx1) * scale / domInfo->dx[i];
			num_elements++;
			// east
			if(i<nx-1)												// no east coefficient for the rightmost column
			{
				L.row_indices[num_elements] = I;
				L.column_indices[num_elements] = I + 1;
				L.values[num_elements] = Cx0 * scale / domInfo->dx[i+1];
				num_elements++;
			}
			// north
			if(j<ny-2)
			{
				L.row_indices[num_elements] = I;
				L.column_indices[num_elements] = I + nx;
				L.values[num_elements] = Cy0 * scale / domInfo->dx[i];
				num_elements++;
			}
		}
	}
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::generateA()
{
	printf("entered generateA\n");
	generateM();
	printf("finished generateM\n");
	generateL();
printf("finished generateL\n");
	cusp::wrapped::subtract(M, L, A);
}

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::generateBN()
{
	BN = Minv;
}
/*
template <typename Matrix, typename Vector>
template <>
void NavierStokesSolver<Matrix, Vector>::generateBN<3>()
{
	Matrix	temp1, temp2;
	cusp::multiply(Minv, L, temp1);
	cusp::multiply(temp1, Minv, BN);
	cusp::add(Minv, BN, BN);
	cusp::multiply(temp1, BN, temp2);
	cusp::add(Minv, temp2, BN);
}*/