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
void TairaColoniusSolver<cooH, vecH>::generateQT()
{
}

template <>
void TairaColoniusSolver<cooD, vecD>::generateQT()
{
	std::cout << "Entered generateQT" << std::endl;
	int  nx = NavierStokesSolver<Matrix, Vector>::domInfo->nx,
	     ny = NavierStokesSolver<Matrix, Vector>::domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  numP  = numU + ny;
	
	/// QT is an (np + 2*nb) x nuv matrix
	NavierStokesSolver<Matrix, Vector>::QT.resize(numP + 2*B.totalPoints, numUV, 4*numP-2*(nx+ny) + 24*B.totalPoints);
	
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
				NavierStokesSolver<Matrix, Vector>::QT.row_indices[num_elements] = row;
				NavierStokesSolver<Matrix, Vector>::QT.column_indices[num_elements] = Iu - 1;
				NavierStokesSolver<Matrix, Vector>::QT.values[num_elements] = 1;
				num_elements++;
			}
			if(i<nx-1)
			{
				NavierStokesSolver<Matrix, Vector>::QT.row_indices[num_elements] = row;
				NavierStokesSolver<Matrix, Vector>::QT.column_indices[num_elements] = Iu;
				NavierStokesSolver<Matrix, Vector>::QT.values[num_elements] = -1;
				num_elements++;
			}
			if(j>0)
			{
				NavierStokesSolver<Matrix, Vector>::QT.row_indices[num_elements] = row;
				NavierStokesSolver<Matrix, Vector>::QT.column_indices[num_elements] = Iv - nx;
				NavierStokesSolver<Matrix, Vector>::QT.values[num_elements] = 1;
				num_elements++;
			}
			if(j<ny-1)
			{
				NavierStokesSolver<Matrix, Vector>::QT.row_indices[num_elements] = row;
				NavierStokesSolver<Matrix, Vector>::QT.column_indices[num_elements] = Iv;
				NavierStokesSolver<Matrix, Vector>::QT.values[num_elements] = -1;
				num_elements++;
			}
			row++;
		}
	}
	std::cout << "Generated GT!" << std::endl;
	/// Generate E
	int  I, J, col;
	real xB, yB, x, y, alpha, Dx;
	int cur = 0;
	for(int k=0; k<B.numBodies; k++)
	{
		Dx = NavierStokesSolver<Matrix, Vector>::domInfo->dx[B.I[B.offsets[k]]];
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
					x   = NavierStokesSolver<Matrix, Vector>::domInfo->x[i+1];
					y   = 0.5*(NavierStokesSolver<Matrix, Vector>::domInfo->y[j] + NavierStokesSolver<Matrix, Vector>::domInfo->y[j+1]);
					NavierStokesSolver<Matrix, Vector>::QT.row_indices[num_elements] = row;
					NavierStokesSolver<Matrix, Vector>::QT.column_indices[num_elements] = col;
					NavierStokesSolver<Matrix, Vector>::QT.values[num_elements] = alpha*delta(x-xB, y-yB, Dx)/Dx; // Dx = Dy near the body
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
		Dx = NavierStokesSolver<Matrix, Vector>::domInfo->dx[B.I[B.offsets[k]]];
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
					x	= 0.5*(NavierStokesSolver<Matrix, Vector>::domInfo->x[i] + NavierStokesSolver<Matrix, Vector>::domInfo->x[i+1]);
					y	= NavierStokesSolver<Matrix, Vector>::domInfo->y[j+1];
					NavierStokesSolver<Matrix, Vector>::QT.row_indices[num_elements] = row;
					NavierStokesSolver<Matrix, Vector>::QT.column_indices[num_elements] = col;
					NavierStokesSolver<Matrix, Vector>::QT.values[num_elements] = alpha*delta(x-xB, y-yB, Dx)/Dx;
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
	cusp::transpose(NavierStokesSolver<Matrix, Vector>::QT, NavierStokesSolver<Matrix, Vector>::Q);
}

template <typename Matrix, typename Vector>
void TairaColoniusSolver<Matrix, Vector>::updateQT()
{	
}