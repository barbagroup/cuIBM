template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::generateQT()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	int  numU  = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  numP  = numU + ny;
//	int  numB  = 0;	
//	for(int k=0; k<flowDesc->numBodies; k++)
//		numB += flowDesc->B[k].numPoints;
	
//	QT.resize(numP + 2*numB, numUV, 4*numP-2*(nx+ny) + 24*numB)
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

/*	int			I, J, col, k, l;
	PRECISION	xB, yB, x, y, alpha, Dx;
	
	/// Generate the E part
	for(int k=0; k<B.no_of_bodies; k++)
	{
		Dx = D.dx[B.Body[k].bp[0].I];
		alpha = Dx*Dx;
		for(int l=0; l<B.Body[k].no_of_points; l++)
		{
			xB	= B.Body[k].bp[l].x;
			yB	= B.Body[k].bp[l].y;
			I	= B.Body[k].bp[l].I;
			J	= B.Body[k].bp[l].J;
			
			for(int j=J-1; j<=J+1; j++)
			{
				for(i=I-2; i<=I+1; i++)
				{
					col	= j*(nx-1) + i;
					x	= D.x[i+1];
					y	= 0.5*(D.y[j] + D.y[j+1]);
					QT.row_indices[num_elements] = row;
					QT.column_indices[num_elements] = col;
					QT.values[num_elements] = alpha*delta(x-xB, y-yB, Dx)/Dx; // Dx = Dy near the body
					num_elements++;
				}
			}		
			row++;
		}
	}
	
	for(int k=0; k<B.no_of_bodies; k++)
	{
		Dx = D.dx[B.Body[k].bp[0].I];
		alpha = Dx*Dx;
		for(int l=0; l<B.Body[k].no_of_points; l++)
		{
			xB	= B.Body[k].bp[l].x;
			yB	= B.Body[k].bp[l].y;
			I	= B.Body[k].bp[l].I;
			J	= B.Body[k].bp[l].J;
			
			for(j=J-2; j<=J+1; j++)
			{
				for(i=I-1; i<=I+1; i++)
				{
					col	= j*nx + i + numU;
					x	= 0.5*(D.x[i] + D.x[i+1]);
					y	= D.y[j+1];
					QT.row_indices[num_elements] = row;
					QT.column_indices[num_elements] = col;
					QT.values[num_elements] = alpha*delta(x-xB, y-yB, Dx)/Dx;
					num_elements++;
				}
			}
			row++;
		}
	}*/
	
	cusp::transpose(QT, Q);
}