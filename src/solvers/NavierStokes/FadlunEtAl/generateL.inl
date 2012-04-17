template <>
void FadlunEtAlSolver<host_memory>::generateL()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	real Cx0 = 0.0, Cx1 = 0.0, Cy0 = 0.0, Cy1 = 0.0,
	     scale = 0.0,
	     dx0 = 1.0, dx1 = 1.0, dy0 = 1.0, dy1 = 1.0;
	
	real *dx = thrust::raw_pointer_cast(&(domInfo->dx[0])),
	     *dy = thrust::raw_pointer_cast(&(domInfo->dy[0]));
	
	int  numU = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  I,
	     num_elements = 0;
			
	L.resize(numUV, numUV, 5*numUV-4*(nx+ny)+4);

	real nu = (*paramDB)["flow"]["nu"].get<real>();
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	
	/**
	* DO NOT FORGET
	* When you generate A as M - alpha*L
	* L has values that should not be multplied by alpha
	* These are the interpolation values of Fadlun et al
	*
	* The sign of the interpolation coeffs
	* is the opposite of that in the toy python code
	* since A is obtained by subtracting alpha*L from M
	*/
	
	///x-component
	for (int j=0; j < ny; j++)
	{
		if(j == 0)
		{
			dy0 = 0.5*dy[0];
			dy1 = 0.5*(dy[0]+dy[1]);
		}
		else if(j == ny-1)
		{
			dy0 = 0.5*(dy[ny-2]+dy[ny-1]);
			dy1 = 0.5*dy[ny-1];
		}
		else
		{
			dy0 = 0.5*(dy[j]+dy[j-1]);
			dy1 = 0.5*(dy[j]+dy[j+1]);
		}
		Cy0 = 2.0 * nu / ( dy1*(dy1+dy0) );
		Cy1 = 2.0 * nu / ( dy0*(dy1+dy0) );
		
		for (int i=0; i < nx-1; i++)
		{
			I = j*(nx-1) + i;										///< calculate the row of the matrix
			
			Cx0 = 2.0 * nu / ( dx[i+1]*(dx[i+1]+dx[i]) );
			Cx1 = 2.0 * nu / ( dx[i]*(dx[i+1]+dx[i]) );
			
			scale = 0.5*(dx[i+1]+dx[i]);					///< scaling factor (to obtain the normalised matrix)
		
			// south
			if(j>0)													///< no south coefficient for the bottom row
			{
				L.row_indices[num_elements] = I;
				L.column_indices[num_elements] = I - (nx-1);
				
				if(tags[I] == -1)
					L.values[num_elements] = Cy1 * scale / dy[j-1];
				else if(tags[I]==I-(nx-1))
					L.values[num_elements] = coeffs[I]/dt * scale / dy[j-1];
				else
					L.values[num_elements] = 0.0;
					
				num_elements++;
			}
			// west
			if(i>0)													///< no west coefficient for the leftmost column			
			{
				L.row_indices[num_elements] = I;
				L.column_indices[num_elements] = I - 1;
				
				if(tags[I] == -1)
					L.values[num_elements] = Cx1 * scale / dy[j];
				else if(tags[I]==I-1)
					L.values[num_elements] = coeffs[I]/dt * scale / dy[j];
				else
					L.values[num_elements] = 0.0;
					
				num_elements++;
			}
			// centre
			L.row_indices[num_elements] = I;
			L.column_indices[num_elements] = I;
			if(tags[I] == -1)
				L.values[num_elements] = (-Cy0-Cy1-Cx0-Cx1) * scale / dy[j];
			else
				L.values[num_elements] = 0.0;
			num_elements++;
			// east
			if(i<nx-2)											// no east coefficient for the rightmost column
			{
				L.row_indices[num_elements] = I;
				L.column_indices[num_elements] = I + 1;
				if(tags[I] == -1)
					L.values[num_elements] = Cx0 * scale / dy[j];
				else if(tags[I]=I+1)
					L.values[num_elements] = coeffs[I]/dt * scale / dy[j];
				else
					L.values[num_elements] = 0.0;
				num_elements++;
			}
			// north
			if(j<ny-1)
			{
				L.row_indices[num_elements] = I;
				L.column_indices[num_elements] = I + (nx-1);
				if(tags[I] == -1)	
					L.values[num_elements] = Cy0 * scale / dy[j+1];
				else if(tags[I]==I+(nx-1))
					L.values[num_elements] = coeffs[I]/dt * scale / dy[j+1];
				else
					L.values[num_elements] = 0.0;
				num_elements++;
			}
		}
	}
	
	///y-component
	for (int j=0; j < ny-1; j++)
	{
		Cy0 = 2.0 * nu / ( dy[j+1]*(dy[j+1]+dy[j]) );
		Cy1 = 2.0 * nu / ( dy[j]*(dy[j+1]+dy[j]) );
		
		for (int i=0; i < nx; i++)
		{
			I = j*nx + i + numU;		// calculate the row of the matrix
			
			if((i>0) && (i<nx-1))
			{
				dx0 = 0.5*(dx[i]+dx[i-1]);
				dx1 = 0.5*(dx[i]+dx[i+1]);
			}
			else if(i==0)
			{
				dx0 = 0.5*dx[i];
				dx1 = 0.5*(dx[i]+dx[i+1]);
			}
			else if(i==nx-1)
			{
				dx0 = 0.5*(dx[i]+dx[i-1]);
				dx1 = 0.5*dx[i];
			}
			Cx0 = 2.0 * nu /( dx1*(dx1+dx0) );
			Cx1 = 2.0 * nu /( dx0*(dx1+dx0) );
			
			scale = 0.5*(dy[j+1]+dy[j]);	// scaling factor (to obtain the normalised matrix)
			
			// south
			if(j>0)													// no south coefficient for the bottom row
			{
				L.row_indices[num_elements] = I;
				L.column_indices[num_elements] = I - nx;
				if(tags[I] == -1)
					L.values[num_elements] = Cy1 * scale / dx[i];
				else if(tags[I]==I-nx)
					L.values[num_elements] = coeffs[I]/dt * scale / dx[i];
				else
					L.values[num_elements] = 0.0;
				num_elements++;
			}
			// west
			if(i>0)													// no west coefficient for the leftmost column			
			{
				L.row_indices[num_elements] = I;
				L.column_indices[num_elements] = I - 1;
				if(tags[I] == -1)
					L.values[num_elements] = Cx1 * scale / dx[i-1];
				else if(tags[I]==I-1)
					L.values[num_elements] = coeffs[I]/dt * scale / dx[i-1];
				else
					L.values[num_elements] = 0.0;
				num_elements++;
			}
			// centre
			L.row_indices[num_elements] = I;
			L.column_indices[num_elements] = I;
			if(tags[I] == -1)
				L.values[num_elements] = (-Cy0-Cy1-Cx0-Cx1) * scale / dx[i];
			else
				L.values[num_elements] = 0.0;
			num_elements++;
			// east
			if(i<nx-1)												// no east coefficient for the rightmost column
			{
				L.row_indices[num_elements] = I;
				L.column_indices[num_elements] = I + 1;
				if(tags[I] == -1)
					L.values[num_elements] = Cx0 * scale / dx[i+1];
				else if(tags[I]==I+1)
					L.values[num_elements] = coeffs[I]/dt * scale / dx[i+1];
				else
					L.values[num_elements] = 0.0;
				num_elements++;
			}
			// north
			if(j<ny-2)
			{
				L.row_indices[num_elements] = I;
				L.column_indices[num_elements] = I + nx;
				if(tags[I] == -1)
					L.values[num_elements] = Cy0 * scale / dx[i];
				else if(tags[I] == I+nx)
					L.values[num_elements] = coeffs[I]/dt * scale / dx[i];
				else
					L.values[num_elements] = 0.0;
				num_elements++;
			}
		}
	}
}

template <>
void FadlunEtAlSolver<device_memory>::generateL()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	real Cx0 = 0.0, Cx1 = 0.0, Cy0 = 0.0, Cy1 = 0.0,
	     scale = 0.0,
	     dx0 = 1.0, dx1 = 1.0, dy0 = 1.0, dy1 = 1.0;
	
	real *dx = thrust::raw_pointer_cast(&(domInfo->dx[0])),
	     *dy = thrust::raw_pointer_cast(&(domInfo->dy[0]));
	
	int  numU = (nx-1)*ny;
	int  numUV = numU + nx*(ny-1);
	int  I,
	     num_elements = 0;
			
	cooH LHost(numUV, numUV, 5*numUV-4*(nx+ny)+4);

	real nu = (*paramDB)["flow"]["nu"].get<real>();
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	
	/**
	* DO NOT FORGET
	* When you generate A as M - alpha*L
	* L has values that should not be multplied by alpha
	* These are the interpolation values of Fadlun et al
	*
	* The sign of the interpolation coeffs
	* is the opposite of that in the toy python code
	* since A is obtained by subtracting alpha*L from M
	*/
	
	///x-component
	for (int j=0; j < ny; j++)
	{
		if(j == 0)
		{
			dy0 = 0.5*dy[0];
			dy1 = 0.5*(dy[0]+dy[1]);
		}
		else if(j == ny-1)
		{
			dy0 = 0.5*(dy[ny-2]+dy[ny-1]);
			dy1 = 0.5*dy[ny-1];
		}
		else
		{
			dy0 = 0.5*(dy[j]+dy[j-1]);
			dy1 = 0.5*(dy[j]+dy[j+1]);
		}
		Cy0 = 2.0 * nu / ( dy1*(dy1+dy0) );
		Cy1 = 2.0 * nu / ( dy0*(dy1+dy0) );
		
		for (int i=0; i < nx-1; i++)
		{
			I = j*(nx-1) + i;										///< calculate the row of the matrix
			
			Cx0 = 2.0 * nu / ( dx[i+1]*(dx[i+1]+dx[i]) );
			Cx1 = 2.0 * nu / ( dx[i]*(dx[i+1]+dx[i]) );
			
			scale = 0.5*(dx[i+1]+dx[i]);					///< scaling factor (to obtain the normalised matrix)
		
			// south
			if(j>0)													///< no south coefficient for the bottom row
			{
				LHost.row_indices[num_elements] = I;
				LHost.column_indices[num_elements] = I - (nx-1);
				
				if(tags[I] == -1)
					LHost.values[num_elements] = Cy1 * scale / dy[j-1];
				else if(tags[I]==I-(nx-1))
					LHost.values[num_elements] = coeffs[I]/dt * scale / dy[j-1];
				else
					LHost.values[num_elements] = 0.0;
					
				num_elements++;
			}
			// west
			if(i>0)													///< no west coefficient for the leftmost column			
			{
				LHost.row_indices[num_elements] = I;
				LHost.column_indices[num_elements] = I - 1;
				
				if(tags[I] == -1)
					LHost.values[num_elements] = Cx1 * scale / dy[j];
				else if(tags[I]==I-1)
					LHost.values[num_elements] = coeffs[I]/dt * scale / dy[j];
				else
					LHost.values[num_elements] = 0.0;
					
				num_elements++;
			}
			// centre
			LHost.row_indices[num_elements] = I;
			LHost.column_indices[num_elements] = I;
			if(tags[I] == -1)
				LHost.values[num_elements] = (-Cy0-Cy1-Cx0-Cx1) * scale / dy[j];
			else
				LHost.values[num_elements] = 0.0;
			num_elements++;
			// east
			if(i<nx-2)											// no east coefficient for the rightmost column
			{
				LHost.row_indices[num_elements] = I;
				LHost.column_indices[num_elements] = I + 1;
				if(tags[I] == -1)
					LHost.values[num_elements] = Cx0 * scale / dy[j];
				else if(tags[I]=I+1)
					LHost.values[num_elements] = coeffs[I]/dt * scale / dy[j];
				else
					LHost.values[num_elements] = 0.0;
				num_elements++;
			}
			// north
			if(j<ny-1)
			{
				LHost.row_indices[num_elements] = I;
				LHost.column_indices[num_elements] = I + (nx-1);
				if(tags[I] == -1)	
					LHost.values[num_elements] = Cy0 * scale / dy[j+1];
				else if(tags[I]==I+(nx-1))
					LHost.values[num_elements] = coeffs[I]/dt * scale / dy[j+1];
				else
					LHost.values[num_elements] = 0.0;
				num_elements++;
			}
		}
	}
	
	///y-component
	for (int j=0; j < ny-1; j++)
	{
		Cy0 = 2.0 * nu / ( dy[j+1]*(dy[j+1]+dy[j]) );
		Cy1 = 2.0 * nu / ( dy[j]*(dy[j+1]+dy[j]) );
		
		for (int i=0; i < nx; i++)
		{
			I = j*nx + i + numU;		// calculate the row of the matrix
			
			if((i>0) && (i<nx-1))
			{
				dx0 = 0.5*(dx[i]+dx[i-1]);
				dx1 = 0.5*(dx[i]+dx[i+1]);
			}
			else if(i==0)
			{
				dx0 = 0.5*dx[i];
				dx1 = 0.5*(dx[i]+dx[i+1]);
			}
			else if(i==nx-1)
			{
				dx0 = 0.5*(dx[i]+dx[i-1]);
				dx1 = 0.5*dx[i];
			}
			Cx0 = 2.0 * nu /( dx1*(dx1+dx0) );
			Cx1 = 2.0 * nu /( dx0*(dx1+dx0) );
			
			scale = 0.5*(dy[j+1]+dy[j]);	// scaling factor (to obtain the normalised matrix)
			
			// south
			if(j>0)													// no south coefficient for the bottom row
			{
				LHost.row_indices[num_elements] = I;
				LHost.column_indices[num_elements] = I - nx;
				if(tags[I] == -1)
					LHost.values[num_elements] = Cy1 * scale / dx[i];
				else if(tags[I]==I-nx)
					LHost.values[num_elements] = coeffs[I]/dt * scale / dx[i];
				else
					LHost.values[num_elements] = 0.0;
				num_elements++;
			}
			// west
			if(i>0)													// no west coefficient for the leftmost column			
			{
				LHost.row_indices[num_elements] = I;
				LHost.column_indices[num_elements] = I - 1;
				if(tags[I] == -1)
					LHost.values[num_elements] = Cx1 * scale / dx[i-1];
				else if(tags[I]==I-1)
					LHost.values[num_elements] = coeffs[I]/dt * scale / dx[i-1];
				else
					LHost.values[num_elements] = 0.0;
				num_elements++;
			}
			// centre
			LHost.row_indices[num_elements] = I;
			LHost.column_indices[num_elements] = I;
			if(tags[I] == -1)
				LHost.values[num_elements] = (-Cy0-Cy1-Cx0-Cx1) * scale / dx[i];
			else
				LHost.values[num_elements] = 0.0;
			num_elements++;
			// east
			if(i<nx-1)												// no east coefficient for the rightmost column
			{
				LHost.row_indices[num_elements] = I;
				LHost.column_indices[num_elements] = I + 1;
				if(tags[I] == -1)
					LHost.values[num_elements] = Cx0 * scale / dx[i+1];
				else if(tags[I]==I+1)
					LHost.values[num_elements] = coeffs[I]/dt * scale / dx[i+1];
				else
					LHost.values[num_elements] = 0.0;
				num_elements++;
			}
			// north
			if(j<ny-2)
			{
				LHost.row_indices[num_elements] = I;
				LHost.column_indices[num_elements] = I + nx;
				if(tags[I] == -1)
					LHost.values[num_elements] = Cy0 * scale / dx[i];
				else if(tags[I] == I+nx)
					LHost.values[num_elements] = coeffs[I]/dt * scale / dx[i];
				else
					LHost.values[num_elements] = 0.0;
				num_elements++;
			}
		}
	}
	L = LHost;
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
