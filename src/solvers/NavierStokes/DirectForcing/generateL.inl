/***************************************************************************//**
 * \file generateL.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c DirectForcingSolver to generate
 *        the Laplacian matrix.
 */


/*template <typename memoryType>
void DirectForcingSolver<memoryType>::generateL()
{
	
}*/

/**
 * \brief Sets up the Laplacian matrix for the implicit diffusion term. 
 *        Uses second order central differences on the non-uniform grid.
 */
template <>
void DirectForcingSolver<host_memory>::generateL() 
{
}

/**
 * \brief Sets up the Laplacian matrix for the implicit diffusion term. 
 *        Uses second order central differences on the non-uniform grid.
 */
template <>
void DirectForcingSolver<device_memory>::generateL() 
{
#if 1
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
	interpolationType interpType = (*paramDB)["simulation"]["interpolationType"].get<interpolationType>();
	
	bool notBdryNode = true;
	
	/**
	* DO NOT FORGET (IMPLEMENTED)
	* A is generated as M-alpha*LHost
	* But some values of L are the interpolation values of the direct forcing method
	* These should not be multplied by alpha
	*
	* The sign of the interpolation coeffs
	* is the opposite of that in the toy python code
	* since A is obtained by subtracting alpha*L from M
	* 
	* (NOT IMPLEMENTED YET)
	* Also remember to take alpha into account in BC1 for moving bodies
	*/

	// 1-D interpolation is performed
	// along the x-direction for the u-velocities
	// along the y-direction for the v-velocities
	
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
			// calculate the row of the matrix
			I = j*(nx-1) + i;
	
			// check if the node is on the boundary or not
			notBdryNode = (tags[I]==-1);
			
			Cx0 = 2.0 * nu / ( dx[i+1]*(dx[i+1]+dx[i]) );
			Cx1 = 2.0 * nu / ( dx[i]*(dx[i+1]+dx[i]) );
			
			// scaling factor (to obtain the normalised matrix)
			scale = 0.5*(dx[i+1]+dx[i]);
		
			// south
			if(j>0)  // no south coefficient for the bottom row
			{
				LHost.row_indices[num_elements] = I;
				LHost.column_indices[num_elements] = I-(nx-1);
				LHost.values[num_elements] = 0.0;
				
				if(notBdryNode)	// if the interpolation is definitely not taking place
					LHost.values[num_elements] = Cy1 * scale / dy[j-1];
				else if(tags[I]==I-(nx-1) && interpType==LINEAR)
					LHost.values[num_elements] = coeffs[I] /dt * scale / dy[j-1];
				else if((tags2[I]==I-2 || tags2[I]==I-2*(nx-1)) && interpType==QUADRATIC)
				{
					LHost.column_indices[num_elements] = tags2[I];
					LHost.values[num_elements] = coeffs2[I] /dt * scale / dy[j-1];
				}
					
				num_elements++;
			}
			// west
			if(i>0)  // no west coefficient for the leftmost column			
			{
				LHost.row_indices[num_elements] = I;
				LHost.column_indices[num_elements] = I-1;
				LHost.values[num_elements] = 0.0;
				
				if(notBdryNode)
					LHost.values[num_elements] = Cx1 * scale / dy[j];
				else if(tags[I]==I-1)
					LHost.values[num_elements] = coeffs[I] /dt * scale / dy[j];
				else if(tags[I]==I-(nx-1) && interpType==QUADRATIC)
				{
					LHost.column_indices[num_elements] = I-(nx-1);
					LHost.values[num_elements] = coeffs[I] /dt * scale / dy[j];
				}
					
				num_elements++;
			}
			// centre
			LHost.row_indices[num_elements] = I;
			LHost.column_indices[num_elements] = I;
			if(notBdryNode)
				LHost.values[num_elements] = (-Cy0-Cy1-Cx0-Cx1) * scale / dy[j];
			else
				LHost.values[num_elements] = 0.0;
			num_elements++;
			// east
			if(i<nx-2)  // no east coefficient for the rightmost column
			{
				LHost.row_indices[num_elements] = I;
				LHost.column_indices[num_elements] = I+1;
				LHost.values[num_elements] = 0.0;

				if(notBdryNode)
					LHost.values[num_elements] = Cx0 * scale / dy[j];
				else if(tags[I]==I+1)
					LHost.values[num_elements] = coeffs[I] /dt * scale / dy[j];
				else if(tags[I]==I+(nx-1) && interpType==QUADRATIC)
				{
					LHost.column_indices[num_elements] = I+(nx-1);
					LHost.values[num_elements] = coeffs[I] /dt * scale / dy[j];
				}

				num_elements++;
			}
			// north
			if(j<ny-1)
			{
				LHost.row_indices[num_elements] = I;
				LHost.column_indices[num_elements] = I+(nx-1);
				LHost.values[num_elements] = 0.0;

				if(notBdryNode)	
					LHost.values[num_elements] = Cy0 * scale / dy[j+1];
				else if(tags[I]==I+(nx-1) && interpType==LINEAR)
					LHost.values[num_elements] = coeffs[I] /dt * scale / dy[j+1];
				else if((tags2[I]==I+2 || tags2[I]==I+2*(nx-1)) && interpType==QUADRATIC)
				{
					LHost.column_indices[num_elements] = tags2[I];
					LHost.values[num_elements] = coeffs2[I] /dt * scale / dy[j+1];
				}
				
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
			I = j*nx + i + numU;  // calculate the row of the matrix
			notBdryNode = (tags[I]==-1);
			
			if((i>0) && (i<nx-1))
			{
				dx0 = 0.5*(dx[i]+dx[i-1]);
				dx1 = 0.5*(dx[i]+dx[i+1]);
			}
			else if(i == 0)
			{
				dx0 = 0.5*dx[i];
				dx1 = 0.5*(dx[i]+dx[i+1]);
			}
			else if(i == nx-1)
			{
				dx0 = 0.5*(dx[i]+dx[i-1]);
				dx1 = 0.5*dx[i];
			}
			Cx0 = 2.0 * nu /( dx1*(dx1+dx0) );
			Cx1 = 2.0 * nu /( dx0*(dx1+dx0) );
			
			scale = 0.5*(dy[j+1]+dy[j]);  // scaling factor (to obtain the normalised matrix)
			
			// south
			if(j>0)  // no south coefficient for the bottom row
			{
				LHost.row_indices[num_elements] = I;
				LHost.column_indices[num_elements] = I-nx;
				LHost.values[num_elements] = 0.0;

				if(notBdryNode)
					LHost.values[num_elements] = Cy1 * scale / dx[i];
				else if(tags[I]==I-nx && interpType==LINEAR)
					LHost.values[num_elements] = coeffs[I] /dt * scale / dx[i];
				else if((tags2[I]==I-2 || tags2[I]==I-2*nx) && interpType==QUADRATIC)
				{
					LHost.column_indices[num_elements] = tags2[I];
					LHost.values[num_elements] = coeffs2[I] /dt * scale / dx[i];
				}

				num_elements++;
			}
			// west
			if(i>0)  // no west coefficient for the leftmost column			
			{
				LHost.row_indices[num_elements] = I;
				LHost.column_indices[num_elements] = I-1;
				LHost.values[num_elements] = 0.0;

				if(notBdryNode)
					LHost.values[num_elements] = Cx1 * scale / dx[i-1];
				else if(tags[I]==I-1)
					LHost.values[num_elements] = coeffs[I] /dt * scale / dx[i-1];
				else if(tags[I]==I-nx && interpType==QUADRATIC)
				{
					LHost.column_indices[num_elements] = I-nx;
					LHost.values[num_elements] = coeffs[I] /dt * scale / dx[i-1];
				}

				num_elements++;
			}
			// centre
			LHost.row_indices[num_elements] = I;
			LHost.column_indices[num_elements] = I;
			if(notBdryNode)
				LHost.values[num_elements] = (-Cy0-Cy1-Cx0-Cx1) * scale / dx[i];
			else
				LHost.values[num_elements] = 0.0;
			num_elements++;
			// east
			if(i<nx-1)  // no east coefficient for the rightmost column
			{
				LHost.row_indices[num_elements] = I;
				LHost.column_indices[num_elements] = I+1;
				LHost.values[num_elements] = 0.0;

				if(notBdryNode)
					LHost.values[num_elements] = Cx0 * scale / dx[i+1];
				else if(tags[I]==I+1)
					LHost.values[num_elements] = coeffs[I] /dt * scale / dx[i+1];
				else if(tags[I]==I+nx && interpType==QUADRATIC)
				{
					LHost.column_indices[num_elements] = I+nx;
					LHost.values[num_elements] = coeffs[I] /dt * scale / dx[i+1];
				}

				num_elements++;
			}
			// north
			if(j<ny-2)
			{
				LHost.row_indices[num_elements] = I;
				LHost.column_indices[num_elements] = I+nx;
				LHost.values[num_elements] = 0.0;

				if(notBdryNode)
					LHost.values[num_elements] = Cy0 * scale / dx[i];
				else if(tags[I]==I+nx && interpType==LINEAR)
					LHost.values[num_elements] = coeffs[I] /dt * scale / dx[i];
				else if((tags2[I]==I+2 || tags2[I]==I+2*nx) && interpType==QUADRATIC)
				{
					LHost.column_indices[num_elements] = tags2[I];
					LHost.values[num_elements] = coeffs2[I] /dt * scale / dx[i];
				}
					
				num_elements++;
			}
		}
	}
	L = LHost;
#endif
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
