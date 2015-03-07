/***************************************************************************//**
 * \file tagPoints.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c DirectForcingSolver to tag
 *        points near the immersed boundary using a ray-tracing algorithm.
 */


/**
 * \brief Tags the forcing nodes among the velocity nodes, i.e. the nodes at 
 *        which the velocity interpolation is performed.
 */
template <>
void DirectForcingSolver<host_memory>::tagPoints()
{
	logger.startTimer("tagPoints");

	real *bx = thrust::raw_pointer_cast(&(B.x[0])),
	     *by = thrust::raw_pointer_cast(&(B.y[0])),
	     *uB = thrust::raw_pointer_cast(&(B.uB[0])),
	     *vB = thrust::raw_pointer_cast(&(B.vB[0]));
	     
	tagPoints(bx, by, uB, vB);

	logger.stopTimer("tagPoints");
}

/**
 * \brief Tags the forcing nodes among the velocity nodes, i.e. the nodes at 
 *        which the velocity interpolation is performed.
 */
template <>
void DirectForcingSolver<device_memory>::tagPoints()
{
	logger.startTimer("tagPoints");
	
	// transferring boundary point coordinates to the host
	vecH bxH(B.totalPoints), byH(B.totalPoints), uBH(B.totalPoints), vBH(B.totalPoints);
	bxH = B.x;
	byH = B.y;
	uBH = B.uB;
	vBH = B.vB;

	// creating raw pointers
	real *bx = thrust::raw_pointer_cast(&(bxH[0])),
	     *by = thrust::raw_pointer_cast(&(byH[0])),
	     *uB = thrust::raw_pointer_cast(&(uBH[0])),
	     *vB = thrust::raw_pointer_cast(&(vBH[0]));
	
	//tagPoints(bx, by);
	tagPoints(bx, by, uB, vB);
	
	// transferring tag and coeffs data to the device
	tagsD    = tags;
	tags2D   = tags2;
	coeffsD  = coeffs;
	coeffs2D = coeffs2;
	uvD      = uv;
	
	logger.stopTimer("tagPoints");
}

// Bilinear Fadlun1c-type interpolation outside the body, for a moving body.
/**
 * \brief Tags all the forcing nodes required for the type of linear
 *        interpolation explained in the paper by Fadlun et al. (2000).
 *
 * It uses a raytracing algorithm to detect points that are near the boundary,
 * and just outside it. For more information about the algorithm, read the 
 * section on ray-crossings in the Search and Intersection chapter of the 
 * book Computational Geometry in C by Joseph O'Rourke.
 *
 * \param bx host array of the x-coordinates of the boundary points
 * \param by host array of the y-coordinates of the boundary points
 * \param uB host array of the x-components of the boundary velocities
 * \param vB host array of the y-components of the boundary velocities
 */
template <typename memoryType>
void DirectForcingSolver<memoryType>::tagPoints(real *bx, real *by, real *uB, real *vB)
{
	int  nx = NavierStokesSolver<memoryType>::domInfo->nx,
	     ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	real *xu = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->xu[0])),
	     *yu = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->yu[0]));

	parameterDB &db = *NavierStokesSolver<memoryType>::paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	interpolationType interpType = db["simulation"]["interpolationType"].get<interpolationType>();

	int  I;
	int  bottom, top, left, right;
	real eps = 1.e-10;
	bool outsideX, outsideY;
	int  bdryFlagX, bdryFlagY, bdryFlag2X, bdryFlag2Y;
	real a, b;
	real cfX, cfY, cf2X, cf2Y;
	real etaX, etaY;
	real uvX, uvY;
	bool flag;
	real x, y;
	int  k, l;
	int  totalPoints = NSWithBody<memoryType>::B.totalPoints;

	std::ofstream mask((folder+"/mask.txt").c_str());
	std::ofstream eta((folder+"/eta_u.txt").c_str());
	
	// tag points at which u is evaluated
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			// index of the current point on the u-grid
			I = j*(nx-1)+i;
			
			// tags and coefficients
			tags[I] = -1;
			tags2[I] = -1;
			coeffs[I] = 0.0;
			coeffs2[I] = 0.0;
			
			// initial indices of the points on the body that define the segment under consideration
			k = totalPoints-1;
			l = 0;
			
			outsideX = true;
			outsideY = true;
			bdryFlagX = -1;  // stores if a point is near the boundary
			bdryFlagY = -1;
			bdryFlag2X = -1;
			bdryFlag2Y = -1;
			uvX = 0.0;
			uvY = 0.0;
			cfX = 0.0;
			cfY = 0.0;
			cf2X = 0.0;
			cf2Y = 0.0;
			flag = false;
			
			// cycle through all the segments on the body surface
			while(l<totalPoints && !flag)
			{
				// figure out which of the two end points of the segment are at the bottom and the left
				if (by[k] > by[l])
				{
					bottom = l;
					top = k;
				}
				else
				{
					bottom = k;
					top = l;
				}
				if (bx[k] > bx[l])
				{
					left = l;
					right = k;
				}
				else
				{
					left = k;
					right = l;
				}
				
				// consider rays along the x-direction
				// if the ray intersects the boundary segment (top endpoint must be strictly above the ray; bottom can be on or below the ray)
				if (by[bottom]-eps < yu[j] && by[top]-eps > yu[j])
				{
					// if the segment is not parallel to the x-direction
					if (fabs(by[l]-by[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						x = bx[k] + (bx[l]-bx[k]) * (yu[j]-by[k])/(by[l]-by[k]);
						
						// calculate the body velocity at the point of intersection
						uvX = uB[k] + (uB[l]-uB[k]) * (yu[j]-by[k])/(by[l]-by[k]);
					
						// if the point of intersection coincides with the grid point
						if (fabs(x-xu[i]) < eps)
						{
							outsideX   = true;
							bdryFlagX = I;
							bdryFlag2X= I;
							etaX      = 0.0;
							cfX       = 0.0;
							cf2X      = 0.0;
							flag      = true; // flag is true when the point of intersection coincides with the grid point
						}
						// if the point of intersection lies to the right of the grid point (right-facing ray intersects the boundary)
				 		else if (x > xu[i]+eps)
							outsideX = !outsideX;
						
						// if the point of intersection is in the cell to the immediate left of the grid point
						if (x>xu[i-1]+eps && x<xu[i]-eps)
						{
							bdryFlagX  = I+1;
							bdryFlag2X = I+2;
							a = xu[i]-x;
							b = xu[i+1]-xu[i];
							etaX = a/b;
							switch(interpType)
							{
								case CONSTANT : cfX = 0.0; cf2X = 0.0; break;
								case LINEAR   : cfX = a/(a+b); cf2X = 0.0; break;
								case QUADRATIC: cfX = 2*a/(a+b); cf2X = -a/(a+2*b); break;
								default: cfX = a/(a+b); cf2X = 0.0; break;
							}
						}
						// if the point of intersection is in the cell to the immediate right of the grid point
						else if (x>xu[i]+eps && x<xu[i+1]-eps)
						{
							bdryFlagX  = I-1;
							bdryFlag2X = I-2;
							a = x-xu[i];
							b = xu[i]-xu[i-1];
							etaX = a/b;
							switch(interpType)
							{
								case CONSTANT : cfX = 0.0; cf2X = 0.0; break;
								case LINEAR   : cfX = a/(a+b); cf2X = 0.0; break;
								case QUADRATIC: cfX = 2*a/(a+b); cf2X = -a/(a+2*b); break;
								default: cfX = a/(a+b); cf2X = 0.0; break;
							}
						}
					}
				}
				// consider rays along the y-direction
				// if the ray intersects the boundary segment (right endpoint must be strictly to the right of ray; left can be on or to the left of the ray)
				if ( (bx[left]-eps < xu[i]) && (bx[right]-eps > xu[i]) && ( !flag ) ) // no need to do this part if flag is false
				{
					// if the segment is not parallel to the y-direction
					if (fabs(bx[l]-bx[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						y = by[k] + (by[l]-by[k]) * (xu[i]-bx[k]) / (bx[l]-bx[k]);
						
						// calculate the body velocity at the point of intersection
						uvY = uB[k] + (uB[l]-uB[k]) * (xu[i]-bx[k])/(bx[l]-bx[k]);
						
						// if the point of intersection coincides with the grid point
						if (fabs(y-yu[j]) < eps)
						{
							outsideY  = true; // then the point is considered to be outside the grid
							bdryFlagY = I;    // the point is considered to be a forcing point, with index I
							bdryFlag2Y= I;
							etaY      = 0.0;
							cfY       = 0.0;  // the coefficient for the linear interpolation during forcing
							cf2Y      = 0.0;
							flag      = true; // flag is true when the point of intersection coincides with the grid point
						}
						// if the point of intersection lies to the top of the grid point
				 		else if (y > yu[j]+eps)
							outsideY = !outsideY; // then flip if inside or outside (start with true, i.e. outside)
							
						// if point of intersection is just below the concerned grid point
						if (y>yu[j-1]+eps && y<yu[j]-eps)
						{
							bdryFlagY = I+(nx-1);
							bdryFlag2Y= I+2*(nx-1);
							a = yu[j]-y;
							b = yu[j+1]-yu[j];
							etaY = a/b;
							switch(interpType)
							{
								case CONSTANT : cfY = 0.0; cf2Y = 0.0; break;
								case LINEAR   : cfY = a/(a+b); cf2Y = 0.0; break;
								case QUADRATIC: cfY = 2*a/(a+b); cf2Y = -a/(a+2*b); break;
								default: cfY = a/(a+b); cf2Y = 0.0; break;
							}
						}
						// if point of intersection is just above the concerned grid point
						else if (y>yu[j]+eps && y<yu[j+1]-eps)
						{
							bdryFlagY = I-(nx-1);
							bdryFlag2Y= I-2*(nx-1);
							a = y-yu[j];
							b = yu[j]-yu[j-1];
							etaY = a/b;
							switch(interpType)
							{
								case CONSTANT : cfY = 0.0; cf2Y = 0.0; break;
								case LINEAR   : cfY = a/(a+b); cf2Y = 0.0; break;
								case QUADRATIC: cfY = 2*a/(a+b); cf2Y = -a/(a+2*b); break;
								default: cfY = a/(a+b); cf2Y = 0.0; break;
							}
						}
					}
				}
				k = l;
				l = l+1;
			}
			if (outsideX && bdryFlagX>=0)
			{
				tags[I]    = bdryFlagX;
				tags2[I]   = bdryFlag2X;
				coeffs[I]  = cfX;
				coeffs2[I] = cf2X;
				uv[I]      = uvX;
				eta << etaX << std::endl;
			}
			else if (outsideY && bdryFlagY>=0)
			{					
				tags[I]    = bdryFlagY;
				tags2[I]   = bdryFlag2Y;
				coeffs[I]  = cfY;
				coeffs2[I] = cf2Y;
				uv[I]      = uvY;
				eta << etaY << std::endl;
			}
			mask << ((outsideX || outsideY)? 1 : 0) << std::endl;
		}
	}
	eta.close();
	
	std::ofstream file((folder+"/tagx.txt").c_str());
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			I = j*(nx-1)+i;
			if (tags[I] >= 0)
			{
				file << xu[i] << '\t' << yu[j] << std::endl;
				if (tags[I] == I-1)
					file << xu[i-1] << '\t' << yu[j] << std::endl;
				else if (tags[I] == I+1)
					file << xu[i+1] << '\t' << yu[j] << std::endl;
				else if (tags[I] == I-(nx-1))
					file << xu[i] << '\t' << yu[j-1] << std::endl;
				else if (tags[I] == I+(nx-1))
					file << xu[i] << '\t' << yu[j+1] << std::endl;
				file << std::endl;
			}
		}
	}
	file.close();

	eta.open((folder+"/eta_v.txt").c_str());

	real *xv = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->xv[0])),
	     *yv = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->yv[0]));
	
	// tag points at which v is evaluated
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			// index of the current point on the u-grid
			I = j*nx+i + (nx-1)*ny;

			// tags and coefficients
			tags[I] = -1;
			coeffs[I] = 0.0;

			// initial indices of the points on the body that define the segment under consideration
			k = totalPoints-1;
			l = 0;
			
			outsideX = true;
			outsideY = true;
			bdryFlagX = -1;
			bdryFlagY = -1;
			uvX = 0.0;
			uvY = 0.0;
			cfX = 0.0;
			cfY = 0.0;
			flag = false;
			
			while(l<totalPoints)
			{
				if (by[k] > by[l])
				{
					bottom = l;
					top = k;
				}
				else
				{
					bottom = k;
					top = l;
				}
				if (bx[k] > bx[l])
				{
					left = l;
					right = k;
				}
				else
				{
					left = k;
					right = l;
				}
				// consider rays along the x-direction
				/**
				* if the ray intersects the boundary segment
				* top endpoint must be strictly above the ray
				* bottom can be on or below the ray
				*/
				if (by[bottom]-eps < yv[j] && by[top]-eps > yv[j] && !flag)
				{
					// if the segment is not parallel to the ray
					if (fabs(by[l]-by[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						x = bx[k] + (bx[l]-bx[k]) * (yv[j]-by[k])/(by[l]-by[k]);
						// calculate the body velocity at the point of intersection
						uvX = vB[k] + (vB[l]-vB[k]) * (yv[j]-by[k])/(by[l]-by[k]);

						// if the point of intersection coincides with the grid point
						if (fabs(x-xv[i]) < eps)
						{
							outsideX   = true;							
							bdryFlagX = I;
							bdryFlag2X= I;
							etaX      = 0.0;
							cfX       = 0.0;
							cf2X      = 0.0;
							flag      = true;
						}
						// if the point of intersection lies to the right of the grid point
				 		else if (x > xv[i]+eps)
							outsideX = !outsideX;

						// if the point of intersection is in the cell to the immediate right of the grid point
						if (x>xv[i-1]+eps && x<xv[i]-eps)
						{
							bdryFlagX  = I+1;
							bdryFlag2X = I+2;
							a = xv[i]-x;
							b = xv[i+1]-xv[i];
							etaX = a/b;
							switch(interpType)
							{
								case CONSTANT : cfX = 0.0; cf2X = 0.0; break;
								case LINEAR   : cfX = a/(a+b); cf2X = 0.0; break;
								case QUADRATIC: cfX = 2*a/(a+b); cf2X = -a/(a+2*b); break;
								default: cfX = a/(a+b); cf2X = 0.0; break;
							}
						}
						// if the point of intersection is in the cell to the immediate left of the grid point
						else if (x>xv[i]+eps && x<xv[i+1]-eps)
						{
							bdryFlagX  = I-1;
							bdryFlag2X = I-2;
							a = x-xv[i];
							b = xv[i]-xv[i-1];
							etaX = a/b;
							switch(interpType)
							{
								case CONSTANT : cfX = 0.0; cf2X = 0.0; break;
								case LINEAR   : cfX = a/(a+b); cf2X = 0.0; break;
								case QUADRATIC: cfX = 2*a/(a+b); cf2X = -a/(a+2*b); break;
								default: cfX = a/(a+b); cf2X = 0.0; break;
							}
						}
					}
				}
				// consider rays along the y-direction
				if (bx[left]-eps < xv[i] && bx[right]-eps > xv[i] && !flag)
				{
					// if the segment is not parallel to the ray
					if (fabs(bx[l]-bx[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						y = by[k] + (by[l]-by[k]) * (xv[i]-bx[k])/(bx[l]-bx[k]);
						// calculate the body velocity at the point of intersectioin
						uvY = vB[k] + (vB[l]-vB[k]) * (xv[i]-bx[k])/(bx[l]-bx[k]);

						// if the point of intersection coincides with the grid point
						if (fabs(y-yv[j]) < eps)
						{
							outsideY   = true;
							bdryFlagY = I;
							bdryFlag2Y= I;
							etaY      = 0.0;
							cfY       = 0.0;
							cf2Y      = 0.0;
							flag      = true;
						}
						// if the point of intersection lies to the top of the grid point
				 		else if (y > yv[j]+eps)
							outsideY = !outsideY;

						if (y>yv[j-1]+eps && y<yv[j]-eps)
						{
							bdryFlagY  = I+nx;
							bdryFlag2Y = I+2*nx;
							a = yv[j]-y;
							b = yv[j+1]-yv[j];
							etaY = a/b;
							switch(interpType)
							{
								case CONSTANT : cfY = 0.0; cf2Y = 0.0; break;
								case LINEAR   : cfY = a/(a+b); cf2Y = 0.0; break;
								case QUADRATIC: cfY = 2*a/(a+b); cf2Y = -a/(a+2*b); break;
								default: cfY = a/(a+b); cf2Y = 0.0; break;
							}
						}
						else if (y>yv[j]+eps && y<yv[j+1]-eps)
						{
							bdryFlagY  = I-nx;
							bdryFlag2Y = I-2*nx;
							a = y-yv[j];
							b = yv[j]-yv[j-1];
							etaY = a/b;
							switch(interpType)
							{
								case CONSTANT : cfY = 0.0; cf2Y = 0.0; break;
								case LINEAR   : cfY = a/(a+b); cf2Y = 0.0; break;
								case QUADRATIC: cfY = 2*a/(a+b); cf2Y = -a/(a+2*b); break;
								default: cfY = a/(a+b); cf2Y = 0.0; break;
							}
						}
					}
				}
				k = l;
				l = l+1;
			}
			if (outsideY && bdryFlagY>=0)
			{					
				tags[I]    = bdryFlagY;
				tags2[I]   = bdryFlag2Y;
				coeffs[I]  = cfY;
				coeffs2[I] = cf2Y;
				uv[I]      = uvY;
				eta << etaY << std::endl;
			}
			else if (outsideX && bdryFlagX>=0)
			{
				tags[I]    = bdryFlagX;
				tags2[I]   = bdryFlag2X;
				coeffs[I]  = cfX;
				coeffs2[I] = cf2X;
				uv[I]      = uvX;
				eta << etaX << std::endl;
			}
			mask << ((outsideX || outsideY)? 1 : 0) << std::endl;
		}
	}
	eta.close();
	mask.close();
	
	file.open((folder+"/tagy.txt").c_str());
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			I = j*nx+i+(nx-1)*ny;
			if (tags[I] >= 0 || tags[I] >= 0)
			{
				file << xv[i] << '\t' << yv[j] << std::endl;
				if (tags[I] == I-1)
					file << xv[i-1] << '\t' << yv[j] << std::endl;
				else if (tags[I] == I+1)
					file << xv[i+1] << '\t' << yv[j] << std::endl;
				else if (tags[I] == I-nx)
					file << xv[i] << '\t' << yv[j-1] << std::endl;
				else if (tags[I] == I+nx)
					file << xv[i] << '\t' << yv[j+1] << std::endl;
				file << std::endl;
			}
		}
	}
	file.close();
}
