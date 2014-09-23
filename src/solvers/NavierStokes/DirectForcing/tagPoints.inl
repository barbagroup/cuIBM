/**
* @file tagPoints.inl
* @brief Tag points near the immsersed boundary using a ray-tracing algorithm
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
	tagsXD   = tagsX;
	tagsYD   = tagsY;
	coeffsXD = coeffsX;
	coeffsYD = coeffsY;
	uvXD     = uvX;
	uvYD     = uvY;
	
	logger.stopTimer("tagPoints");
}
#if 0
// Bilinear Fadlun1c-type interpolation inside the body.
template <typename memoryType>
void DirectForcingSolver<memoryType>::tagPoints(real *bx, real *by)
{
	int  nx = NavierStokesSolver<memoryType>::domInfo->nx,
	     ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	real *xu = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->xu[0])),
	     *yu = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->yu[0]));

	bool insideX, insideY, flag;
	int  bdryFlagX, bdryFlagY, k, l, I;
	int  bottom, top, left, right;
	real eps = 1.e-8;
	real cfX = 0.0, cfY = 0.0, x, y;
	
	// tag points at which u is evaluated
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			I = j*(nx-1)+i;
			tagsX[I] = -1;
			tagsY[I] = -1;
			coeffsX[I] = 0.0;
			coeffsY[I] = 1.0;
			insideX = false;
			insideY = false;
			k = NSWithBody<memoryType>::B.totalPoints - 1;
			l = 0;
			flag = false;
			bdryFlagX = -1;
			bdryFlagY = -1;
				
			while( l < NSWithBody<memoryType>::B.totalPoints)
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
				if (by[bottom]-eps < yu[j] && by[top]-eps > yu[j])
				{
					if (fabs(by[l]-by[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						x = bx[k] + (bx[l]-bx[k]) * (yu[j]-by[k]) / (by[l]-by[k]);
						// if the point of intersection coincides with the grid point
						if (fabs(x-xu[i]) < eps)
						{
							insideX   = true;
							bdryFlagX = I;
							cfX       = 0.0;
							flag      = true;
						}
						// if the point of intersection lies to the right of the grid point
				 		else if (x > xu[i]+eps)
							insideX = !insideX;
						// if the point of intersection is in the cell to the immediate left of the grid point
						if (x>xu[i-1]+eps && x<xu[i]-eps)
						{
							bdryFlagX = I+1;
							cfX = (xu[i]-x)/(xu[i+1]-x);
						}
						// if the point of intersection is in the cell to the immediate right of the grid point
						else if (x>xu[i]+eps && x<xu[i+1]-eps)
						{
							bdryFlagX = I-1;
							cfX = (xu[i]-x)/(xu[i-1]-x);
						}
					}
				}
				// consider rays along the y-direction
				if (bx[left]-eps < xu[i] && bx[right]-eps > xu[i] && !flag)
				{
					if (fabs(bx[l]-bx[k]) > eps)
					{
						y = by[k] + (by[l]-by[k]) * (xu[i]-bx[k]) / (bx[l]-bx[k]);
						if (fabs(y-yu[j]) < eps)
						{
							insideY   = true;
							bdryFlagY = I;
							cfY       = 0.0;
						}
				 		else if (y > yu[j]+eps)
							insideY = !insideY;
						// if point of intersection is below the concerned grid point
						if (y>yu[j-1]+eps && y<yu[j]-eps)
						{
							bdryFlagY = I+(nx-1);
							cfY = (yu[j]-y)/(yu[j+1]-y);
						}
						// if point of intersection is above the concerned grid point
						else if (y>yu[j]+eps && y<yu[j+1]-eps)
						{
							bdryFlagY = I-(nx-1);
							cfY = (yu[j]-y)/(yu[j-1]-y);
						}
					}
				}
				k = l;
				l = l+1;
			}
			if (insideX && bdryFlagX>=0)
			{
				tagsX[I]   = bdryFlagX;
				coeffsX[I] = cfX;
			}
			if (insideY && bdryFlagY>=0)
			{					
				tagsY[I]   = bdryFlagY;
				coeffsY[I] = cfY;
			}
		}
	}
	
	std::ofstream file("tagx.txt");
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			I = j*(nx-1)+i;
			if (tagsX[I] >= 0 || tagsY[I] >= 0 )
			{
				file << xu[i] << '\t' << yu[j] << std::endl;
				if (tagsX[I] == I-1)
					file << xu[i-1] << '\t' << yu[j] << std::endl;
				else if (tagsX[I] == I+1)
					file << xu[i+1] << '\t' << yu[j] << std::endl;
				if (tagsY[I] == I-(nx-1))
					file << xu[i] << '\t' << yu[j-1] << std::endl;
				else if (tagsY[I] == I+(nx-1))
					file << xu[i] << '\t' << yu[j+1] << std::endl;
				file << xu[i] << '\t' << yu[j] << std::endl;
			}
			file << std::endl;
		}
	}
	file.close();

	real *xv = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->xv[0])),
	     *yv = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->yv[0]));
	
	// tag points at which v is evaluated
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			I = j*nx+i + (nx-1)*ny;
			tagsX[I] = -1;
			tagsY[I] = -1;
			coeffsX[I] = 0.0;
			coeffsY[I] = 1.0;
			insideX = false;
			insideY = false;
			k = NSWithBody<memoryType>::B.totalPoints - 1;
			l = 0;
			flag = false;
			bdryFlagX = -1;
			bdryFlagY = -1;
			
			while( l < NSWithBody<memoryType>::B.totalPoints)
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
				if (by[bottom]-eps < yv[j] && by[top]-eps > yv[j])
				{
					if (fabs(by[l]-by[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						x = bx[k] + (bx[l]-bx[k]) * (yv[j]-by[k]) / (by[l]-by[k]);
						// if the point of intersection coincides with the grid point
						if (fabs(x-xv[i]) < eps)
						{
							insideX   = true;
							flag      = true;
							bdryFlagX = I;
							cfX       = 0.0;
						}
						// if the point of intersection lies to the right of the grid point
				 		else if (x > xv[i]+eps)
							insideX = !insideX;
						// if the point of intersection is in the cell to the immediate right of the grid point
						if (x>xv[i-1]+eps && x<xv[i]-eps)
						{
							bdryFlagX = I+1;
							cfX = (xv[i]-x)/(xv[i+1]-x);
						}
						// if the point of intersection is in the cell to the immediate left of the grid point
						else if (x>xv[i]+eps && x<xv[i+1]-eps)
						{
							bdryFlagX = I-1;
							cfX = (xv[i]-x)/(xv[i-1]-x);
						}
					}
				}
				// consider rays along the y-direction
				if (bx[left]-eps < xv[i] && bx[right]-eps > xv[i] && !flag)
				{
					if (fabs(bx[l]-bx[k]) > eps)
					{
						y = by[k] + (by[l]-by[k]) * (xv[i]-bx[k]) / (bx[l]-bx[k]);
						if (fabs(y-yv[j]) < eps)
						{
							insideY   = true;
							bdryFlagY = I;
							cfY       = 0.0;
							flag      = true;
						}
				 		else if (y > yv[j]+eps)
							insideY = !insideY;

						if (y>yv[j-1]+eps && y<yv[j]-eps)
						{
							bdryFlagY = I+nx;
							cfY = (yv[j]-y)/(yv[j+1]-y);
						}
						else if (y>yv[j]+eps && y<yv[j+1]-eps)
						{
							bdryFlagY = I-nx;
							cfY = (yv[j]-y)/(yv[j-1]-y);
						}
					}
				}
				k = l;
				l = l+1;
			}
			if (insideY && bdryFlagY>=0)
			{					
				tagsY[I]   = bdryFlagY;
				coeffsY[I] = cfY;
			}
			if (insideX && bdryFlagX>=0)
			{
				tagsX[I]   = bdryFlagX;
				coeffsX[I] = cfX;
			} 
		}
	}
	
	file.open("tagy.txt");
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			I = j*nx+i+(nx-1)*ny;
			if (tagsY[I] >= 0 || tagsX[I] >= 0)
			{
				file << xv[i] << '\t' << yv[j] << std::endl;
				if (tagsX[I] == I-1)
					file << xv[i-1] << '\t' << yv[j] << std::endl;
				else if (tagsX[I] == I+1)
					file << xv[i+1] << '\t' << yv[j] << std::endl;
				if (tagsY[I] == I-nx)
					file << xv[i] << '\t' << yv[j-1] << std::endl;
				else if (tagsY[I] == I+nx)
					file << xv[i] << '\t' << yv[j+1] << std::endl;
				file << xv[i] << '\t' << yv[j] << std::endl;
			}
			file << std::endl;
		}
	}
	file.close();
	
	file.open("body.txt");
	for(int k=0; k<NSWithBody<memoryType>::B.totalPoints; k++)
	{
		file << bx[k] << '\t' << by[k] << std::endl;
	}
	file.close();
}
#endif

/********************* ABOVE IS FOR INSIDE THE BODY ***************************/

/********************* BELOW IS FOR OUTSIDE THE BODY **************************/
#if 0
// Bilinear Fadlun1c-type interpolation outside the body.
template <typename memoryType>
void DirectForcingSolver<memoryType>::tagPoints(real *bx, real *by)
{
	int  nx = NavierStokesSolver<memoryType>::domInfo->nx,
	     ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	real *xu = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->xu[0])),
	     *yu = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->yu[0]));

	bool outsideX, outsideY, flag;
	int  bdryFlagX, bdryFlagY, k, l, I;
	int  bottom, top, left, right;
	real eps = 1.e-8;
	real cfX = 0.0, cfY = 0.0, x, y;
	
	// tag points at which u is evaluated
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			I = j*(nx-1)+i;
			tagsX[I] = -1;
			tagsY[I] = -1;
			coeffsX[I] = 0.0;
			coeffsY[I] = 1.0;
			outsideX = true;
			outsideY = true;
			k = NSWithBody<memoryType>::B.totalPoints - 1;
			l = 0;
			flag = false;
			bdryFlagX = -1;
			bdryFlagY = -1;
				
			while( l < NSWithBody<memoryType>::B.totalPoints)
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
				if (by[bottom]-eps < yu[j] && by[top]-eps > yu[j])
				{
					if (fabs(by[l]-by[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						x = bx[k] + (bx[l]-bx[k]) * (yu[j]-by[k]) / (by[l]-by[k]);
						// if the point of intersection coincides with the grid point
						if (fabs(x-xu[i]) < eps)
						{
							outsideX   = true;
							bdryFlagX = I;
							cfX       = 0.0;
							flag      = true;
						}
						// if the point of intersection lies to the right of the grid point
				 		else if (x > xu[i]+eps)
							outsideX = !outsideX;
						// if the point of intersection is in the cell to the immediate left of the grid point
						if (x>xu[i-1]+eps && x<xu[i]-eps)
						{
							bdryFlagX = I+1;
							cfX = (xu[i]-x)/(xu[i+1]-x);
						}
						// if the point of intersection is in the cell to the immediate right of the grid point
						else if (x>xu[i]+eps && x<xu[i+1]-eps)
						{
							bdryFlagX = I-1;
							cfX = (xu[i]-x)/(xu[i-1]-x);
						}
					}
				}
				// consider rays along the y-direction
				if ( ( bx[left]-eps < xu[i] ) && ( bx[right]-eps > xu[i] ) && ( !flag ) )
				{
					if (fabs(bx[l]-bx[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						y = by[k] + (by[l]-by[k]) * (xu[i]-bx[k]) / (bx[l]-bx[k]);
						// if the point of intersection coincides with the grid point
						if (fabs(y-yu[j]) < eps)
						{
							outsideY   = true;
							bdryFlagY = I;
							cfY       = 0.0;
						}
						// if the point of intersection lies to the top of the grid point
				 		else if (y > yu[j]+eps)
							outsideY = !outsideY;
						// if point of intersection is just below the concerned grid point
						if (y>yu[j-1]+eps && y<yu[j]-eps)
						{
							bdryFlagY = I+(nx-1);
							cfY = (yu[j]-y)/(yu[j+1]-y);
						}
						// if point of intersection is just above the concerned grid point
						else if (y>yu[j]+eps && y<yu[j+1]-eps)
						{
							bdryFlagY = I-(nx-1);
							cfY = (yu[j]-y)/(yu[j-1]-y);
						}
					}
				}
				k = l;
				l = l+1;
			}
			if (outsideX && bdryFlagX>=0)
			{
				tagsX[I]   = bdryFlagX;
				coeffsX[I] = cfX;
			}
			if (outsideY && bdryFlagY>=0)
			{					
				tagsY[I]   = bdryFlagY;
				coeffsY[I] = cfY;
			}
		}
	}
	
	std::ofstream file("tagx.txt");
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			I = j*(nx-1)+i;
			if (tagsX[I] >= 0 || tagsY[I] >= 0 )
			{
				file << xu[i] << '\t' << yu[j] << std::endl;
				if (tagsX[I] == I-1)
					file << xu[i-1] << '\t' << yu[j] << std::endl;
				else if (tagsX[I] == I+1)
					file << xu[i+1] << '\t' << yu[j] << std::endl;
				if (tagsY[I] == I-(nx-1))
					file << xu[i] << '\t' << yu[j-1] << std::endl;
				else if (tagsY[I] == I+(nx-1))
					file << xu[i] << '\t' << yu[j+1] << std::endl;
				file << xu[i] << '\t' << yu[j] << std::endl;
			}
			file << std::endl;
		}
	}
	file.close();

	real *xv = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->xv[0])),
	     *yv = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->yv[0]));
	
	// tag points at which v is evaluated
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			I = j*nx+i + (nx-1)*ny;
			tagsX[I] = -1;
			tagsY[I] = -1;
			coeffsX[I] = 0.0;
			coeffsY[I] = 1.0;
			outsideX = true;
			outsideY = true;
			k = NSWithBody<memoryType>::B.totalPoints - 1;
			l = 0;
			flag = false;
			bdryFlagX = -1;
			bdryFlagY = -1;
			
			while( l < NSWithBody<memoryType>::B.totalPoints)
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
				if (by[bottom]-eps < yv[j] && by[top]-eps > yv[j])
				{
					// if the segment is not parallel to the ray
					if (fabs(by[l]-by[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						x = bx[k] + (bx[l]-bx[k]) * (yv[j]-by[k]) / (by[l]-by[k]);
						// if the point of intersection coincides with the grid point
						if (fabs(x-xv[i]) < eps)
						{
							outsideX   = true;
							flag      = true;
							bdryFlagX = I;
							cfX       = 0.0;
						}
						// if the point of intersection lies to the right of the grid point
				 		else if (x > xv[i]+eps)
							outsideX = !outsideX;
						// if the point of intersection is in the cell to the immediate right of the grid point
						if (x>xv[i-1]+eps && x<xv[i]-eps)
						{
							bdryFlagX = I+1;
							cfX = (xv[i]-x)/(xv[i+1]-x);
						}
						// if the point of intersection is in the cell to the immediate left of the grid point
						else if (x>xv[i]+eps && x<xv[i+1]-eps)
						{
							bdryFlagX = I-1;
							cfX = (xv[i]-x)/(xv[i-1]-x);
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
						y = by[k] + (by[l]-by[k]) * (xv[i]-bx[k]) / (bx[l]-bx[k]);
						// if the point of intersection coincides with the grid point
						if (fabs(y-yv[j]) < eps)
						{
							outsideY   = true;
							bdryFlagY = I;
							cfY       = 0.0;
							flag      = true;
						}
						// if the point of intersection lies to the top of the grid point
				 		else if (y > yv[j]+eps)
							outsideY = !outsideY;

						if (y>yv[j-1]+eps && y<yv[j]-eps)
						{
							bdryFlagY = I+nx;
							cfY = (yv[j]-y)/(yv[j+1]-y);
						}
						else if (y>yv[j]+eps && y<yv[j+1]-eps)
						{
							bdryFlagY = I-nx;
							cfY = (yv[j]-y)/(yv[j-1]-y);
						}
					}
				}
				k = l;
				l = l+1;
			}
			if (outsideY && bdryFlagY>=0)
			{					
				tagsY[I]   = bdryFlagY;
				coeffsY[I] = cfY;
			}
			if (outsideX && bdryFlagX>=0)
			{
				tagsX[I]   = bdryFlagX;
				coeffsX[I] = cfX;
			} 
		}
	}
	
	file.open("tagy.txt");
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			I = j*nx+i+(nx-1)*ny;
			if (tagsY[I] >= 0 || tagsX[I] >= 0)
			{
				file << xv[i] << '\t' << yv[j] << std::endl;
				if (tagsX[I] == I-1)
					file << xv[i-1] << '\t' << yv[j] << std::endl;
				else if (tagsX[I] == I+1)
					file << xv[i+1] << '\t' << yv[j] << std::endl;
				if (tagsY[I] == I-nx)
					file << xv[i] << '\t' << yv[j-1] << std::endl;
				else if (tagsY[I] == I+nx)
					file << xv[i] << '\t' << yv[j+1] << std::endl;
				file << xv[i] << '\t' << yv[j] << std::endl;
			}
			file << std::endl;
		}
	}
	file.close();
	
	file.open("body.txt");
	for(int k=0; k<NSWithBody<memoryType>::B.totalPoints; k++)
	{
		file << bx[k] << '\t' << by[k] << std::endl;
	}
	file.close();
}
#endif

/*********************** Fadlun-1c for MOVING BODY ***************************/
#if 1
// Bilinear Fadlun1c-type interpolation outside the body, for a moving body.
template <typename memoryType>
//void DirectForcingSolver<memoryType>::tagPoints(real *bx, real *by)
void DirectForcingSolver<memoryType>::tagPoints(real *bx, real *by, real *uB, real *vB)
{
	int  nx = NavierStokesSolver<memoryType>::domInfo->nx,
	     ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	real *xu = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->xu[0])),
	     *yu = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->yu[0]));

	parameterDB &db = *NavierStokesSolver<memoryType>::paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();

	bool outsideX, outsideY, flag;
	int  bdryFlagX, bdryFlagY, k, l, I;
	int  bottom, top, left, right;
	real eps = 1.e-8;
	real cfX = 0.0, cfY = 0.0, x, y;
	
	// tag points at which u is evaluated
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			// index of the current point on the u-grid
			I = j*(nx-1)+i;
			
			// tags and coefficients
			tagsX[I] = -1;
			tagsY[I] = -1;
			coeffsX[I] = 0.0;
			coeffsY[I] = 1.0;
			
			// initial indices of the points on the body that define the segment under consideration
			k = NSWithBody<memoryType>::B.totalPoints - 1;
			l = 0;
			
			outsideX = true;
			outsideY = true;
			bdryFlagX = -1;  // stores if a point is near the boundary
			bdryFlagY = -1;
			flag = false;
			
			// cycle through all the segments on the body surface
			while( l < NSWithBody<memoryType>::B.totalPoints)
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
					// if the segments is not parallel to the x-direction
					if (fabs(by[l]-by[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						x = bx[k] + (bx[l]-bx[k]) * (yu[j]-by[k])/(by[l]-by[k]);
						
						// calculate the body velocity at the point of intersection
						uvX[I] = uB[k] + (uB[l]-uB[k]) * (yu[j]-by[k])/(by[l]-by[k]);
						
						// if the point of intersection coincides with the grid point
						if (fabs(x-xu[i]) < eps)
						{
							outsideX   = true;
							bdryFlagX = I;
							cfX       = 0.0;
							flag      = true; // flag is true when the point of intersection coincides with the grid point
						}
						// if the point of intersection lies to the right of the grid point (right-facing ray intersects the boundary)
				 		else if (x > xu[i]+eps)
							outsideX = !outsideX;
						
						// if the point of intersection is in the cell to the immediate left of the grid point
						if (x>xu[i-1]+eps && x<xu[i]-eps)
						{
							bdryFlagX = I+1;
							cfX = (xu[i]-x)/(xu[i+1]-x);
						}
						// if the point of intersection is in the cell to the immediate right of the grid point
						else if (x>xu[i]+eps && x<xu[i+1]-eps)
						{
							bdryFlagX = I-1;
							cfX = (xu[i]-x)/(xu[i-1]-x);
						}
					}
				}
				// consider rays along the y-direction
				// if the ray intersects the boundary segment (right endpoint must be strictly to the right of ray; left can be on or to the left of the ray)
				if ( ( bx[left]-eps < xu[i] ) && ( bx[right]-eps > xu[i] ) && ( !flag ) ) // no need to do this part if flag is false
				{
					if (fabs(bx[l]-bx[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						y = by[k] + (by[l]-by[k]) * (xu[i]-bx[k]) / (bx[l]-bx[k]);
						
						// calculate the body velocity at the point of intersection
						uvY[I] = uB[k] + (uB[l]-uB[k]) * (xu[i]-bx[k])/(bx[l]-bx[k]);
						
						// if the point of intersection coincides with the grid point
						if (fabs(y-yu[j]) < eps)
						{
							outsideY   = true;
							bdryFlagY = I;
							cfY       = 0.0;
						}
						// if the point of intersection lies to the top of the grid point
				 		else if (y > yu[j]+eps)
							outsideY = !outsideY;
							
						// if point of intersection is just below the concerned grid point
						if (y>yu[j-1]+eps && y<yu[j]-eps)
						{
							bdryFlagY = I+(nx-1);
							cfY = (yu[j]-y)/(yu[j+1]-y);
						}
						// if point of intersection is just above the concerned grid point
						else if (y>yu[j]+eps && y<yu[j+1]-eps)
						{
							bdryFlagY = I-(nx-1);
							cfY = (yu[j]-y)/(yu[j-1]-y);
						}
					}
				}
				k = l;
				l = l+1;
			}
			if (outsideX && bdryFlagX>=0)
			{
				tagsX[I]   = bdryFlagX;
				coeffsX[I] = cfX;
			}
			if (outsideY && bdryFlagY>=0)
			{					
				tagsY[I]   = bdryFlagY;
				coeffsY[I] = cfY;
			}
		}
	}
	
	std::ofstream file((folder+"/tagx.txt").c_str());
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			I = j*(nx-1)+i;
			if (tagsX[I] >= 0 || tagsY[I] >= 0 )
			{
				file << xu[i] << '\t' << yu[j] << std::endl;
				if (tagsX[I] == I-1)
					file << xu[i-1] << '\t' << yu[j] << std::endl;
				else if (tagsX[I] == I+1)
					file << xu[i+1] << '\t' << yu[j] << std::endl;
				if (tagsY[I] == I-(nx-1))
					file << xu[i] << '\t' << yu[j-1] << std::endl;
				else if (tagsY[I] == I+(nx-1))
					file << xu[i] << '\t' << yu[j+1] << std::endl;
				file << xu[i] << '\t' << yu[j] << std::endl;
			}
			file << std::endl;
		}
	}
	file.close();

	real *xv = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->xv[0])),
	     *yv = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->yv[0]));
	
	// tag points at which v is evaluated
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			I = j*nx+i + (nx-1)*ny;
			tagsX[I] = -1;
			tagsY[I] = -1;
			coeffsX[I] = 0.0;
			coeffsY[I] = 1.0;
			outsideX = true;
			outsideY = true;
			k = NSWithBody<memoryType>::B.totalPoints - 1;
			l = 0;
			flag = false;
			bdryFlagX = -1;
			bdryFlagY = -1;
			
			while( l < NSWithBody<memoryType>::B.totalPoints)
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
				if (by[bottom]-eps < yv[j] && by[top]-eps > yv[j])
				{
					// if the segment is not parallel to the ray
					if (fabs(by[l]-by[k]) > eps)
					{
						// calculate the point of intersection of the double ray with the boundary
						x = bx[k] + (bx[l]-bx[k]) * (yv[j]-by[k])/(by[l]-by[k]);
						// calculate the body velocity at the point of intersection
						uvX[I] = vB[k] + (vB[l]-vB[k]) * (yv[j]-by[k])/(by[l]-by[k]);

						// if the point of intersection coincides with the grid point
						if (fabs(x-xv[i]) < eps)
						{
							outsideX   = true;
							flag      = true;
							bdryFlagX = I;
							cfX       = 0.0;
						}
						// if the point of intersection lies to the right of the grid point
				 		else if (x > xv[i]+eps)
							outsideX = !outsideX;
						// if the point of intersection is in the cell to the immediate right of the grid point
						if (x>xv[i-1]+eps && x<xv[i]-eps)
						{
							bdryFlagX = I+1;
							cfX = (xv[i]-x)/(xv[i+1]-x);
						}
						// if the point of intersection is in the cell to the immediate left of the grid point
						else if (x>xv[i]+eps && x<xv[i+1]-eps)
						{
							bdryFlagX = I-1;
							cfX = (xv[i]-x)/(xv[i-1]-x);
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
						uvY[I] = vB[k] + (vB[l]-vB[k]) * (xv[i]-bx[k])/(bx[l]-bx[k]);
						// if the point of intersection coincides with the grid point
						if (fabs(y-yv[j]) < eps)
						{
							outsideY   = true;
							bdryFlagY = I;
							cfY       = 0.0;
							flag      = true;
						}
						// if the point of intersection lies to the top of the grid point
				 		else if (y > yv[j]+eps)
							outsideY = !outsideY;

						if (y>yv[j-1]+eps && y<yv[j]-eps)
						{
							bdryFlagY = I+nx;
							cfY = (yv[j]-y)/(yv[j+1]-y);
						}
						else if (y>yv[j]+eps && y<yv[j+1]-eps)
						{
							bdryFlagY = I-nx;
							cfY = (yv[j]-y)/(yv[j-1]-y);
						}
					}
				}
				k = l;
				l = l+1;
			}
			if (outsideY && bdryFlagY>=0)
			{					
				tagsY[I]   = bdryFlagY;
				coeffsY[I] = cfY;
			}
			if (outsideX && bdryFlagX>=0)
			{
				tagsX[I]   = bdryFlagX;
				coeffsX[I] = cfX;
			} 
		}
	}
	
	file.open((folder+"/tagy.txt").c_str());
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			I = j*nx+i+(nx-1)*ny;
			if (tagsY[I] >= 0 || tagsX[I] >= 0)
			{
				file << xv[i] << '\t' << yv[j] << std::endl;
				if (tagsX[I] == I-1)
					file << xv[i-1] << '\t' << yv[j] << std::endl;
				else if (tagsX[I] == I+1)
					file << xv[i+1] << '\t' << yv[j] << std::endl;
				if (tagsY[I] == I-nx)
					file << xv[i] << '\t' << yv[j-1] << std::endl;
				else if (tagsY[I] == I+nx)
					file << xv[i] << '\t' << yv[j+1] << std::endl;
				file << xv[i] << '\t' << yv[j] << std::endl;
			}
			file << std::endl;
		}
	}
	file.close();
	
	file.open((folder+"/body.txt").c_str());
	for(int k=0; k<NSWithBody<memoryType>::B.totalPoints; k++)
	{
		file << bx[k] << '\t' << by[k] << std::endl;
	}
	file.close();
}
#endif
