/**
* \brief This now works only for bodies that are not inside each other
*/
template <>
void FadlunEtAlSolver<host_memory>::tagPoints()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	real *xu = thrust::raw_pointer_cast(&(domInfo->xu[0])),
	     *yu = thrust::raw_pointer_cast(&(domInfo->yu[0]));
	
	real *bx = thrust::raw_pointer_cast(&(B.x[0])),
	     *by = thrust::raw_pointer_cast(&(B.y[0]));
	
	bool insideX, insideY, flag;
	int  bdryFlagX, bdryFlagY, k, l, I;
	int  bottom, top, left, right;
	real eps = 1.e-6;
	real cf = 0, x, y;

	// tag points at which u is evaluated
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			for(int n=0; n<B.numBodies; n++)
			{
				insideX = false;
				insideY = false;
				k = B.offsets[n+1]-1;
				l = B.offsets[n];
				I = j*(nx-1)+i;
				flag = false;
				tags[I] = -1;
				bdryFlagX = -1;
				bdryFlagY = -1;
				while( l < B.offsets[n+1] && !flag)
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
								flag      = true;
								bdryFlagX = I;
								cf        = 0.0;
							}
							// if the point of intersection lies to the right of the grid point
					 		else if (x > xu[i]+eps)
								insideX = !insideX;
							// if the point of intersection is in the cell to the immediate right of the grid point
							if (x>xu[i-1]+eps && x<xu[i]-eps)
							{
								bdryFlagX = I-1;
								cf = (xu[i]-x)/(xu[i+1]-x);
							}
							// if the point of intersection is in the cell to the immediate left of the grid point
							else if (x>xu[i]+eps && x<xu[i+1]-eps)
							{
								bdryFlagX = I+1;
								cf = (xu[i]-x)/(xu[i-1]-x);
							}
						}
					}
					// consider rays along the y-direction
					if (bx[left]-eps < xu[i] && bx[right]-eps > xu[i])
					{
						if (fabs(bx[l]-bx[k]) > eps)
						{
							y = by[k] + (by[l]-by[k]) * (xu[i]-bx[k]) / (bx[l]-bx[k]);
							if (fabs(y-yu[j]) < eps)
							{
								insideY   = true;
								flag      = true;
								bdryFlagY = I;
								cf        = 0.0;
							}
					 		else if (y > yu[j]+eps)
								insideY = !insideY;

							if (y>yu[j-1]+eps && y<yu[j]-eps)
							{
								bdryFlagY = I-(nx-1);
								cf = (yu[j]-y)/(yu[j+1]-y);
							}
							else if (y>yu[j]+eps && y<yu[j+1]-eps)
							{
								bdryFlagY = I+(nx-1);
								cf = (yu[j]-y)/(yu[j-1]-y);
							}
						}
					}
					k = l;
					l = l+1;
				}
				if (insideX && bdryFlagX>=0)
				{
					tags[I]   = bdryFlagX;
					coeffs[I] = cf;
				}
				else if (insideY && bdryFlagY>=0)
				{					
					tags[I]   = bdryFlagY;
					coeffs[I] = cf;
				}
			}
		}
	}
	
	std::ofstream file("tagged.txt");
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			I = j*(nx-1)+i;
			if (tags[I] >= 0)
				file << xu[i] << '\t' << yu[i] << std::endl;
		}
	}
	file.close();
	
	real *xv = thrust::raw_pointer_cast(&(domInfo->xv[0])),
	     *yv = thrust::raw_pointer_cast(&(domInfo->yv[0]));
	
	// tag points at which v is evaluated
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			for(int n=0; n<B.numBodies; n++)
			{
				insideX = false;
				insideY = false;
				k = B.offsets[n+1]-1;
				l = B.offsets[n];
				I = j*nx+i + (nx-1)*ny;
				flag = false;
				tags[I] = -1;
				bdryFlagX = -1;
				bdryFlagY = -1;
				while( l < B.offsets[n+1] && !flag)
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
								cf        = 0.0;
							}
							// if the point of intersection lies to the right of the grid point
					 		else if (x > xv[i]+eps)
								insideX = !insideX;
							// if the point of intersection is in the cell to the immediate right of the grid point
							if (x>xv[i-1]+eps && x<xv[i]-eps)
							{
								bdryFlagX = I-1;
								cf = (xv[i]-x)/(xv[i+1]-x);
							}
							// if the point of intersection is in the cell to the immediate left of the grid point
							else if (x>xv[i]+eps && x<xv[i+1]-eps)
							{
								bdryFlagX = I+1;
								cf = (xv[i]-x)/(xv[i-1]-x);
							}
						}
					}
					// consider rays along the y-direction
					if (bx[left]-eps < xv[i] && bx[right]-eps > xv[i])
					{
						if (fabs(bx[l]-bx[k]) > eps)
						{
							y = by[k] + (by[l]-by[k]) * (xv[i]-bx[k]) / (bx[l]-bx[k]);
							if (fabs(y-yv[j]) < eps)
							{
								insideY   = true;
								flag      = true;
								bdryFlagY = I;
								cf        = 0.0;
							}
					 		else if (y > yv[j]+eps)
								insideY = !insideY;

							if (y>yv[j-1]+eps && y<yv[j]-eps)
							{
								bdryFlagY = I-(nx-1);
								cf = (yv[j]-y)/(yv[j+1]-y);
							}
							else if (y>yv[j]+eps && y<yv[j+1]-eps)
							{
								bdryFlagY = I+(nx-1);
								cf = (yv[j]-y)/(yv[j-1]-y);
							}
						}
					}
					k = l;
					l = l+1;
				}
				if (insideX && bdryFlagX>=0)
				{
					tags[I]   = bdryFlagX;
					coeffs[I] = cf;
				}
				else if (insideY && bdryFlagY>=0)
				{					
					tags[I]   = bdryFlagY;
					coeffs[I] = cf;
				}
			}
		}
	}
}

/**
* \brief This now works only for bodies that are not inside each other
*/
template <>
void FadlunEtAlSolver<device_memory>::tagPoints()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	std::cout << B.totalPoints << "\t" << B.numBodies << std::endl;
	vecH bx(B.totalPoints), by(B.totalPoints);
	bx = B.x;
	by = B.y;
	
	cusp::array1d<int, host_memory>
	     numPoints(B.numBodies), offsets(B.numBodies);

	offsets   = B.offsets;
	numPoints = B.numPoints;
	std::cout << "Transferred" << std::endl;
	
	real *xu = thrust::raw_pointer_cast(&(domInfo->xu[0])),
	     *yu = thrust::raw_pointer_cast(&(domInfo->yu[0]));
	
	bool insideX, insideY, flag;
	int  bdryFlagX, bdryFlagY, k, l, I;
	int  bottom, top, left, right;
	real eps = 1.e-6;
	real cf = 0, x, y;

	std::cout << "Tagging points... ";
	// tag points at which u is evaluated
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			I = j*(nx-1)+i;
			tags[I] = -1;
			insideX = false;
			insideY = false;
			k = B.totalPoints - 1;
			l = 0;
			flag = false;
			bdryFlagX = -1;
			bdryFlagY = -1;
				
			while( l < B.totalPoints && !flag)
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
							flag      = true;
							bdryFlagX = I;
							cf        = 0.0;
						}
						// if the point of intersection lies to the right of the grid point
				 		else if (x > xu[i]+eps)
							insideX = !insideX;
						// if the point of intersection is in the cell to the immediate right of the grid point
						if (x>xu[i-1]+eps && x<xu[i]-eps)
						{
							bdryFlagX = I-1;
							cf = (xu[i]-x)/(xu[i+1]-x);
						}
						// if the point of intersection is in the cell to the immediate left of the grid point
						else if (x>xu[i]+eps && x<xu[i+1]-eps)
						{
							bdryFlagX = I+1;
							cf = (xu[i]-x)/(xu[i-1]-x);
						}
					}
				}
				// consider rays along the y-direction
				if (bx[left]-eps < xu[i] && bx[right]-eps > xu[i])
				{
					if (fabs(bx[l]-bx[k]) > eps)
					{
						y = by[k] + (by[l]-by[k]) * (xu[i]-bx[k]) / (bx[l]-bx[k]);
						if (fabs(y-yu[j]) < eps)
						{
							insideY   = true;
							flag      = true;
							bdryFlagY = I;
							cf        = 0.0;
						}
				 		else if (y > yu[j]+eps)
							insideY = !insideY;

						if (y>yu[j-1]+eps && y<yu[j]-eps)
						{
							bdryFlagY = I-(nx-1);
							cf = (yu[j]-y)/(yu[j+1]-y);
						}
						else if (y>yu[j]+eps && y<yu[j+1]-eps)
						{
							bdryFlagY = I+(nx-1);
							cf = (yu[j]-y)/(yu[j-1]-y);
						}
					}
				}
				k = l;
				l = l+1;
			}
			if (insideX && bdryFlagX>=0)
			{
				tags[I]   = bdryFlagX;
				coeffs[I] = cf;
			}
			else if (insideY && bdryFlagY>=0)
			{					
				tags[I]   = bdryFlagY;
				coeffs[I] = cf;
			}
		}
	}
	
	std::ofstream file("tagx.txt");
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			I = j*(nx-1)+i;
			if (tags[I] >= 0)
				file << xu[i] << '\t' << yu[j] << std::endl;
		}
	}
	file.close();

	real *xv = thrust::raw_pointer_cast(&(domInfo->xv[0])),
	     *yv = thrust::raw_pointer_cast(&(domInfo->yv[0]));
	
	// tag points at which v is evaluated
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			I = j*nx+i + (nx-1)*ny;
			tags[I] = -1;
			insideX = false;
			insideY = false;
			k = B.totalPoints-1;
			l = 0;
			flag = false;
			bdryFlagX = -1;
			bdryFlagY = -1;
			while( l < B.totalPoints && !flag)
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
							cf        = 0.0;
						}
						// if the point of intersection lies to the right of the grid point
				 		else if (x > xv[i]+eps)
							insideX = !insideX;
						// if the point of intersection is in the cell to the immediate right of the grid point
						if (x>xv[i-1]+eps && x<xv[i]-eps)
						{
							bdryFlagX = I-1;
							cf = (xv[i]-x)/(xv[i+1]-x);
						}
						// if the point of intersection is in the cell to the immediate left of the grid point
						else if (x>xv[i]+eps && x<xv[i+1]-eps)
						{
							bdryFlagX = I+1;
							cf = (xv[i]-x)/(xv[i-1]-x);
						}
					}
				}
				// consider rays along the y-direction
				if (bx[left]-eps < xv[i] && bx[right]-eps > xv[i])
				{
					if (fabs(bx[l]-bx[k]) > eps)
					{
						y = by[k] + (by[l]-by[k]) * (xv[i]-bx[k]) / (bx[l]-bx[k]);
						if (fabs(y-yv[j]) < eps)
						{
							insideY   = true;
							flag      = true;
							bdryFlagY = I;
							cf        = 0.0;
						}
				 		else if (y > yv[j]+eps)
							insideY = !insideY;

						if (y>yv[j-1]+eps && y<yv[j]-eps)
						{
							bdryFlagY = I-nx;
							cf = (yv[j]-y)/(yv[j+1]-y);
						}
						else if (y>yv[j]+eps && y<yv[j+1]-eps)
						{
							bdryFlagY = I+nx;
							cf = (yv[j]-y)/(yv[j-1]-y);
						}
					}
				}
				k = l;
				l = l+1;
			}
			if (insideX && bdryFlagX>=0)
			{
				tags[I]   = bdryFlagX;
				coeffs[I] = cf;
			}
			else if (insideY && bdryFlagY>=0)
			{					
				tags[I]   = bdryFlagY;
				coeffs[I] = cf;
			}
		}
	}
	std::cout << "DONE!" << std::endl;
	
	file.open("tagy.txt");
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			I = j*nx+i+(nx-1)*ny;
			if (tags[I] >= 0)
				file << xv[i] << '\t' << yv[j] << std::endl;
		}
	}
	file.close();
	
	file.open("body.txt");
	for(int k=0; k<B.totalPoints; k++)
	{
		file << bx[k] << '\t' << by[k] << std::endl;
	}
	file.close();
	
	tagsD = tags;
	coeffsD = coeffs;
}

#if 0
template <>
void FadlunEtAlSolver<device_memory>::tagPoints()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	vecH bx(B.totalPoints), by(B.totalPoints);
	bx = B.x;
	by = B.y;
	
	cusp::array1d<int, host_memory>
	     numPoints(B.numBodies), offsets(B.numBodies);
	
	offsets   = B.offsets;
	numPoints = B.numPoints;
	
	real *xu = thrust::raw_pointer_cast(&(domInfo->xu[0])),
	     *yu = thrust::raw_pointer_cast(&(domInfo->yu[0]));
	
	bool insideX, insideY, flag;
	int  bdryFlagX, bdryFlagY, k, l, I;
	int  bottom, top, left, right;
	real eps = 1.e-6;
	real cf = 0, x, y;

	// tag points at which u is evaluated
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			I = j*(nx-1)+i;
			tags[I] = -1;
			for(int n=0; n<B.numBodies; n++)
			{
				insideX = false;
				insideY = false;
				k = offsets[n] + numPoints[n] - 1;
				l = offsets[n];
				flag = false;
				bdryFlagX = -1;
				bdryFlagY = -1;
				while( l < offsets[n]+numPoints[n] && !flag)
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
								flag      = true;
								bdryFlagX = I;
								cf        = 0.0;
							}
							// if the point of intersection lies to the right of the grid point
					 		else if (x > xu[i]+eps)
								insideX = !insideX;
							// if the point of intersection is in the cell to the immediate right of the grid point
							if (x>xu[i-1]+eps && x<xu[i]-eps)
							{
								bdryFlagX = I-1;
								cf = (xu[i]-x)/(xu[i+1]-x);
							}
							// if the point of intersection is in the cell to the immediate left of the grid point
							else if (x>xu[i]+eps && x<xu[i+1]-eps)
							{
								bdryFlagX = I+1;
								cf = (xu[i]-x)/(xu[i-1]-x);
							}
						}
					}
					// consider rays along the y-direction
					if (bx[left]-eps < xu[i] && bx[right]-eps > xu[i])
					{
						if (fabs(bx[l]-bx[k]) > eps)
						{
							y = by[k] + (by[l]-by[k]) * (xu[i]-bx[k]) / (bx[l]-bx[k]);
							if (fabs(y-yu[j]) < eps)
							{
								insideY   = true;
								flag      = true;
								bdryFlagY = I;
								cf        = 0.0;
							}
					 		else if (y > yu[j]+eps)
								insideY = !insideY;

							if (y>yu[j-1]+eps && y<yu[j]-eps)
							{
								bdryFlagY = I-(nx-1);
								cf = (yu[j]-y)/(yu[j+1]-y);
							}
							else if (y>yu[j]+eps && y<yu[j+1]-eps)
							{
								bdryFlagY = I+(nx-1);
								cf = (yu[j]-y)/(yu[j-1]-y);
							}
						}
					}
					k = l;
					l = l+1;
				}
				if (insideX && bdryFlagX>=0)
				{
					tags[I]   = bdryFlagX;
					coeffs[I] = cf;
				}
				else if (insideY && bdryFlagY>=0)
				{					
					tags[I]   = bdryFlagY;
					coeffs[I] = cf;
				}
			}
		}
	}
	
	std::ofstream file("tagx.txt");
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			I = j*(nx-1)+i;
			if (tags[I] >= 0)
				file << xu[i] << '\t' << yu[j] << std::endl;
		}
	}
	file.close();

	real *xv = thrust::raw_pointer_cast(&(domInfo->xv[0])),
	     *yv = thrust::raw_pointer_cast(&(domInfo->yv[0]));
	
	// tag points at which v is evaluated
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			I = j*nx+i + (nx-1)*ny;
			tags[I] = -1;
			for(int n=0; n<B.numBodies; n++)
			{
				insideX = false;
				insideY = false;
				k = offsets[n]+numPoints[n]-1;
				l = offsets[n];
				flag = false;
				bdryFlagX = -1;
				bdryFlagY = -1;
				while( l < offsets[n]+numPoints[n] && !flag)
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
								cf        = 0.0;
							}
							// if the point of intersection lies to the right of the grid point
					 		else if (x > xv[i]+eps)
								insideX = !insideX;
							// if the point of intersection is in the cell to the immediate right of the grid point
							if (x>xv[i-1]+eps && x<xv[i]-eps)
							{
								bdryFlagX = I-1;
								cf = (xv[i]-x)/(xv[i+1]-x);
							}
							// if the point of intersection is in the cell to the immediate left of the grid point
							else if (x>xv[i]+eps && x<xv[i+1]-eps)
							{
								bdryFlagX = I+1;
								cf = (xv[i]-x)/(xv[i-1]-x);
							}
						}
					}
					// consider rays along the y-direction
					if (bx[left]-eps < xv[i] && bx[right]-eps > xv[i])
					{
						if (fabs(bx[l]-bx[k]) > eps)
						{
							y = by[k] + (by[l]-by[k]) * (xv[i]-bx[k]) / (bx[l]-bx[k]);
							if (fabs(y-yv[j]) < eps)
							{
								insideY   = true;
								flag      = true;
								bdryFlagY = I;
								cf        = 0.0;
							}
					 		else if (y > yv[j]+eps)
								insideY = !insideY;

							if (y>yv[j-1]+eps && y<yv[j]-eps)
							{
								bdryFlagY = I-nx;
								cf = (yv[j]-y)/(yv[j+1]-y);
							}
							else if (y>yv[j]+eps && y<yv[j+1]-eps)
							{
								bdryFlagY = I+nx;
								cf = (yv[j]-y)/(yv[j-1]-y);
							}
						}
					}
					k = l;
					l = l+1;
				}
				if (insideX && bdryFlagX>=0)
				{
					tags[I]   = bdryFlagX;
					coeffs[I] = cf;
				}
				else if (insideY && bdryFlagY>=0)
				{					
					tags[I]   = bdryFlagY;
					coeffs[I] = cf;
				}
			}
		}
	}
	
	file.open("tagy.txt");
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			I = j*nx+i+(nx-1)*ny;
			if (tags[I] >= 0)
				file << xv[i] << '\t' << yv[j] << std::endl;
		}
	}
	file.close();
	
	tagsD = tags;
	coeffsD = coeffs;
}
#endif
