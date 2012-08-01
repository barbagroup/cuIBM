/**
*  Copyright (C) 2011 by Anush Krishnan, Simon Layton, Lorena Barba
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in
*  all copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*  THE SOFTWARE.
*/

template <>
void FadlunEtAlSolver<host_memory>::tagPoints()
{
	logger.startTimer("tagPoints");

	real *bx = thrust::raw_pointer_cast(&(B.x[0])),
	     *by = thrust::raw_pointer_cast(&(B.y[0]));
	     
	tagPoints(bx, by);

	logger.stopTimer("tagPoints");
}

template <>
void FadlunEtAlSolver<device_memory>::tagPoints()
{
	logger.startTimer("tagPoints");
	
	// transferring boundary point coordinates to the host
	vecH bxH(B.totalPoints), byH(B.totalPoints);
	bxH = B.x;
	byH = B.y;

	// creating raw pointers
	real *bx = thrust::raw_pointer_cast(&(bxH[0])),
	     *by = thrust::raw_pointer_cast(&(byH[0]));
	
	tagPoints(bx, by);
	
	// transferring tag and coeffs data to the device
	tagsXD = tagsX;
	tagsYD = tagsY;
	coeffsXD = coeffsX;
	coeffsYD = coeffsY;
	
	logger.stopTimer("tagPoints");
}

// Bilinear Fadlun1c-type interpolation inside the body.
template <typename memoryType>
void FadlunEtAlSolver<memoryType>::tagPoints(real *bx, real *by)
{
	int  nx = NavierStokesSolver<memoryType>::domInfo->nx,
	     ny = NavierStokesSolver<memoryType>::domInfo->ny;
	
	real *xu = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->xu[0])),
	     *yu = thrust::raw_pointer_cast(&(NavierStokesSolver<memoryType>::domInfo->yu[0]));

	bool insideX, insideY, flag;
	int  bdryFlagX, bdryFlagY, k, l, I;
	int  bottom, top, left, right;
	real eps = 1.e-6;
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
						// if the point of intersection is in the cell to the immediate left of the grid point (was earlier typed as 'right')
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
