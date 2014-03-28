/*
*  Copyright (C) 2012 by Anush Krishnan, Simon Layton, Lorena Barba
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

#include <solvers/NavierStokes/kernels/calculateForce.h>

template <typename memoryType>
void TairaColoniusSolver<memoryType>::calculateForce()
{
	typedef typename cusp::coo_matrix<int, real, memoryType>   Matrix;
	typedef typename cusp::array1d<int, memoryType>::iterator  IndexIterator;
	typedef typename cusp::array1d<real, memoryType>::iterator ValueIterator;
	typedef typename cusp::array1d_view<IndexIterator>         IndexView;
	typedef typename cusp::array1d_view<ValueIterator>         VecView;
	typedef typename cusp::coo_matrix_view<IndexView, IndexView, VecView> MatView;
	
	int  nx = NavierStokesSolver<memoryType>::domInfo->nx,
	     ny = NavierStokesSolver<memoryType>::domInfo->ny,
	     ETRows = (nx-1)*ny + nx*(ny-1), 
	     ETCols,
	     numBodies = NSWithBody<memoryType>::B.numBodies,
	     totalPoints = NSWithBody<memoryType>::B.totalPoints;
	
	cusp::array1d<real, memoryType> F(ETRows);
	
	VecView  fView;
	MatView  ETView;
	real     dx, dy;
	
	// loop through bodies
	for(int l=0; l < numBodies; l++)
	{
		dx = NavierStokesSolver<memoryType>::domInfo->dx[ NSWithBody<memoryType>::B.I[l] ],
		dy = dx;
	    
		// x-component of the force
		if(l < numBodies-1)
		{
			fView  = VecView(NavierStokesSolver<memoryType>::lambda.begin() + nx*ny + NSWithBody<memoryType>::B.offsets[l],
			                 NavierStokesSolver<memoryType>::lambda.begin() + nx*ny + NSWithBody<memoryType>::B.offsets[l+1]
			                );
			ETCols = NSWithBody<memoryType>::B.offsets[l+1] - NSWithBody<memoryType>::B.offsets[l];
			ETView = MatView(ETRows, ETCols, 12*ETCols,
			                 cusp::make_array1d_view(ET.row_indices.begin(), ET.row_indices.end()),
			                 cusp::make_array1d_view(ET.column_indices.begin() + NSWithBody<memoryType>::B.offsets[l], ET.column_indices.begin() + NSWithBody<memoryType>::B.offsets[l+1]),
			                 cusp::make_array1d_view(ET.values.begin() + NSWithBody<memoryType>::B.offsets[l], ET.values.begin() + NSWithBody<memoryType>::B.offsets[l+1])
			                );
		}
		else
		{
			fView  = VecView(NavierStokesSolver<memoryType>::lambda.begin() + nx*ny + NSWithBody<memoryType>::B.offsets[l],
			                 NavierStokesSolver<memoryType>::lambda.begin() + nx*ny + totalPoints
			                );
			ETCols = totalPoints - NSWithBody<memoryType>::B.offsets[l];
			ETView = MatView(ETRows, ETCols, 12*ETCols,
			                 cusp::make_array1d_view(ET.row_indices.begin(), ET.row_indices.end()),
			                 cusp::make_array1d_view(ET.column_indices.begin() + NSWithBody<memoryType>::B.offsets[l], ET.column_indices.begin() + totalPoints),
			                 cusp::make_array1d_view(ET.values.begin() + NSWithBody<memoryType>::B.offsets[l], ET.values.begin() + totalPoints)
			                );
		}
		
		cusp::multiply(ETView, fView, F);
		NSWithBody<memoryType>::B.forceX[l] = (dx*dy)/dx*thrust::reduce( F.begin(), F.begin()+(nx-1)*ny );
		
		// y-component of the force
		if(l < numBodies-1)
		{
			fView  = VecView(NavierStokesSolver<memoryType>::lambda.begin() + nx*ny + totalPoints + NSWithBody<memoryType>::B.offsets[l], 
			                 NavierStokesSolver<memoryType>::lambda.begin() + nx*ny + totalPoints + NSWithBody<memoryType>::B.offsets[l+1]
			                );
			ETCols = NSWithBody<memoryType>::B.offsets[l+1] - NSWithBody<memoryType>::B.offsets[l];
			ETView = MatView(ETRows, ETCols, 12*ETCols,
			                 cusp::make_array1d_view(ET.row_indices.begin(), ET.row_indices.end()),
			                 cusp::make_array1d_view(ET.column_indices.begin() + totalPoints + NSWithBody<memoryType>::B.offsets[l], ET.column_indices.begin() + totalPoints + NSWithBody<memoryType>::B.offsets[l+1]),
			                 cusp::make_array1d_view(ET.values.begin() + totalPoints + NSWithBody<memoryType>::B.offsets[l], ET.values.begin() + totalPoints + NSWithBody<memoryType>::B.offsets[l+1])
			                );
		}
		else
		{
			fView  = VecView(NavierStokesSolver<memoryType>::lambda.begin() + nx*ny + totalPoints + NSWithBody<memoryType>::B.offsets[l],
			                 NavierStokesSolver<memoryType>::lambda.end()
			                );
			ETCols = totalPoints - NSWithBody<memoryType>::B.offsets[l];
			ETView = MatView(ETRows, ETCols, 12*ETCols,
			                 cusp::make_array1d_view(ET.row_indices.begin(), ET.row_indices.end()),
			                 cusp::make_array1d_view(ET.column_indices.begin() + totalPoints + NSWithBody<memoryType>::B.offsets[l], ET.column_indices.end()),
			                 cusp::make_array1d_view(ET.values.begin() + totalPoints + NSWithBody<memoryType>::B.offsets[l], ET.values.end())
			                );
		}
		
		cusp::multiply(ETView, fView, F);
		NSWithBody<memoryType>::B.forceY[l] = (dx*dy)/dy*thrust::reduce( F.begin()+(nx-1)*ny, F.end() );
	}
}
