#include <solvers/NavierStokes/NavierStokesSolver.h>

template <typename Matrix, typename Vector>
void NavierStokesSolver<Matrix, Vector>::generateBC2()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	//int  numP  = nx*ny;
	
	real *dx    = thrust::raw_pointer_cast(&(domInfo->dx[0])),
	     *dy    = thrust::raw_pointer_cast(&(domInfo->dy[0]));
	
	/// rhs2 of size np+2nb
	
	/// boundary conditions from the continuity equations
		
	cusp::blas::fill(bc2, 0.0);
	
	for(int i=0; i<nx; i++)
	{
		bc2[i] -= bcN[YMINUS][i+nx-1]*dx[i]; // v[0][i+1]*dx;
		bc2[(ny-1)*nx + i] += bcN[YPLUS][i+nx-1]*dx[i]; // v[ny][i+1]*dx;
	}

	for(int j=0; j<ny; j++)
	{
		// left
		bc2[j*nx] -= bcN[XMINUS][j]*dy[j]; // u[j+1][0]*dy;
		// right
		bc2[j*nx+nx-1] += bcN[XPLUS][j]*dy[j]; // u[j+1][nx]*dy;
	}	
/*	int	row = np;
	
	/// no-slip condition on the body surface
	for(k=0; k<B.no_of_bodies; k++)
	{
		for(int l=0; l<B.Body[k].no_of_points; l++)
		{
			rhs2[row] = B.Body[k].uB[l];
			row++;
		}
	}
	for(k=0; k<B.no_of_bodies; k++)
	{
		for(int l=0; l<B.Body[k].no_of_points; l++)
		{
			rhs2[row] = B.Body[k].vB[l];
			row++;
		}
	}*/
}