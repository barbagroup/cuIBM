template <typename Matrix, typename Vector>
void TairaColoniusSolver<Matrix, Vector>::generateBC2()
{
	int  nx = NavierStokesSolver<Matrix, Vector>::domInfo->nx,
	     ny = NavierStokesSolver<Matrix, Vector>::domInfo->ny;
	//int  numP  = nx*ny;
	
	real *dx    = thrust::raw_pointer_cast(&(NavierStokesSolver<Matrix, Vector>::domInfo->dx[0])),
	     *dy    = thrust::raw_pointer_cast(&(NavierStokesSolver<Matrix, Vector>::domInfo->dy[0]));
	
	/// rhs2 of size np+2nb
	
	/// boundary conditions from the continuity equations
		
	cusp::blas::fill(NavierStokesSolver<Matrix, Vector>::bc2, 0.0);
	
	for(int i=0; i<nx; i++)
	{
		NavierStokesSolver<Matrix, Vector>::bc2[i] -= NavierStokesSolver<Matrix, Vector>::bcN[YMINUS][i+nx-1]*dx[i]; // v[0][i+1]*dx;
		NavierStokesSolver<Matrix, Vector>::bc2[(ny-1)*nx + i] += NavierStokesSolver<Matrix, Vector>::bcN[YPLUS][i+nx-1]*dx[i]; // v[ny][i+1]*dx;
	}

	for(int j=0; j<ny; j++)
	{
		// left
		NavierStokesSolver<Matrix, Vector>::bc2[j*nx] -= NavierStokesSolver<Matrix, Vector>::bcN[XMINUS][j]*dy[j]; // u[j+1][0]*dy;
		// right
		NavierStokesSolver<Matrix, Vector>::bc2[j*nx+nx-1] += NavierStokesSolver<Matrix, Vector>::bcN[XPLUS][j]*dy[j]; // u[j+1][nx]*dy;
	}	

	int	row = nx*ny;
	/// no-slip condition on the body surface
	for(int k=0; k<B.totalPoints; k++)
	{
		NavierStokesSolver<Matrix, Vector>::bc2[row] = B.uB[k];
		row++;
	}
	for(int k=0; k<B.totalPoints; k++)
	{
		NavierStokesSolver<Matrix, Vector>::bc2[row] = B.vB[k];
		row++;
	}
}