/***************************************************************************//**
* \file generateRN.cu
* \author Krishnan, A. (anush@bu.edu)
* \brief Definition of the kernels required to generate vector rn
*/

#include <solvers/NavierStokes/kernels/generateRN.h>

#define BSZ 16

/********************//**
* \namespace kernels
* \brief Contain all custom-written CUDA kernels
*/
namespace kernels
{

/**
\brief To be documented
*/
// CUDA kernel to generate Hx and rn
__global__ void convectionTermU(real *rn, real *H, real *q,
                                  int nx, int ny, real *dx, real *dy, 
                                  real dt, real gamma, real zeta, real alpha, real nu)
{
	int bx	= blockIdx.x,
		by	= blockIdx.y,
		i	= threadIdx.x,
		j	= threadIdx.y;
	
	// work out global index of first point in block
	//int I = (BSZ-2)*bx + i,
	//	J = (BSZ-2)*by + j;
	int I = (BSZ-2)*bx + i,
	    J = (BSZ-2)*by + j;
	
	if (I >= nx-1 || J >= ny) {
		return;
	}

	int  N_u = (nx-1)*ny,
	     Gidx_x = J*(nx-1) + I,
	     Gidx_y = (J-1)*nx + I + N_u;
			
	real cTerm, dTerm, Hxn;
	
	__shared__ real u[BSZ][BSZ],
	                  v[BSZ][BSZ],
	                  Hx[BSZ][BSZ];
						
	__shared__ real		Dx[BSZ][BSZ], Dy[BSZ][BSZ];
	
	Dy[j][i] = dy[J];
	Dx[j][i] = dx[I];
	
	/// transfer from global to shared memory
	u[j][i] = q[Gidx_x]/Dy[j][i];
	v[j][i] = q[Gidx_y]/Dx[j][i];
	__syncthreads();
	
	/// check bounds for convective term in the x-direction
	int		global_check = ( I==0 || I==(nx-2) || J==0 || J==(ny-1) ),		///< check if we compute globally
			block_check  = ( i==0 || i==(BSZ-1) || j==0 || j==(BSZ-1) );	///< check if element within block computes
	
	/// X-component
	if( !(global_check || block_check) )
	{
		// copy Hx -> Hxn, Hy -> Hyn
		Hxn = H[Gidx_x];
		
		// apply stencil
		Hx[j][i] = - ( (u[j][i+1]+u[j][i])*(u[j][i+1]+u[j][i]) - (u[j][i]+u[j][i-1])*(u[j][i]+u[j][i-1]) )/( 2.0 * (Dx[j][i]+Dx[j][i+1]) ) \
						- ( (u[j+1][i]+u[j][i])*(v[j+1][i]+v[j+1][i+1]) - (u[j][i]+u[j-1][i])*(v[j][i]+v[j][i+1]) )/( 4.0 * Dy[j][i] );

		H[Gidx_x] = Hx[j][i];
		
		// rN for u
		cTerm = gamma*Hx[j][i] + zeta*Hxn;
		//cTerm = Hx[j][i]; // 1st order Euler
		dTerm = alpha*nu*2.0*( \
						 ( Dx[j][i]*u[j][i+1] - (Dx[j][i]+Dx[j][i+1])*u[j][i] + Dx[j][i+1]*u[j][i-1] )/( Dx[j][i]*Dx[j][i+1]*(Dx[j][i]+Dx[j][i+1]) ) \
					   + 4.0*( (Dy[j][i]+Dy[j-1][i])*u[j+1][i] - (Dy[j-1][i] + 2.0*Dy[j][i] + Dy[j+1][i])*u[j][i] + (Dy[j][i]+Dy[j+1][i])*u[j-1][i] ) \
							/( (Dy[j][i]+Dy[j-1][i]) * (Dy[j-1][i] + 2.0*Dy[j][i] + Dy[j+1][i]) * (Dy[j][i]+Dy[j+1][i]) ) \
					     );
		rn[Gidx_x] = (u[j][i]/dt + cTerm + dTerm) * 0.5*(Dx[j][i]+Dx[j][i+1]);
	}
}

/**
* \brief To be documented
*/
// CUDA kernel to generate Hy and rn
__global__ void convectionTermV(real *rn, real *H, real *q,
                                  int nx, int ny, real *dx, real *dy,
                                  real dt, real gamma, real zeta, real alpha, real nu)
{
	int bx	= blockIdx.x,
		by	= blockIdx.y,
		i	= threadIdx.x,
		j	= threadIdx.y;
	
	// work out global index of first point in block
	int I = (BSZ-2)*bx + i,
		J = (BSZ-2)*by + j;
	
	if (I >= nx || J >= ny-1) {
		return;
	}

	int		N_u		= (nx-1)*ny,
			Gidx_x	= J*(nx-1) + I,
			Gidx_y	= J*nx + I + N_u;
			
	real	cTerm, dTerm, Hyn;
	
	__shared__ real		u[BSZ][BSZ],
						v[BSZ][BSZ],
						Hy[BSZ][BSZ];
						
	__shared__ real		Dx[BSZ][BSZ], Dy[BSZ][BSZ];
	
	Dy[j][i] = dy[J];
	Dx[j][i] = dx[I];
	
	/// transfer from global to shared memory
	u[j][i] = q[Gidx_x]/Dy[j][i];
	v[j][i] = q[Gidx_y]/Dx[j][i];
	__syncthreads();
	
	/// check bounds for convective term in the x-direction
	int		global_check = ( I==0 || I==(nx-1) || J==0 || J==(ny-2) ),		///< check if we compute globally
			block_check  = ( i==0 || i==(BSZ-1) || j==0 || j==(BSZ-1) );	///< check if element within block computes
	
	/// Y-component
	if( !(global_check || block_check) )
	{	
		// Y-component
		// copy global data to the register (to store previous value)
		Hyn = H[Gidx_y];
		
		// apply stencil
		Hy[j][i] = - ( (u[j+1][i]+u[j][i])*(v[j][i]+v[j][i+1]) - (u[j+1][i-1]+u[j][i-1])*(v[j][i]+v[j][i-1]) )/(4.0 * Dx[j][i]) \
						- ( (v[j+1][i]+v[j][i])*(v[j+1][i]+v[j][i]) - (v[j-1][i]+v[j][i])*(v[j-1][i]+v[j][i]) )/( 2.0 * (Dy[j+1][i]+Dy[j][i]) );

		// store newly calculated value in global memory (to be used in the next time step)
		H[Gidx_y] = Hy[j][i];
		
		// rN for v
		cTerm = gamma*Hy[j][i] + zeta*Hyn;
		//cTerm = Hy[j][i]; // 1st order Euler
		dTerm = alpha*nu*2.0*( \
						4.0*( (Dx[j][i-1]+Dx[j][i])*v[j][i+1] - (Dx[j][i-1]+2.0*Dx[j][i]+Dx[j][i+1])*v[j][i] + (Dx[j][i]+Dx[j][i+1])*v[j][i-1] ) \
							/( (Dx[j][i-1]+Dx[j][i]) * (Dx[j][i]+Dx[j][i+1]) * (Dx[j][i-1]+2.0*Dx[j][i]+Dx[j][i+1]) ) \
						+ ( Dy[j][i]*v[j+1][i] - (Dy[j][i]+Dy[j+1][i])*v[j][i] + Dy[j+1][i]*v[j-1][i] )/( Dy[j][i]*Dy[j+1][i]*(Dy[j][i]+Dy[j+1][i]) ) \
					);
		rn[Gidx_y] = (v[j][i]/dt + cTerm + dTerm) * 0.5*(Dy[j][i]+Dy[j+1][i]);
	}
}

/**
* \brief To be documented
*/
__global__
void convectionTermVBottomTop(real *rn, real *H, real *q, \
                                      int nx, int ny, real *dx, real *dy, \
                                      real dt, real gamma, real zeta, real alpha, real nu,\
                                      real *bcBottom, real *bcTop)
{
	int	I = blockIdx.x*blockDim.x + threadIdx.x;
	
	if(I==0 || I >= nx-1) return;
	/// Boundary Conditions for v at the Bottom ***********************************************************************************
	int		Iu, Iv = I + (nx-1)*ny;	
	real	Hn, cTerm, dTerm;
	/// Convection Term
	Hn = H[Iv];
	H[Iv] = -( \
				( 0.5*(q[Iv]/dx[I]+q[Iv+1]/dx[I+1]) ) * ( 0.5*(q[I]/dy[0]+q[I+(nx-1)]/dy[1]) ) \
				- ( 0.5*(q[Iv-1]/dx[I-1]+q[Iv]/dx[I]) ) * ( 0.5*(q[I-1]/dy[0]+q[I-1+(nx-1)]/dy[1]) ) \
			 )/dx[I] \
	 		-( \
				(0.5*(q[Iv] + q[Iv+nx])/dx[I]) * (0.5*(q[Iv] + q[Iv+nx])/dx[I]) \
				- (0.5*(bcBottom[I+nx-1] + q[Iv]/dx[I])) * (0.5*(bcBottom[I+nx-1] + q[Iv]/dx[I])) \
			 )/(0.5*(dy[0] + dy[1]));
	cTerm = gamma*H[Iv] + zeta*Hn;	/// 2nd order Adams-Bashforth
	//cTerm = H[Iv];	/// 1st order Euler
	/// Diffusion Term
	dTerm = alpha*nu*2.0*( \
					4.0 * ( (dx[I-1]+dx[I])*q[Iv+1]/dx[I+1] - (dx[I-1]+2.0*dx[I]+dx[I+1])*q[Iv]/dx[I] + (dx[I]+dx[I+1])*q[Iv-1]/dx[I-1] ) \
						/( (dx[I-1]+dx[I]) * (dx[I]+dx[I+1]) * (dx[I-1]+2.0*dx[I]+dx[I+1]) ) \
					+ ( dy[0]*q[Iv+nx]/dx[I] - (dy[0]+dy[1])*q[Iv]/dx[I] + dy[1]*bcBottom[I+nx-1] )/( dy[0]*dy[1]*(dy[0]+dy[1]) ) \
				);
	/// Calculate rn
	rn[Iv] = ( q[Iv]/(dx[I]*dt) + cTerm + dTerm ) * 0.5*(dy[0] + dy[1]);
	/// Boundary conditions for v at the Top **************************************************************************************
	Iu = I + (ny-2)*(nx-1);
	Iv = I + (nx-1)*ny + (ny-2)*nx;
	/// Convection Term
	Hn = H[Iv];
	H[Iv] = -(
				( 0.5*(q[Iv]/dx[I] + q[Iv+1]/dx[I+1]) )*( 0.5*(q[Iu]/dy[ny-2] + q[Iu+(nx-1)]/dy[ny-1]) ) \
				- ( 0.5*(q[Iv-1]/dx[I-1] + q[Iv]/dx[I]) )*( 0.5*(q[Iu-1]/dy[ny-2] + q[Iu-1+(nx-1)]/dy[ny-1]) )
			 )/(dx[I])
			-(
				( 0.5*(q[Iv]/dx[I] + bcTop[I+nx-1]) )*( 0.5*(q[Iv]/dx[I] + bcTop[I+nx-1]) ) \
				- ( 0.5*(q[Iv-nx] + q[Iv])/dx[I] )*( 0.5*(q[Iv-nx] + q[Iv])/dx[I] )
			 )/(0.5*(dy[ny-2] + dy[ny-1]));
	cTerm = gamma*H[Iv] + zeta*Hn;	/// 2nd order Adams-Bashforth
	//cTerm = H[Iv];	/// 1st order Euler
	/// Diffusion Term
	dTerm = alpha*nu*2.0*( \
					4.0*( (dx[I-1]+dx[I])*q[Iv+1]/dx[I+1] - (dx[I-1]+2.0*dx[I]+dx[I+1])*q[Iv]/dx[I] + (dx[I]+dx[I+1])*q[Iv-1]/dx[I-1] ) \
						/( (dx[I-1]+dx[I]) * (dx[I]+dx[I+1]) * (dx[I-1]+2.0*dx[I]+dx[I+1]) ) \
					+ ( dy[ny-2]*bcTop[I+nx-1] - (dy[ny-1]+dy[ny-2])*q[Iv]/dx[I] + dy[ny-1]*q[Iv-nx]/dx[I] )/( dy[ny-2]*dy[ny-1]*(dy[ny-2]+dy[ny-1]) ) \
				);
	/// Calculate rn
	rn[Iv] = ( q[Iv]/(dx[I]*dt) + cTerm + dTerm ) * 0.5*(dy[ny-2] + dy[ny-1]);
}

/**
* \brief To be documented
*/
__global__
void convectionTermUBottomTop(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcBottom, real *bcTop, real *bcLeft, real *bcRight)
{
	int	I = blockIdx.x*blockDim.x + threadIdx.x;
	
	if(I >= nx-1) return;
	/// Boundary Conditions for u at the Bottom ***********************************************************************************
	int		Iu,
			Iv	= I + (nx-1)*ny;	
	real	u = q[I]/dy[0],
			u0x, u1x,
			ul, ur,
			Hn, cTerm, dTerm;
	/// Calculate the Diffusion Term and Velocity Components for the convection term
	if(I==0){
		u0x	= 0.5*(bcLeft[0] + u);
		u1x	= 0.5*(u + q[I+1]/dy[0]);
		ul	= bcLeft[0];
		ur	= q[I+1]/dy[0];
	}
	else if(I==nx-2){
		u0x	= 0.5*(q[I-1]/dy[0] + u);
		u1x	= 0.5*(u + bcRight[0]);
		ul	= q[I-1]/dy[0];
		ur	= bcRight[0];
	}
	else{
		u0x	= 0.5*(q[I-1]/dy[0] + u);
		u1x	= 0.5*(u + q[I+1]/dy[0]);
		ul	= q[I-1]/dy[0];
		ur	= q[I+1]/dy[0];
	}
	/// Convection Term
	Hn = H[I];
	H[I] = - ( u1x*u1x - u0x*u0x )/( 0.5*(dx[I]+dx[I+1]) ) \
			- ( 0.5*(q[Iv]/dx[I] + q[Iv+1]/dx[I+1]) * 0.5*(u + q[I+(nx-1)]/dy[1]) - 0.5*(bcBottom[I+(nx-1)] + bcBottom[I+1+(nx-1)])*bcBottom[I] )/(dy[0]);
	cTerm = gamma*H[I] + zeta*Hn;
	//cTerm = H[I]; // 1st order Euler **************************** DID NOT CHANGE I TO Iu HERE
	// Diffusion Term
	dTerm = alpha*nu*2.0*( \
					( dx[I]*ur - (dx[I]+dx[I+1])*u + dx[I+1]*ul )/( dx[I] * (dx[I]+dx[I+1]) * dx[I+1] ) \
					+ 4.0*( dy[0]*q[I+(nx-1)]/dy[1] - (2.0*dy[0]+dy[1])*u + (dy[0]+dy[1])*bcBottom[I] ) \
						/( dy[0] * (2.0*dy[0]+dy[1]) * (dy[0]+dy[1]) )
				);
	/// Calculate rn
	rn[I] = ( u/dt + cTerm + dTerm ) * 0.5*(dx[I] + dx[I+1]);
	/// Boundary conditions for u at the Top **************************************************************************************
	Iu	= I + (ny-1)*(nx-1);
	Iv	= I + (ny-2)*nx + (nx-1)*ny;
	u = q[Iu]/dy[ny-1];
	/// Calculate the Diffusion Term and Velocity Components for the convection term
	if(I==0){
		u0x	= 0.5*(bcLeft[ny-1] + u);
		u1x	= 0.5*(u + q[Iu+1]/dy[ny-1]);
		ul	= bcLeft[ny-1];
		ur	= q[Iu+1]/dy[ny-1];
	}
	else if(I==nx-2){
		u0x	= 0.5*(q[Iu-1]/dy[ny-1] + u);
		u1x	= 0.5*(u + bcRight[ny-1]);
		ul	= q[Iu-1]/dy[ny-1];
		ur	= bcRight[ny-1];
	}
	else{
		u0x	= 0.5*(q[Iu-1]/dy[ny-1] + u);
		u1x	= 0.5*(u + q[Iu+1]/dy[ny-1]);
		ul	= q[Iu-1]/dy[ny-1];
		ur	= q[Iu+1]/dy[ny-1];
	}
	/// Convection Term
	Hn = H[Iu];
	H[Iu] = - ( u1x*u1x - u0x*u0x )/( 0.5*(dx[I]+dx[I+1]) ) \
			- ( bcTop[I]*0.5*(bcTop[I+(nx-1)]+bcTop[I+1+(nx-1)]) - 0.5*(q[Iv]/dx[I] + q[Iv+1]/dx[I+1])*0.5*(u + q[Iu-(nx-1)]/dy[ny-2]) )/(dy[ny-1]);
	cTerm = gamma*H[Iu] + zeta*Hn;
	//cTerm = H[Iu];	/// 1st order Euler  ************************* CHANGED I TO Iu HERE
	// Diffusion Term
	dTerm = alpha*nu*2.0*( \
					( dx[I]*ur - (dx[I]+dx[I+1])*u + dx[I+1]*ul )/( dx[I] * (dx[I]+dx[I+1]) * dx[I+1] ) \
					+ 4.0*( (dy[ny-2]+dy[ny-1])*bcTop[I] - (dy[ny-2]+2.0*dy[ny-1])*u + (dy[ny-1])*q[Iu-(nx-1)]/dy[ny-2] ) \
						/( (dy[ny-2]+dy[ny-1])*(dy[ny-1])*(dy[ny-2]+2.0*dy[ny-1]) )
				);
	/// Calculate rn
	rn[Iu] = ( u/dt + cTerm + dTerm) * 0.5*(dx[I] + dx[I+1]);
}

/**
* \brief To be documented
*/
__global__
void convectionTermULeftRight(real *rn, real *H, real *q, \
                              int nx, int ny, real *dx, real *dy, \
                              real dt, real gamma, real zeta, real alpha, real nu, \
                              real *bcLeft, real *bcRight)
{
	int	I = blockIdx.x*blockDim.x + threadIdx.x;
	
	if(I==0 || I >= ny-1) return;
	/// Boundary Conditions for u at the Left *************************************************************************************
	int		Iu = I*(nx-1),
			Iv = I*(nx) + (nx-1)*ny;
	real	Hn, cTerm, dTerm;
	/// Convection Term
	Hn = H[Iu];
	H[Iu] = -( \
				( 0.5*(q[Iu] + q[Iu+1])/dy[I] ) * ( 0.5*(q[Iu] + q[Iu+1])/dy[I] ) \
				- ( 0.5*(bcLeft[I] + q[Iu]/dy[I]) ) * ( 0.5*(bcLeft[I] + q[Iu]/dy[I]) ) \
			 )/(0.5*(dx[0]+dx[1])) \
	 		-( \
				( 0.5*(q[Iv]/dx[0] + q[Iv+1]/dx[1]) ) * ( 0.5*(q[Iu]/dy[I] + q[Iu+(nx-1)]/dy[I+1]) ) \
				- ( 0.5*(q[Iv-nx]/dx[0] + q[Iv+1-nx]/dx[1]) ) * ( 0.5*(q[Iu-(nx-1)]/dy[I-1] + q[Iu]/dy[I]) ) \
			 )/(dy[I]);
	cTerm = gamma*H[Iu] + zeta*Hn;
	//cTerm = H[Iu];	/// 1st order Euler
	/// Diffusion Term
	dTerm = alpha*nu*2.0*( \
					( dx[0]*q[Iu+1]/dy[I] - (dx[0]+dx[1])*q[Iu]/dy[I] + dx[1]*bcLeft[I] )/( dx[0]*dx[1]*(dx[0]+dx[1]) ) \
					+ 4.0 * ( (dy[I-1]+dy[I])*q[Iu+(nx-1)]/dy[I+1] - (dy[I-1]+2.0*dy[I]+dy[I+1])*q[Iu]/dy[I] + (dy[I]+dy[I+1])*q[Iu-(nx-1)]/dy[I-1] ) \
						/( (dy[I-1]+dy[I]) * (dy[I]+dy[I+1]) * (dy[I-1]+2.0*dy[I]+dy[I+1]) ) \
				);
	/// Calculate rn
	rn[Iu] = ( q[Iu]/(dy[I]*dt) + cTerm + dTerm ) * 0.5*(dx[0] + dx[1]);
	/// Boundary conditions for u at the Right ************************************************************************************
	Iu = I*(nx-1) + (nx-2);
	Iv = I*(nx) + (nx-1)*ny + (nx-2);
	/// Convection Term
	H[Iu] = -( \
				( 0.5*(q[Iu]/dy[I] + bcRight[I]) ) * ( 0.5*(q[Iu]/dy[I] + bcRight[I]) ) \
				- ( 0.5*(q[Iu-1] + q[Iu])/dy[I] ) * ( 0.5*(q[Iu-1] + q[Iu])/dy[I] ) \
			 )/(0.5*(dx[nx-2]+dx[nx-1])) \
	 		-( \
				( 0.5*(q[Iv]/dx[nx-2] + q[Iv+1]/dx[nx-1]) ) * ( 0.5*(q[Iu]/dy[I] + q[Iu+(nx-1)]/dy[I+1]) ) \
				- ( 0.5*(q[Iv-nx]/dx[nx-2] + q[Iv+1-nx]/dx[nx-1]) ) * ( 0.5*(q[Iu-(nx-1)]/dy[I-1] + q[Iu]/dy[I]) ) \
			 )/(dy[I]);
	cTerm = gamma*H[Iu] + zeta*Hn;
	//cTerm = H[Iu];	/// 1st order Euler
	/// Diffusion Term
	dTerm = alpha*nu*2.0*( \
					( dx[nx-2]*bcRight[I] - (dx[nx-2]+dx[nx-1])*q[Iu]/dy[I] + dx[nx-1]*q[Iu-1]/dy[I] )/( dx[nx-2]*dx[nx-1]*(dx[nx-2]+dx[nx-1]) ) \
					+ 4.0 * ( (dy[I-1]+dy[I])*q[Iu+(nx-1)]/dy[I+1] - (dy[I-1]+2.0*dy[I]+dy[I+1])*q[Iu]/dy[I] + (dy[I]+dy[I+1])*q[Iu-(nx-1)]/dy[I-1] ) \
						/( (dy[I-1]+dy[I]) * (dy[I]+dy[I+1]) * (dy[I-1]+2.0*dy[I]+dy[I+1]) ) \
				);
	/// Calculate rn
	rn[Iu] = ( q[Iu]/(dy[I]*dt) + cTerm + dTerm ) * 0.5*(dx[nx-2] + dx[nx-1]);
}

/**
* \brief To be documented
*/
__global__
void convectionTermVLeftRight(real *rn, real *H, real *q, \
                             int nx, int ny, real *dx, real *dy, \
                             real dt, real gamma, real zeta, real alpha, real nu,\
                             real *bcBottom, real *bcTop, real *bcLeft, real *bcRight)
{
	int	I = blockIdx.x*blockDim.x + threadIdx.x;
	
	if(I > ny-2) return;
	/// Boundary Conditions for v at the Left *************************************************************************************
	int		Iu = I*(nx-1),
			Iv = I*nx + (nx-1)*ny;
	real	vb, vt, v0y, v1y, v,
			Hn, cTerm, dTerm;
	v = q[Iv]/dx[0];
	if(I==0){
		v0y = 0.5*(bcBottom[nx-1] + v);
		v1y = 0.5*(v + q[Iv+nx]/dx[0]);
		vb = bcBottom[nx-1];
		vt = q[Iv+nx]/dx[0];
	}
	else if(I==ny-2){
		v0y = 0.5*(q[Iv-nx]/dx[0] + v);
		v1y = 0.5*(q[Iv]/dx[0] + bcTop[nx-1]);
		vb = q[Iv-nx]/dx[0];
		vt = bcTop[nx-1];
	}
	else{
		v0y = 0.5*(q[Iv-nx]/dx[0] + v);
		v1y = 0.5*(v + q[Iv+nx]/dx[0]);
		vb = q[Iv-nx]/dx[0];
		vt = q[Iv+nx]/dx[0];
	}
	__syncthreads();
	/// Convection Term
	Hn = H[Iv];
	H[Iv] = -( 0.5*(v + q[Iv+1]/dx[1])*0.5*(q[Iu]/dy[I] + q[Iu+(nx-1)]/dy[I+1]) - bcLeft[I+ny]*0.5*(bcLeft[I] + bcLeft[I+1]) )/dx[0] \
	 		-( v1y*v1y - v0y*v0y )/(0.5*(dy[I] + dy[I+1]));
	cTerm = gamma*H[Iv] + zeta*Hn;
	//cTerm = H[Iv];	/// 1st order Euler
	/// Diffusion Term
	dTerm = alpha*nu*2.0*( \
					4.0 * ( (dx[0])*q[Iv+1]/dx[1] - (2.0*dx[0]+dx[1])*v + (dx[0]+dx[1])*bcLeft[I+ny] ) \
						/( (dx[0]) * (dx[0]+dx[1]) * (2.0*dx[0]+dx[1]) ) \
					+ ( dy[I]*vt - (dy[I]+dy[I+1])*v + dy[I+1]*vb )/( dy[I]*dy[I+1]*(dy[I]+dy[I+1]) ) \
				);
	/// Calculate rn
	rn[Iv] = ( v/dt + cTerm + dTerm ) * 0.5*(dy[I] + dy[I+1]);
	/// Boundary Conditions for v at the Right ************************************************************************************
	Iu = I*(nx-1) + (nx-1);
	Iv = I*nx + (nx-1)*ny + (nx-1);
	v = q[Iv]/dx[nx-1];
	if(I==0){
		v0y = 0.5*(bcBottom[nx-1+(nx-1)] + v);
		v1y = 0.5*(v + q[Iv+nx]/dx[nx-1]);
		vb = bcBottom[nx-1+(nx-1)];
		vt = q[Iv+nx]/dx[nx-1];
	}
	else if(I==ny-2){
		v0y = 0.5*(q[Iv-nx]/dx[nx-1] + v);
		v1y = 0.5*(q[Iv]/dx[nx-1] + bcTop[nx-1+(nx-1)]);
		vb = q[Iv-nx]/dx[nx-1];
		vt = bcTop[nx-1+(nx-1)];
	}
	else{
		v0y = 0.5*(q[Iv-nx]/dx[nx-1] + v);
		v1y = 0.5*(v + q[Iv+nx]/dx[nx-1]);
		vb = q[Iv-nx]/dx[nx-1];
		vt = q[Iv+nx]/dx[nx-1];
	}
	__syncthreads();
	/// Convection Term
	Hn = H[Iv];
	H[Iv] = -( bcRight[I+ny]*0.5*(bcRight[I]+bcRight[I+1]) - 0.5*(q[Iv-1]/dx[nx-2] + v)*0.5*(q[Iu-1]/dy[I] + q[Iu-1+(nx-1)]/dy[I+1]) )/dx[nx-1] \
	 		-( v1y*v1y - v0y*v0y )/(0.5*(dy[I] + dy[I+1]));
	cTerm = gamma*H[Iv] + zeta*Hn;
	//cTerm = H[Iv];	/// 1st order Euler
	/// Diffusion Term
	dTerm = alpha*nu*2.0*( \
					4.0 * ( (dx[nx-1]+dx[nx-2])*bcRight[I+ny] - (2.0*dx[nx-1]+dx[nx-2])*v + (dx[nx-1])*q[Iv-1]/dx[nx-2] ) \
						/( (dx[nx-1]) * (dx[nx-1]+dx[nx-2]) * (2.0*dx[nx-1]+dx[nx-2]) ) \
					+ ( dy[I]*vt - (dy[I]+dy[I+1])*v + dy[I+1]*vb )/( dy[I]*dy[I+1]*(dy[I]+dy[I+1]) ) \
				);
	/// Calculate rn
	rn[Iv] = ( v/dt + cTerm + dTerm ) * 0.5*(dy[I] + dy[I+1]);
}

} // end of namespace kernels
