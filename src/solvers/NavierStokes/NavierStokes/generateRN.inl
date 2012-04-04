#define BSZ 16
// CUDA kernel to generate Hx and rn
__global__ void convection_term_u_cuda(real *rn, real *H, real *q, int nx, int ny, real *dx, real *dy, real dt, real d_coeff)
{
	int bx	= blockIdx.x,
		by	= blockIdx.y,
		i	= threadIdx.x,
		j	= threadIdx.y;
	
	// work out global index of first point in block
	int I = (BSZ-2)*bx + i,
		J = (BSZ-2)*by + j;
	
	if (I >= nx-1 || J >= ny) {
		return;
	}

	int		N_u		= (nx-1)*ny,
			Gidx_x	= J*(nx-1) + I,
			Gidx_y	= (J-1)*nx + I + N_u;
			
	real	c_term, d_term, Hxn;
	
	__shared__ real		u[BSZ][BSZ],
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
		//c_term = 1.5*Hx[j][i] - 0.5*Hxn;
		c_term = Hx[j][i]; // 1st order Euler
		d_term = d_coeff*( \
						 ( Dx[j][i]*u[j][i+1] - (Dx[j][i]+Dx[j][i+1])*u[j][i] + Dx[j][i+1]*u[j][i-1] )/( Dx[j][i]*Dx[j][i+1]*(Dx[j][i]+Dx[j][i+1]) ) \
					   + 4.0*( (Dy[j][i]+Dy[j-1][i])*u[j+1][i] - (Dy[j-1][i] + 2.0*Dy[j][i] + Dy[j+1][i])*u[j][i] + (Dy[j][i]+Dy[j+1][i])*u[j-1][i] ) \
							/( (Dy[j][i]+Dy[j-1][i]) * (Dy[j-1][i] + 2.0*Dy[j][i] + Dy[j+1][i]) * (Dy[j][i]+Dy[j+1][i]) ) \
					     );
		rn[Gidx_x] = (u[j][i]/dt + c_term + d_term) * 0.5*(Dx[j][i]+Dx[j][i+1]);
	}
}

// CUDA kernel to generate Hy and rn
__global__ void convection_term_v_cuda(real *rn, real *H, real *q, int nx, int ny, real *dx, real *dy, real dt, real d_coeff)
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
			
	real	c_term, d_term, Hyn;
	
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
		//c_term = 1.5*Hy[j][i] - 0.5*Hyn;
		c_term = Hy[j][i]; // 1st order Euler
		d_term = d_coeff*( \
						4.0*( (Dx[j][i-1]+Dx[j][i])*v[j][i+1] - (Dx[j][i-1]+2.0*Dx[j][i]+Dx[j][i+1])*v[j][i] + (Dx[j][i]+Dx[j][i+1])*v[j][i-1] ) \
							/( (Dx[j][i-1]+Dx[j][i]) * (Dx[j][i]+Dx[j][i+1]) * (Dx[j][i-1]+2.0*Dx[j][i]+Dx[j][i+1]) ) \
						+ ( Dy[j][i]*v[j+1][i] - (Dy[j][i]+Dy[j+1][i])*v[j][i] + Dy[j+1][i]*v[j-1][i] )/( Dy[j][i]*Dy[j+1][i]*(Dy[j][i]+Dy[j+1][i]) ) \
					);
		rn[Gidx_y] = (v[j][i]/dt + c_term + d_term) * 0.5*(Dy[j][i]+Dy[j+1][i]);
	}
}

__global__ void convection_term_v_bottomtop_cuda(real *rn, real *H, real *q, \
												int nx, int ny, real *dx, real *dy, \
												real dt, real d_coeff, \
												real *bc_bottom, real *bc_top)
{
	int	I = blockIdx.x*blockDim.x + threadIdx.x;
	
	if(I==0 || I >= nx-1) return;
	/// Boundary Conditions for v at the Bottom ***********************************************************************************
	int		Iu, Iv = I + (nx-1)*ny;	
	real	Hn, c_term, d_term;
	/// Convection Term
	Hn = H[Iv];
	H[Iv] = -( \
				( 0.5*(q[Iv]/dx[I]+q[Iv+1]/dx[I+1]) ) * ( 0.5*(q[I]/dy[0]+q[I+(nx-1)]/dy[1]) ) \
				- ( 0.5*(q[Iv-1]/dx[I-1]+q[Iv]/dx[I]) ) * ( 0.5*(q[I-1]/dy[0]+q[I-1+(nx-1)]/dy[1]) ) \
			 )/dx[I] \
	 		-( \
				(0.5*(q[Iv] + q[Iv+nx])/dx[I]) * (0.5*(q[Iv] + q[Iv+nx])/dx[I]) \
				- (0.5*(bc_bottom[I+nx-1] + q[Iv]/dx[I])) * (0.5*(bc_bottom[I+nx-1] + q[Iv]/dx[I])) \
			 )/(0.5*(dy[0] + dy[1]));
	//c_term = 1.5*H[Iv] - 0.5*Hn;	/// 2nd order Adams-Bashforth
	c_term = H[Iv];	/// 1st order Euler
	/// Diffusion Term
	d_term = d_coeff*( \
					4.0 * ( (dx[I-1]+dx[I])*q[Iv+1]/dx[I+1] - (dx[I-1]+2.0*dx[I]+dx[I+1])*q[Iv]/dx[I] + (dx[I]+dx[I+1])*q[Iv-1]/dx[I-1] ) \
						/( (dx[I-1]+dx[I]) * (dx[I]+dx[I+1]) * (dx[I-1]+2.0*dx[I]+dx[I+1]) ) \
					+ ( dy[0]*q[Iv+nx]/dx[I] - (dy[0]+dy[1])*q[Iv]/dx[I] + dy[1]*bc_bottom[I+nx-1] )/( dy[0]*dy[1]*(dy[0]+dy[1]) ) \
				);
	/// Calculate rn
	rn[Iv] = ( q[Iv]/(dx[I]*dt) + c_term + d_term ) * 0.5*(dy[0] + dy[1]);
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
				( 0.5*(q[Iv]/dx[I] + bc_top[I+nx-1]) )*( 0.5*(q[Iv]/dx[I] + bc_top[I+nx-1]) ) \
				- ( 0.5*(q[Iv-nx] + q[Iv])/dx[I] )*( 0.5*(q[Iv-nx] + q[Iv])/dx[I] )
			 )/(0.5*(dy[ny-2] + dy[ny-1]));
	//c_term = 1.5*H[Iv] - 0.5*Hn;	/// 2nd order Adams-Bashforth
	c_term = H[Iv];	/// 1st order Euler
	/// Diffusion Term
	d_term = d_coeff*( \
					4.0*( (dx[I-1]+dx[I])*q[Iv+1]/dx[I+1] - (dx[I-1]+2.0*dx[I]+dx[I+1])*q[Iv]/dx[I] + (dx[I]+dx[I+1])*q[Iv-1]/dx[I-1] ) \
						/( (dx[I-1]+dx[I]) * (dx[I]+dx[I+1]) * (dx[I-1]+2.0*dx[I]+dx[I+1]) ) \
					+ ( dy[ny-2]*bc_top[I+nx-1] - (dy[ny-1]+dy[ny-2])*q[Iv]/dx[I] + dy[ny-1]*q[Iv-nx]/dx[I] )/( dy[ny-2]*dy[ny-1]*(dy[ny-2]+dy[ny-1]) ) \
				);
	/// Calculate rn
	rn[Iv] = ( q[Iv]/(dx[I]*dt) + c_term + d_term ) * 0.5*(dy[ny-2] + dy[ny-1]);
}

__global__ void convection_term_u_bottomtop_cuda(real *rn, real *H, real *q, \
												int nx, int ny, real *dx, real *dy, \
												real dt, real d_coeff, \
												real *bc_bottom, real *bc_top, real *bc_left, real *bc_right)
{
	int	I = blockIdx.x*blockDim.x + threadIdx.x;
	
	if(I >= nx-1) return;
	/// Boundary Conditions for u at the Bottom ***********************************************************************************
	int		Iu,
			Iv	= I + (nx-1)*ny;	
	real	u = q[I]/dy[0],
			u0x, u1x,
			ul, ur,
			Hn, c_term, d_term;
	/// Calculate the Diffusion Term and Velocity Components for the convection term
	if(I==0){
		u0x	= 0.5*(bc_left[0] + u);
		u1x	= 0.5*(u + q[I+1]/dy[0]);
		ul	= bc_left[0];
		ur	= q[I+1]/dy[0];
	}
	else if(I==nx-2){
		u0x	= 0.5*(q[I-1]/dy[0] + u);
		u1x	= 0.5*(u + bc_right[0]);
		ul	= q[I-1]/dy[0];
		ur	= bc_right[0];
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
			- ( 0.5*(q[Iv]/dx[I] + q[Iv+1]/dx[I+1]) * 0.5*(u + q[I+(nx-1)]/dy[1]) - 0.5*(bc_bottom[I+(nx-1)] + bc_bottom[I+1+(nx-1)])*bc_bottom[I] )/(dy[0]);
	//c_term = 1.5*H[I] - 0.5*Hn;	/// 2nd order Adams-Bashforth
	c_term = H[I]; // 1st order Euler **************************** DID NOT CHANGE I TO Iu HERE
	// Diffusion Term
	d_term = d_coeff*( \
					( dx[I]*ur - (dx[I]+dx[I+1])*u + dx[I+1]*ul )/( dx[I] * (dx[I]+dx[I+1]) * dx[I+1] ) \
					+ 4.0*( dy[0]*q[I+(nx-1)]/dy[1] - (2.0*dy[0]+dy[1])*u + (dy[0]+dy[1])*bc_bottom[I] ) \
						/( dy[0] * (2.0*dy[0]+dy[1]) * (dy[0]+dy[1]) )
				);
	/// Calculate rn
	rn[I] = ( u/dt + c_term + d_term ) * 0.5*(dx[I] + dx[I+1]);
	/// Boundary conditions for u at the Top **************************************************************************************
	Iu	= I + (ny-1)*(nx-1);
	Iv	= I + (ny-2)*nx + (nx-1)*ny;
	u = q[Iu]/dy[ny-1];
	/// Calculate the Diffusion Term and Velocity Components for the convection term
	if(I==0){
		u0x	= 0.5*(bc_left[ny-1] + u);
		u1x	= 0.5*(u + q[Iu+1]/dy[ny-1]);
		ul	= bc_left[ny-1];
		ur	= q[Iu+1]/dy[ny-1];
	}
	else if(I==nx-2){
		u0x	= 0.5*(q[Iu-1]/dy[ny-1] + u);
		u1x	= 0.5*(u + bc_right[ny-1]);
		ul	= q[Iu-1]/dy[ny-1];
		ur	= bc_right[ny-1];
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
			- ( bc_top[I]*0.5*(bc_top[I+(nx-1)]+bc_top[I+1+(nx-1)]) - 0.5*(q[Iv]/dx[I] + q[Iv+1]/dx[I+1])*0.5*(u + q[Iu-(nx-1)]/dy[ny-2]) )/(dy[ny-1]);
	//c_term = 1.5*H[Iu] - 0.5*Hn;	/// 2nd order Adams-Bashforth
	c_term = H[Iu];	/// 1st order Euler  ************************* CHANGED I TO Iu HERE
	// Diffusion Term
	d_term = d_coeff*( \
					( dx[I]*ur - (dx[I]+dx[I+1])*u + dx[I+1]*ul )/( dx[I] * (dx[I]+dx[I+1]) * dx[I+1] ) \
					+ 4.0*( (dy[ny-2]+dy[ny-1])*bc_top[I] - (dy[ny-2]+2.0*dy[ny-1])*u + (dy[ny-1])*q[Iu-(nx-1)]/dy[ny-2] ) \
						/( (dy[ny-2]+dy[ny-1])*(dy[ny-1])*(dy[ny-2]+2.0*dy[ny-1]) )
				);
	/// Calculate rn
	rn[Iu] = ( u/dt + c_term + d_term) * 0.5*(dx[I] + dx[I+1]);
}

__global__ void convection_term_u_leftright_cuda(real *rn, real *H, real *q, \
												int nx, int ny, real *dx, real *dy, \
												real dt, real d_coeff, \
												real *bc_left, real *bc_right)
{
	int	I = blockIdx.x*blockDim.x + threadIdx.x;
	
	if(I==0 || I >= ny-1) return;
	/// Boundary Conditions for u at the Left *************************************************************************************
	int		Iu = I*(nx-1),
			Iv = I*(nx) + (nx-1)*ny;
	real	Hn, c_term, d_term;
	/// Convection Term
	Hn = H[Iu];
	H[Iu] = -( \
				( 0.5*(q[Iu] + q[Iu+1])/dy[I] ) * ( 0.5*(q[Iu] + q[Iu+1])/dy[I] ) \
				- ( 0.5*(bc_left[I] + q[Iu]/dy[I]) ) * ( 0.5*(bc_left[I] + q[Iu]/dy[I]) ) \
			 )/(0.5*(dx[0]+dx[1])) \
	 		-( \
				( 0.5*(q[Iv]/dx[0] + q[Iv+1]/dx[1]) ) * ( 0.5*(q[Iu]/dy[I] + q[Iu+(nx-1)]/dy[I+1]) ) \
				- ( 0.5*(q[Iv-nx]/dx[0] + q[Iv+1-nx]/dx[1]) ) * ( 0.5*(q[Iu-(nx-1)]/dy[I-1] + q[Iu]/dy[I]) ) \
			 )/(dy[I]);
	//c_term = 1.5*H[Iu] - 0.5*Hn;	/// 2nd order Adams-Bashforth
	c_term = H[Iu];	/// 1st order Euler
	/// Diffusion Term
	d_term = d_coeff*( \
					( dx[0]*q[Iu+1]/dy[I] - (dx[0]+dx[1])*q[Iu]/dy[I] + dx[1]*bc_left[I] )/( dx[0]*dx[1]*(dx[0]+dx[1]) ) \
					+ 4.0 * ( (dy[I-1]+dy[I])*q[Iu+(nx-1)]/dy[I+1] - (dy[I-1]+2.0*dy[I]+dy[I+1])*q[Iu]/dy[I] + (dy[I]+dy[I+1])*q[Iu-(nx-1)]/dy[I-1] ) \
						/( (dy[I-1]+dy[I]) * (dy[I]+dy[I+1]) * (dy[I-1]+2.0*dy[I]+dy[I+1]) ) \
				);
	/// Calculate rn
	rn[Iu] = ( q[Iu]/(dy[I]*dt) + c_term + d_term ) * 0.5*(dx[0] + dx[1]);
	/// Boundary conditions for u at the Right ************************************************************************************
	Iu = I*(nx-1) + (nx-2);
	Iv = I*(nx) + (nx-1)*ny + (nx-2);
	/// Convection Term
	H[Iu] = -( \
				( 0.5*(q[Iu]/dy[I] + bc_right[I]) ) * ( 0.5*(q[Iu]/dy[I] + bc_right[I]) ) \
				- ( 0.5*(q[Iu-1] + q[Iu])/dy[I] ) * ( 0.5*(q[Iu-1] + q[Iu])/dy[I] ) \
			 )/(0.5*(dx[nx-2]+dx[nx-1])) \
	 		-( \
				( 0.5*(q[Iv]/dx[nx-2] + q[Iv+1]/dx[nx-1]) ) * ( 0.5*(q[Iu]/dy[I] + q[Iu+(nx-1)]/dy[I+1]) ) \
				- ( 0.5*(q[Iv-nx]/dx[nx-2] + q[Iv+1-nx]/dx[nx-1]) ) * ( 0.5*(q[Iu-(nx-1)]/dy[I-1] + q[Iu]/dy[I]) ) \
			 )/(dy[I]);
	//c_term = 1.5*H[Iu] - 0.5*Hn;	/// 2nd order Adams-Bashforth
	c_term = H[Iu];	/// 1st order Euler
	/// Diffusion Term
	d_term = d_coeff*( \
					( dx[nx-2]*bc_right[I] - (dx[nx-2]+dx[nx-1])*q[Iu]/dy[I] + dx[nx-1]*q[Iu-1]/dy[I] )/( dx[nx-2]*dx[nx-1]*(dx[nx-2]+dx[nx-1]) ) \
					+ 4.0 * ( (dy[I-1]+dy[I])*q[Iu+(nx-1)]/dy[I+1] - (dy[I-1]+2.0*dy[I]+dy[I+1])*q[Iu]/dy[I] + (dy[I]+dy[I+1])*q[Iu-(nx-1)]/dy[I-1] ) \
						/( (dy[I-1]+dy[I]) * (dy[I]+dy[I+1]) * (dy[I-1]+2.0*dy[I]+dy[I+1]) ) \
				);
	/// Calculate rn
	rn[Iu] = ( q[Iu]/(dy[I]*dt) + c_term + d_term ) * 0.5*(dx[nx-2] + dx[nx-1]);
}

__global__ void convection_term_v_leftright_cuda(real *rn, real *H, real *q, \
												int nx, int ny, real *dx, real *dy, \
												real dt, real d_coeff, \
												real *bc_bottom, real *bc_top, real *bc_left, real *bc_right)
{
	int	I = blockIdx.x*blockDim.x + threadIdx.x;
	
	if(I > ny-2) return;
	/// Boundary Conditions for v at the Left *************************************************************************************
	int		Iu = I*(nx-1),
			Iv = I*nx + (nx-1)*ny;
	real	vb, vt, v0y, v1y, v,
			Hn, c_term, d_term;
	v = q[Iv]/dx[0];
	if(I==0){
		v0y = 0.5*(bc_bottom[nx-1] + v);
		v1y = 0.5*(v + q[Iv+nx]/dx[0]);
		vb = bc_bottom[nx-1];
		vt = q[Iv+nx]/dx[0];
	}
	else if(I==ny-2){
		v0y = 0.5*(q[Iv-nx]/dx[0] + v);
		v1y = 0.5*(q[Iv]/dx[0] + bc_top[nx-1]);
		vb = q[Iv-nx]/dx[0];
		vt = bc_top[nx-1];
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
	H[Iv] = -( 0.5*(v + q[Iv+1]/dx[1])*0.5*(q[Iu]/dy[I] + q[Iu+(nx-1)]/dy[I+1]) - bc_left[I+ny]*0.5*(bc_left[I] + bc_left[I+1]) )/dx[0] \
	 		-( v1y*v1y - v0y*v0y )/(0.5*(dy[I] + dy[I+1]));
	//c_term = 1.5*H[Iv] - 0.5*Hn;	/// 2nd order Adams-Bashforth
	c_term = H[Iv];	/// 1st order Euler
	/// Diffusion Term
	d_term = d_coeff*( \
					4.0 * ( (dx[0])*q[Iv+1]/dx[1] - (2.0*dx[0]+dx[1])*v + (dx[0]+dx[1])*bc_left[I+ny] ) \
						/( (dx[0]) * (dx[0]+dx[1]) * (2.0*dx[0]+dx[1]) ) \
					+ ( dy[I]*vt - (dy[I]+dy[I+1])*v + dy[I+1]*vb )/( dy[I]*dy[I+1]*(dy[I]+dy[I+1]) ) \
				);
	/// Calculate rn
	rn[Iv] = ( v/dt + c_term + d_term ) * 0.5*(dy[I] + dy[I+1]);
	/// Boundary Conditions for v at the Right ************************************************************************************
	Iu = I*(nx-1) + (nx-1);
	Iv = I*nx + (nx-1)*ny + (nx-1);
	v = q[Iv]/dx[nx-1];
	if(I==0){
		v0y = 0.5*(bc_bottom[nx-1+(nx-1)] + v);
		v1y = 0.5*(v + q[Iv+nx]/dx[nx-1]);
		vb = bc_bottom[nx-1+(nx-1)];
		vt = q[Iv+nx]/dx[nx-1];
	}
	else if(I==ny-2){
		v0y = 0.5*(q[Iv-nx]/dx[nx-1] + v);
		v1y = 0.5*(q[Iv]/dx[nx-1] + bc_top[nx-1+(nx-1)]);
		vb = q[Iv-nx]/dx[nx-1];
		vt = bc_top[nx-1+(nx-1)];
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
	H[Iv] = -( bc_right[I+ny]*0.5*(bc_right[I]+bc_right[I+1]) - 0.5*(q[Iv-1]/dx[nx-2] + v)*0.5*(q[Iu-1]/dy[I] + q[Iu-1+(nx-1)]/dy[I+1]) )/dx[nx-1] \
	 		-( v1y*v1y - v0y*v0y )/(0.5*(dy[I] + dy[I+1]));
	//c_term = 1.5*H[Iv] - 0.5*Hn;	/// 2nd order Adams-Bashforth
	c_term = H[Iv];	/// 1st order Euler
	/// Diffusion Term
	d_term = d_coeff*( \
					4.0 * ( (dx[nx-1]+dx[nx-2])*bc_right[I+ny] - (2.0*dx[nx-1]+dx[nx-2])*v + (dx[nx-1])*q[Iv-1]/dx[nx-2] ) \
						/( (dx[nx-1]) * (dx[nx-1]+dx[nx-2]) * (2.0*dx[nx-1]+dx[nx-2]) ) \
					+ ( dy[I]*vt - (dy[I]+dy[I+1])*v + dy[I+1]*vb )/( dy[I]*dy[I+1]*(dy[I]+dy[I+1]) ) \
				);
	/// Calculate rn
	rn[Iv] = ( v/dt + c_term + d_term ) * 0.5*(dy[I] + dy[I+1]);
}

template <>
void NavierStokesSolver<device_memory>::generateRN()
{
	// raw pointers for cup arrays
	real *H_r  = thrust::raw_pointer_cast(&H[0]),
	     *q_r  = thrust::raw_pointer_cast(&q[0]),
	     *rn_r = thrust::raw_pointer_cast(&rn[0]),
	     *dxD = thrust::raw_pointer_cast(&(domInfo->dxD[0])),
	     *dyD = thrust::raw_pointer_cast(&(domInfo->dyD[0]));
	
	real *xminus = thrust::raw_pointer_cast(&(bc[XMINUS][0])),
	     *xplus  = thrust::raw_pointer_cast(&(bc[XPLUS][0])),
	     *yminus = thrust::raw_pointer_cast(&(bc[YMINUS][0])),
	     *yplus  = thrust::raw_pointer_cast(&(bc[YPLUS][0]));
	
	real d_coeff = 0.0;//(1.0 - D.omega) * 2.0 * F.nu;
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	dim3 dimGridx( int( (nx-1-0.5)/(BSZ-2) ) + 1, int( (ny-0.5)/(BSZ-2) ) + 1 ),
	     dimGridy( int( (nx-0.5)/(BSZ-2) ) + 1, int( (ny-1-0.5)/(BSZ-2) ) + 1 );
	dim3 dimBlock(BSZ, BSZ);
	
	// call the kernel
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	convection_term_u_cuda <<<dimGridx, dimBlock>>> (rn_r, H_r, q_r, nx, ny, dxD, dyD, dt, d_coeff);
	convection_term_v_cuda <<<dimGridy, dimBlock>>> (rn_r, H_r, q_r, nx, ny, dxD, dyD, dt, d_coeff);
	
	dim3 dimGridbc(int((nx+ny-0.5)/(BSZ*BSZ))+1, 1);
	dim3 dimBlockbc(BSZ*BSZ, 1);
	
	convection_term_u_bottomtop_cuda<<<dimGridbc, dimBlockbc>>>(rn_r, H_r, q_r, nx, ny, dxD, dyD, dt, d_coeff, \
																yminus, yplus, xminus, xplus);
	convection_term_v_bottomtop_cuda<<<dimGridbc, dimBlockbc>>>(rn_r, H_r, q_r, nx, ny, dxD, dyD, dt, d_coeff, \
																yminus, yplus);

	convection_term_u_leftright_cuda<<<dimGridbc, dimBlockbc>>>(rn_r, H_r, q_r, nx, ny, dxD, dyD, dt, d_coeff, \
																xminus, xplus);
	convection_term_v_leftright_cuda<<<dimGridbc, dimBlockbc>>>(rn_r, H_r, q_r, nx, ny, dxD, dyD, dt, d_coeff, \
																yminus, yplus, xminus, xplus);
}

template <>
void NavierStokesSolver<host_memory>::generateRN()
{
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	int  numU = (nx-1)*ny;
	int  Iu = 0, Iv = 0;
	real east = 0, west = 0, north = 0, south = 0, Hn = 0, c_term = 0, d_term = 0, u = 0, v = 0;
	
	real *dx = thrust::raw_pointer_cast(&(domInfo->dx[0])),
	     *dy = thrust::raw_pointer_cast(&(domInfo->dy[0]));
	
	real *xminus = thrust::raw_pointer_cast(&(bc[XMINUS][0])),
	     *xplus  = thrust::raw_pointer_cast(&(bc[XPLUS][0])),
	     *yminus = thrust::raw_pointer_cast(&(bc[YMINUS][0])),
	     *yplus  = thrust::raw_pointer_cast(&(bc[YPLUS][0]));
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx-1; i++)
		{
			Iu = j*(nx-1)+i;
			Iv = j*nx + i + numU;
			Hn = H[Iu];
			u  = q[Iu]/dy[j];
			if(i==0)
				west = (xminus[j] + u)/2.0 * (xminus[j] + u)/2.0;
			else
				west = (q[Iu-1]/dy[j] + u)/2.0 * (q[Iu-1]/dy[j] + u)/2.0;
			if(i==nx-2)
				east = (u + xplus[j])/2.0 * (u + xplus[j])/2.0;
			else
				east = (u + q[Iu+1]/dy[j])/2.0 * (u + q[Iu+1]/dy[j])/2.0;
			if(j==0)
				south = yminus[i] * (yminus[i+(nx-1)]+yminus[i+1+(nx-1)])/2.0;
			else
				south = (q[Iu-(nx-1)]/dy[j-1] + u)/2.0 * (q[Iv-nx]/dx[i] + q[Iv-nx+1]/dx[i+1])/2.0;
			if(j==ny-1)
				north = yplus[i] * (yplus[i+(nx-1)]+yplus[i+1+(nx-1)])/2.0;
			else
				north = (u + q[Iu+(nx-1)]/dy[j+1])/2.0 * (q[Iv]/dx[i] + q[Iv+1]/dx[i+1])/2.0;
				
			H[Iu]  = - (east-west)/0.5*(dx[i]+dx[i+1]) - (north-south)/dy[j];
			c_term = 1.5*H[Iu] - 0.5*Hn;
			d_term = 0.0;
			rn[Iu] = (u/dt + c_term + d_term) * 0.5*(dx[i]+dx[i+1]);
		}
	}
	
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			Iv = j*nx+i + numU;
			Iv = j*(nx-1) + i;
			Hn = H[Iv];
			v  = q[Iv]/dx[i];
			if(j==0)
				south = (yminus[i+(nx-1)] + v)/2.0 * (yminus[i+(nx-1)] + v)/2.0;
			else
				south = (q[Iv-nx]/dx[i] + v)/2.0 * (q[Iv-nx]/dx[i] + v)/2.0;
			if(j==ny-2)
				north = (v + yplus[i+(nx-1)])/2.0 * (v + yplus[i+(nx-1)])/2.0;
			else
				north = (v + q[Iv+nx]/dx[i])/2.0 * (v + q[Iv+nx]/dx[i])/2.0;
			if(i==0)
				west = xminus[j+ny]*(xminus[j]+xminus[j+1])/2.0;
			else
				west = (q[Iv-1]/dx[i-1] + v)/2.0 * (q[Iu-1]/dy[j] + q[Iu-1+(nx-1)]/dy[j+1])/2.0;
			if(i==nx-1)
				east = xplus[j+ny]*(xplus[j]+xplus[j+1])/2.0;
			else
				east = (v + q[Iv+1]/dx[i+1])/2.0 * (q[Iu]/dy[j] + q[Iu+(nx-1)]/dy[j+1])/2.0;
				
			H[Iu]  = - (east-west)/dx[i] - (north-south)/0.5*(dy[j]+dy[j+1]);
			c_term = 1.5*H[Iv] - 0.5*Hn;
			d_term = 0.0;
			rn[Iu] = (v/dt + c_term + d_term) * 0.5*(dy[j]+dy[j+1]);
		}
	}
}

