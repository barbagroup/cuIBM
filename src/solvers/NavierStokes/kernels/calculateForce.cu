#include <solvers/NavierStokes/kernels/calculateForce.h>

#define BSZ 16

namespace kernels
{

__global__
void dragLeftRight(real *FxX, real *q, real *lambda, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy)
{
	int  idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx >= ncy)
		return;
	int  Ip = (J+idx)*nx + I,
	     Iu = (J+idx)*(nx-1) + (I-1);
	FxX[idx] = -(
	              // multiply the pressure with the surface area to get p dy
	              (lambda[Ip+ncx]-lambda[Ip-1])*dy[J+idx]
	              +
	              // divide q^2 by dy, so that just u^2 dy is obtained
	              (
	                  0.25*(q[Iu+ncx+1] + q[Iu+ncx])*(q[Iu+ncx+1] + q[Iu+ncx])
	                - 0.25*(q[Iu] + q[Iu-1])*(q[Iu] + q[Iu-1])
	              )/dy[J+idx]
	              -
	              // no multiplication or division since du/dx dy = dq/dx
	              nu*
	              (
	                  (q[Iu+ncx+1] - q[Iu+ncx])/dx[I+ncx]
	                - (q[Iu] - q[Iu-1])/dx[I-1]
	              )
	            );
}

__global__
void dragBottomTop(real *FxY, real *q, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy)
{
	int  idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx > ncx)
		return;
	int  Iu = J*(nx-1) + (I-1+idx),
	     Iv = (nx-1)*ny + (J-1)*nx + I+idx;
	FxY[idx] = -(
	              // multiply by dS
	              (
	                0.25 * ( q[Iu+ncy*(nx-1)]/dy[J+ncy] + q[Iu+(ncy-1)*(nx-1)]/dy[J+ncy-1] )
	                     * ( q[Iv+ncy*nx]/dx[I+idx] + q[Iv+ncy*nx-1]/dx[I+idx-1] )
	                -
	                0.25 * ( q[Iu]/dy[J] + q[Iu-(nx-1)]/dy[J-1] ) 
	                     * ( q[Iv]/dx[I+idx] + q[Iv-1]/dx[I+idx-1] )
	              )
	              -
	              // multiply by dS (cannot use the leftRight trick in this case)
	              nu*
	              (
	                (
	                  (q[Iu+ncy*(nx-1)]/dy[J+ncy] - q[Iu+(ncy-1)*(nx-1)]/dy[J+ncy-1])/2.0/(dy[J+ncy]+dy[J+ncy-1]) +
	                  (q[Iv+ncy*nx]/dx[I+idx] - q[Iv+ncy*nx-1]/dx[I+idx-1])/2.0/(dx[I+idx]+dx[I+idx-1])
	                ) 
	                -
	                (
	                  (q[Iu]/dy[J] - q[Iu-(nx-1)]/dy[J-1])/2.0/(dy[J]+dy[J-1]) +
	                  (q[Iv]/dx[I+idx] - q[Iv-1]/dx[I+idx-1])/2.0/(dx[I+idx]+dx[I+idx-1])
	                )
	              )
	            )*0.5*(dx[I+idx]+dx[I+idx-1]);

}

__global__
void dragUnsteady(real *FxU, real *q, real *qOld, real *dx, real *dy, real dt,
                   int nx, int ny, int I, int J, int ncx, int ncy)
{
	int  idx = threadIdx.x + blockIdx.x*blockDim.x;
	     
	if(idx >= (ncx+1)*ncy)
		return;
	
	int i = idx%(ncx+1),
	    j = idx/(ncx+1);
	    
	int Iu = (J+j)*(nx-1) + (I-1+i);
	
	FxU[idx] = - (q[Iu] - qOld[Iu])/dt * 0.5*(dx[I+i]+dx[I-1+i]);
}

__global__
void liftLeftRight(real *FyX, real *q, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy)
{
	int  idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx > ncy)
		return;
	int  Iu = (J+idx)*(nx-1) + (I-1),
	     Iv = (nx-1)*ny + (J-1+idx)*nx + I;
	FyX[idx] = -(
	              // multiply by dS
	              (
	                0.25 * ( q[Iu+ncx]/dy[J+idx] + q[Iu+ncx-(nx-1)]/dy[J-1+idx] )
	                     * ( q[Iv+ncx]/dx[I+ncx] + q[Iv+ncx-1]/dx[I+ncx-1] )
	                -
	                0.25 * ( q[Iu]/dy[J+idx] + q[Iu-(nx-1)]/dy[J-1+idx] )
	                     * ( q[Iv]/dx[I] + q[Iv-1]/dx[I-1] )
	              )
	              -
	              // multiply by dS (cannot use the leftRight trick in this case)
	              nu*
	              (
	                (
	                  (q[Iu+ncx]/dy[J+idx] - q[Iu+ncx-(nx-1)]/dy[J-1+idx])/2.0/(dy[J+idx]+dy[J-1+idx]) +
	                  (q[Iv+ncx]/dx[I+ncx] - q[Iv+ncx-1]/dx[I+ncx-1])/2.0/(dx[I+ncx]+dx[I+ncx-1])
	                ) 
	                -
	                (
	                  (q[Iu]/dy[J+idx] - q[Iu-(nx-1)]/dy[J-1+idx])/2.0/(dy[J+idx]+dy[J-1+idx]) +
	                  (q[Iv]/dx[I] - q[Iv-1]/dx[I-1])/2.0/(dx[I]+dx[I-1])
	                )
	              )
	            )*0.5*(dy[J+idx]+dy[J-1+idx]);
}

__global__
void liftBottomTop(real *FyY, real *q, real *lambda, real nu, real *dx, real *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy)
{
	int  idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx >= ncx)
		return;
	int  Ip = J*nx + I+idx,
	     Iv = (nx-1)*ny + (J-1)*nx + I+idx;
	FyY[idx] = -(
	              // multiply the pressure with the surface area to get p dx
	              (lambda[Ip+ncy*nx]-lambda[Ip-nx])*dx[I+idx]
	              +
	              // divide q^2 by dx, so that just v^2 dx is obtained
	              (
	                  0.25*(q[Iv+(ncy+1)*nx] + q[Iv+ncy*nx])*(q[Iv+(ncy+1)*nx] + q[Iv+ncy*nx])
	                - 0.25*(q[Iv] + q[Iv-nx])*(q[Iv] + q[Iv-nx])
	              )/dx[I+idx]
	              -
	              // no multiplication or division since dv/dy dx = dq/dy
	              nu*
	              (
	                  (q[Iv+(ncy+1)*nx] - q[Iv+ncy*nx])/dy[J+ncy]
	                - (q[Iv] - q[Iv-nx])/dy[J-1]
	              )
	            );
}
            
__global__
void liftUnsteady(real *FyU, real *q, real *qOld, real *dx, real *dy, real dt,
                  int nx, int ny, int I, int J, int ncx, int ncy)
{
	int  idx = threadIdx.x + blockIdx.x*blockDim.x;
	     
	if( idx >= ncx*(ncy+1) )
		return;
	
	int i = idx%ncx,
	    j = idx/ncx;
	    
	int Iv = (J-1+j)*nx + (I+i) + (nx-1)*ny;

	FyU[idx] = - (q[Iv] - qOld[Iv])/dt * 0.5*(dy[J+j]+dy[J-1+j]);
}

__global__
void forceX(real *f, real *q, real *rn, int *tags,
            int nx, int ny, real *dx, real *dy,
            real dt, real alpha, real nu)
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

	int  Gidx_x = J*(nx-1) + I;

	real dTerm;
	
	__shared__ real u[BSZ][BSZ];
						
	__shared__ real Dx[BSZ][BSZ], Dy[BSZ][BSZ];
	
	Dy[j][i] = dy[J];
	Dx[j][i] = dx[I];
	
	/// transfer from global to shared memory
	u[j][i] = q[Gidx_x]/Dy[j][i];
	__syncthreads();
	
	/// check bounds for convective term in the x-direction
	int global_check = ( I==0 || I==(nx-2) || J==0 || J==(ny-1) ),		///< check if we compute globally
	    block_check  = ( i==0 || i==(BSZ-1) || j==0 || j==(BSZ-1) );	///< check if element within block computes
	
	/// X-component
	if( !(global_check || block_check) )
	{
		dTerm = alpha*nu*2.0*( \
						 ( Dx[j][i]*u[j][i+1] - (Dx[j][i]+Dx[j][i+1])*u[j][i] + Dx[j][i+1]*u[j][i-1] ) / ( Dx[j][i]*Dx[j][i+1]*(Dx[j][i]+Dx[j][i+1]) ) \
					   
					   + 4.0*( (Dy[j][i]+Dy[j-1][i])*u[j+1][i] - (Dy[j-1][i] + 2.0*Dy[j][i] + Dy[j+1][i])*u[j][i] + (Dy[j][i]+Dy[j+1][i])*u[j-1][i] ) \
							/( (Dy[j][i]+Dy[j-1][i]) * (Dy[j-1][i] + 2.0*Dy[j][i] + Dy[j+1][i]) * (Dy[j][i]+Dy[j+1][i]) ) \
					     );
		
		f[Gidx_x] = ( u[j][i]/dt - dTerm - rn[Gidx_x]/(0.5*(Dx[j][i]+Dx[j][i+1])) ) * (!(tags[Gidx_x]==-1));
	}
}

__global__
void forceY(){}

} // end of namespace kernels
