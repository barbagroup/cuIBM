/***************************************************************************//**
 * \file generateRN.inl
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods to generate 
 *        the explicit terms of the momentum equation.
 */


#include <solvers/NavierStokes/kernels/generateRN.h>

#define BSZ 16


/**
 * \brief Doing nothing. Used in methods that use the explicit pressure term 
 *        in the intermediate velocity solve step, such as FadlunEtAlSolver
 *        and DFModifiedSolver.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::calculateExplicitLambdaTerms()
{
	/**
	* Figure this part out for RK3. Not required for 1-step schemes
	*/
/*	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	// rn = rn - zeta * Q.lambda
	if(fabs(intgSchm.zeta[subStep]) > 1e-6)
	{
		// temp = Q.lambda
		cusp::multiply(Q, lambda, temp1);
		// temp = zeta*temp
		cusp::blas::scal(temp1, intgSchm.zeta[subStep]);
		// rn = rn - temp
		cusp::blas::axpy(temp1, rn, -1.0);
	}
*/
}

/**
 * \brief Generates explicit terms that arise from the velocity fluxes (on the device).
 *        Includes the time derivative, convection and diffusion terms.
 */
template <>
void NavierStokesSolver<device_memory>::calculateExplicitQTerms()
{
	real gamma = intgSchm.gamma[subStep],
	     zeta = intgSchm.zeta[subStep],
	     alpha = intgSchm.alphaExplicit[subStep],
	     nu = (*paramDB)["flow"]["nu"].get<real>();
	     
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
	
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	
	//const int blockEdge = 16;
	
	dim3 dimGridx( int( (nx-1-0.5)/(BSZ-2) ) + 1, int( (ny-0.5)/(BSZ-2) ) + 1 ),
	     dimGridy( int( (nx-0.5)/(BSZ-2) ) + 1, int( (ny-1-0.5)/(BSZ-2) ) + 1 );
	dim3 dimBlock(BSZ, BSZ);
	
	// call the kernel

	// convection terms for the interior points
	kernels::convectionTermU <<<dimGridx, dimBlock>>> (rn_r, H_r, q_r, nx, ny, dxD, dyD, dt, gamma, zeta, alpha, nu);
	kernels::convectionTermV <<<dimGridy, dimBlock>>> (rn_r, H_r, q_r, nx, ny, dxD, dyD, dt, gamma, zeta, alpha, nu);
	
	dim3 dimGridbc(int((nx+ny-0.5)/(BSZ*BSZ))+1, 1);
	dim3 dimBlockbc(BSZ*BSZ, 1);
	
	// calculate convection terms for the rows adjoining the top and bottom boundaries
	kernels::convectionTermUBottomTop <<<dimGridbc, dimBlockbc>>> (rn_r, H_r, q_r, nx, ny, dxD, dyD, dt, gamma, zeta, alpha, nu, yminus, yplus, xminus, xplus);
	kernels::convectionTermVBottomTop <<<dimGridbc, dimBlockbc>>> (rn_r, H_r, q_r, nx, ny, dxD, dyD, dt, gamma, zeta, alpha, nu, yminus, yplus);

	// calculate convection terms for the columns adjoining the left and right boundaries
	kernels::convectionTermULeftRight <<<dimGridbc, dimBlockbc>>> (rn_r, H_r, q_r, nx, ny, dxD, dyD, dt, gamma, zeta, alpha, nu, xminus, xplus);
	kernels::convectionTermVLeftRight <<<dimGridbc, dimBlockbc>>> (rn_r, H_r, q_r, nx, ny, dxD, dyD, dt, gamma, zeta, alpha, nu, yminus, yplus, xminus, xplus);
}

/**
 * \brief Generates explicit terms that arise from the velocity fluxes (on the host).
 *        Includes the time derivative, convection and diffusion terms.
 *        Currently incomplete. Need to fill in the diffusion terms.
 */
template <>
void NavierStokesSolver<host_memory>::calculateExplicitQTerms()
{
	real gamma = intgSchm.gamma[subStep],
	     zeta = intgSchm.zeta[subStep],
	     alpha = intgSchm.alphaExplicit[subStep];
	      
	int  nx = domInfo->nx,
	     ny = domInfo->ny;
	     
	int  numU = (nx-1)*ny;
	int  Iu = 0, Iv = 0;
	real east = 0, west = 0, north = 0, south = 0, Hn = 0, cTerm = 0, dTerm = 0, u = 0, v = 0;

	real *dx = thrust::raw_pointer_cast(&(domInfo->dx[0])),
	     *dy = thrust::raw_pointer_cast(&(domInfo->dy[0]));
	     
	real *xminus = thrust::raw_pointer_cast(&(bc[XMINUS][0])),
	     *xplus  = thrust::raw_pointer_cast(&(bc[XPLUS][0])),
	     *yminus = thrust::raw_pointer_cast(&(bc[YMINUS][0])),
	     *yplus  = thrust::raw_pointer_cast(&(bc[YPLUS][0]));
	     
	real dt = (*paramDB)["simulation"]["dt"].get<real>();
	
	for(int j=0; j<ny; j++)
	{
		// convection terms for the x-momentum equation
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
				
			H[Iu]  = -(east-west)/(0.5*(dx[i]+dx[i+1])) -(north-south)/dy[j];
			cTerm = gamma*H[Iu] + zeta*Hn;
			dTerm = alpha*0;
			rn[Iu] = (u/dt + cTerm + dTerm) * 0.5*(dx[i]+dx[i+1]);
		}
	}
	
	// convection terms for the y-momentum equation
	for(int j=0; j<ny-1; j++)
	{
		for(int i=0; i<nx; i++)
		{
			Iv = j*nx+i + numU;
			Iu = j*(nx-1) + i;
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
				
			H[Iv]  = -(east-west)/dx[i] -(north-south)/(0.5*(dy[j]+dy[j+1]));
			cTerm = gamma*H[Iv] +zeta*Hn;
			dTerm = alpha*0;
			rn[Iv] = (v/dt + cTerm + dTerm) * 0.5*(dy[j]+dy[j+1]);
		}
	}
}

/**
 * \brief Generates explicit terms of the momentum equation.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::generateRN()
{
	logger.startTimer("generateRN");
	
	calculateExplicitQTerms();
	calculateExplicitLambdaTerms();
	
	logger.stopTimer("generateRN");
}
