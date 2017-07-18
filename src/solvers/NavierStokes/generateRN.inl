/**
 * \file generateRN.inl
 * \brief Implementation of the methods to generate 
 *        the explicit terms of the momentum equation.
 */


#include <solvers/kernels/generateRN.h>

#define BSZ 16


/**
 * \brief Doing nothing. Used in methods that use the explicit pressure term 
 *        in the intermediate velocity solve step, such as FadlunEtAlSolver
 *        and DFModifiedSolver.
 */
template <typename memoryType>
void NavierStokesSolver<memoryType>::calculateExplicitLambdaTerms()
{
}


/**
 * \brief Generates explicit terms that arise from the velocity fluxes.
 *
 * Includes the time derivative, convection and diffusion terms.
 */
template <>
void NavierStokesSolver<device_memory>::calculateExplicitQTerms()
{
	real gamma = intgSchm.gamma[subStep],
	     zeta = intgSchm.zeta[subStep],
	     alpha = intgSchm.alphaExplicit[subStep],
	     nu = (*paramDB)["flow"]["nu"].get<real>();
	// if first time-step: use Euler explicit for convective terms
	if (timeStep == (*paramDB)["simulation"]["startStep"].get<int>())
	{
		gamma = 1.0;
		zeta = 0.0;
	}
	     
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
	
	int nx = domInfo->nx,
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
} // calculateExplicitQTerms


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
} // generateRN
