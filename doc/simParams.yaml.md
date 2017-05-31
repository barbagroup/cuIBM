The file `simParams.yaml` is required to provide information to the different solvers.

Here, we provide an example of such a file followed by an explanation of its content.

Note: text is case-sensitive.


## Example

    - dt: 0.01
      scaleCV: 5.0
      startStep: 0
      nt: 300
      nsave: 50
      ibmScheme: TAIRA_COLONIUS
      timeScheme: [ADAMS_BASHFORTH_2, CRANK_NICOLSON]
      interpolationType: LINEAR
      linearSolvers:
        - system: velocity
          solver: BICGSTAB
          preconditioner: DIAGONAL
          absTolerance: 1.0E-06
          relTolerance: 0.0
          maxIterations: 10000
        - system: Poisson
          solver: CG
          preconditioner: SMOOTHED_AGGREGATION
          absTolerance: 1.0E-06
          relTolerance: 0.0
          maxIterations: 20000


## File options

* `dt`: (mandatory) the time-increment. Can be any number greater than zero (as long as it satisfies the stability criterion for the numerical scheme used).
* `scaleCV`: (optional, default: `2.0`) parameter to help defining the box surrounding the immersed boundary.
* `startStep`: (optional, default: `0`) starting time-step. If different than the default value `0`, then the program reads the numerical solution of the stating time-step given. This parameters should be used to restart a simulation.
* `nt`: (mandatory) number of time-steps to execute.
* `nsave`: (mandatory) time-step interval at which Eulerian and Lagrangian quantities are saved into files.
* `ibmScheme`: (optional) specifies the immersed boundary method used in the simulation. There are several immersed-boundary methods implemented in cuIBM.
  - `TAIRA_COLONIUS` is an immersed-boundary projection method where the pressure field and the Lagrangian forces are coupled together and a modified Poisson system is solved at each time step.
  - `FADLUN_ET_AL` implements a direct-forcing method proposed by Fadlun et al. (2000), a pressure-correction scheme where the no-slip condition is enforced on the intermediate velocity field.
  - `DIRECT_FORCING` is based on the method first proposed by Fadlun et al (2000) with modifications by Kim et al (2001). It does not follow the same equations that they used, but use a fractional step method starting with the discretized equations, based on the idea by Perot (1993). This method does not use an explicit pressure term in the step where the intermediate velocity is calculated, and the pressure is directly obtained from the Poisson equation.
  - `DIFFUSION` solves unsteady diffusion equation with an immersed boundary using the direct forcing method.
  - `DF_IMPROVED` is a second-order fully-discrete direct forcing method.
* `timeScheme`: (optional, default: `[EULER_EXPLICIT, EULER_IMPLICIT]`) specifies the time-schemes to use for the convective terms and diffusive terms of the momentum equation. In cuIBM, the convective terms can be temporally discretized using an explicit Euler method (`EULER_EXPLICIT`, default value) or a second-order Adams-Bashforth scheme (`ADAMS_BASHFORTH_2`). The diffusive terms can be  treated explicitly (`EULER_EXPLICIT`), implicitly (`EULER_IMPLICIT`, default), or using a second-order Crank-Nicolson scheme (`CRANK_NICOLSON`).
* For the linear systems, the user can control the type of preconditioner to use, the type of solver, and the exit criterion of the solver.
  - `preconditioner`: (optional, default: `DIAGONAL`) the type of preconditioner. Choices are: `NONE`, `DIAGONAL`, `SMOOTH_AGGREGATION`, and `AINV`.
  - `solver`: (optional, default: `CG`) Krylov solver to use. Currently, conjugate-gradient (`CG`), bi-conjugate gradient stabilized (`BICGSTAB`), and GMRES (`GMRES`) are supported.
  - `restart`: (optional, default: `50`) restart value for GMRES.
  - `absTolerance`: (optional, default: `1.0E-50`) absolute tolerance in the reduction of the residuals as an exit criterion.
  - `relTolerance`: (optional, default: `1.0E-05`) relative tolerance in the reduction of the residuals as an exit criterion.
  - `maxIterations`: (optional, default: `10000`) maximum number of iterations allowed for the solver.
