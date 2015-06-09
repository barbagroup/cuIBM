
# Re-designing cuIBM:

**Objective**: re-think the design of `cuIBM` to make it more modular, more readable, so that it will be easier to export the code to multi-gpus and to deal with AMR and to move to 3D.

Barbagroup members involved: Pi-Yueh, Gil and Olivier.

---

## General naming rules

* the name of a class starts with an uppercase
* the name of an instance of a class starts with a lower case
* use camelCase convention for variable names
* methods and functions start with a lowercase

---

## Code modification workflow

Everything related to the re-design of `cuIBM` is located in the branch `redesign` of its child branches.

A modification has to be provided, accepted and merged.

To provide a modification, create a branch from base-branch `redesign`.

Modifications are submitted via pull-request. Once the pull-request has been reviewed and accepted by all members listed above, it can be merge into the branch `redesign`.

---

## Modifications proposed

* Remove the namespace `io` and parse the input files directly in the method (called `parseInputFile` for example) of the appropriate class.

* Remove the class `property` and stop storing parameters into the container `database` (or `DB`); instead, store directly the parameters in the appropriate class.

* Remove the CPU functions/methods. This code is solving PDEs on GPU.

* Command-line parameters should overwrite parameters parsed from input files.

---

## Non-exhaustive list of classes with brief description

### Domain
* methods:
  - `parseInputFile`: parses the file `domain.yaml`.
* attributes:
  - list of Cartesian grids, each member being an instance of the class `CartesianGrid`.

### CartesianGrid
* methods:
  - `generate`
* attributes:
  - list of coordinates in each direction.
  - mesh-spacings in each direction.
  - number of cells in each direction.

### Bodies
* methods:
  - `parseInputFile`: parses the file `bodies.yaml` (that contains the path of each body files with coordinates).
* attributes:
  - list of immersed bodies (each body being an instance of the class Body).

### Body
* methods:
  - `readFile`: reads the file containing the coordinates of the body.
  - `prescribedMotion`: implements the kinematics of the body.
* attributes:
  - Lagrangian coordinates.
  - Lagrangian velocity.
  - Lagrangian forces.
  - parameters related to the kinematics.

### FlowDescription
* methods:
  - `parseInputFile`: parses the file `flowDescription.yaml`.
* attributes:
  - viscosity
  - list of boundary condition descriptions (each member being an instance of the class `BoundaryCondition`): xMinus, xPlus, yMinus, yPlus, zMinus, zPlus.

### SimulationParameters
* methods:
  - `parseInputFile`: parses the file `simulationParameters.yaml`
* attributes:
  - `dt`: the time-increment.
  - `nt`: number of time-steps to perform.
  - `start`: the initial time-step.
  - `nsave`: the time-step saving-interval.
  - name of the immersed boundary method to use.
  - temporal discretization of the advective and diffusive term.
  - list of information about the linear solvers (tolerance, maximum number of iterations, type of preconditioner, type of solver).

### BoundaryCondition
* methods:
* attributes:
  - type: e.g. DIRICHLET, NEUMANN, CONVECTIVE, PERIODIC.
  - value: value associated to the boundary condition.

### LinearSolver
* methods:
  - `solve`: solves the linearSystem given the initial solution and the RHS of the system.
  - `setup`: grabs the information about the linear solver that were stored in the list of linear solvers of the class `SimulationParameters`.
  - `generate`: generates the matrix and the preconditioner.
* attributes:
  - type of solver
  - the matrix
  - the preconditioner
  - tolerance
  - maximum number of iterations

### Data
* methods:
  - `read`: reads the field data from a file at a certain time-step.
  - `write`: writes the field data into a file at a certain time-step.
* attributes:
  - field (data structure to be defined)
  - time-step

### NavierStokesSolver
Solves the Navier-Stokes equation using a projection method.

### Perot
Inherits from `NavierStokesSolver`.
Solves the Navier-Stokes equations using a projection method seen as a block-LU decomposition.
* [Reference]((http://www.ecs.umass.edu/mie/tcfd/Papers/FSM0.pdf)): J.B. Perot. An Analysis of the Fractional Step Method. *J. Comp. Phys.* 108, 51-58 (1993)

### TairaColonius
Inherits from `Perot`.
The discrete pressure and the Lagrangian forces are grouped together.
* [Reference](http://colonius.caltech.edu/pdfs/TairaColonius2007.pdf): K. Taira and T. Colonius. The immersed boundary method: A projection approach. *J. Comp. Phys.* 225, 2118-2137 (2007).

### DirectForcing
Inherits from `Perot`.
Fully discrete formulation of the direct forcing method.
At forcing nodes, the discretized momentum equation is replaced by some kind of interpolation.

### FadlunEtAl
Inherits from `DirectForcing`.
Includes an explicit pressure-gradient term in the discretized momentum equation.
The solution of the Poisson system is a pressure correction, that is added to the pressure at the previous time-step to get the new pressure.
The Poisson system is solved on the entire domain, like if there was no immersed boundary in the flow.
* [Reference](http://www.researchgate.net/profile/Paolo_Orlandi/publication/222523772_Combined_Immersed-Boundary_Finite-Difference_Methods_for_Three-Dimensional_Complex_Flow_Simulations/links/540ec5690cf2d8daaacd6e2b.pdf): E.A. Fadlun, R. Verzicco, P. Orlandi and J. Mohd-Yusof. Combined Immersed-Boundary Finite-Difference Methods for Three-Dimensional Complex Flow Simulations. *J. Comp. Phys.* 161, 35-60 (2000).

### Krishnan
Inherits from `DirectForcing`.
Second-order fully-discrete direct forcing method.
The pressure gradient is interpolated at forcing nodes.
This corresponds to imposing Neumann boundary condition on the pressure at the physical location of the immersed body, hence, decoupling the pressure inside and outside the body.
* [Reference](https://github.com/barbagroup/Theses/tree/master/Anush): A. Krishnan. Towards the study of flying snake aerodynamics, and An Analysis of the Direct Forcing Method. Phd thesis (2015).

### Chorin
Inherits from `NavierStokesSolver`.
Time-splitting method.
