The file `flow.yaml` is required to prescribe the characteristics of the fluid, as well as the initial and boundary conditions of the flow.

Here, we provide an example of such a file followed by an explanation of its content.

Note: text is case-sensitive.


## Example
The following file could be used to simulate external flows over bodies:

    - nu: 0.025
      initialVelocity: [1.0, 0.0]
      initialPerturbation: [0.0, 0.0]
      boundaryConditions:
        - location: xMinus
          u: [DIRICHLET, 1.0]
          v: [DIRICHLET, 0.0]
        - location: xPlus
          u: [CONVECTIVE, 1.0]
          v: [CONVECTIVE, 0.0]
        - location: yMinus
          u: [DIRICHLET, 1.0]
          v: [DIRICHLET, 0.0]
        - location: yPlus
          u: [DIRICHLET, 1.0]
          v: [DIRICHLET, 0.0]


## File options
* `nu`: (mandatory) the kinematic viscosity of the fluid. Can be any positive number.
* `initialVelocity`: (mandatory) this is the initial velocity of the fluid set throughout the computational domain. For 2d flows, two components are required within square brackets, separated by a comma.
* `initialPerturbation`: (optional, default: `[0.0, 0.0]`) initial sinusoidal perturbation in the x- and y- directions.
* `boundaryConditions`: (mandatory) specifies the velocity boundary condition on each edge of the domain. These are listed in subsections, which need to be indented using two spaces.
  * Each boundary edge is represented by its relative `location` on the Cartesian axes - `xMinus`, `xPlus`, `yMinus` and `yPlus`.
  * For each component of velocity (`u` and `v`) on the boundary, you can specify the type of boundary condition and an associated value. There are four types currently available: `DIRICHLET`, `NEUMANN`, and `CONVECTIVE`.
  * The number following the type of boundary conditions is interpreted as:
    * for a `DIRICHLET` type: the value of that component of velocity at the boundary.
    * for a `NEUMANN` type: the value of the normal derivative of that component at the boundary.
    * for a `CONVECTIVE` type: the speed at which the fluid is convected out of the boundary. The values for the components normal to the boundary are ignored.
  * Specific types of boundary conditions can be created by using a suitable combination of the above. For example, a slip boundary condition on the `yMinus` boundary would use homogeneous Neumann for `u` and zero Dirichlet for `v`.
* In the above 2d example, the kinematic viscosity of the fluid is set to `0.025`. `xMinus` is the inlet, `yMinus` and `yPlus` are the outer boundaries and `xPlus` represents the outlet. Dirichlet boundary conditions are used on each edge of the computational domain (velocity fixed to `(1.0, 0.0)`), except at the outlet where a convective condition is used (the flow is convected in the x-direction with speed `1.0`). The initial velocity field is set to `(1.0, 0.0)` everywhere in the domain.
