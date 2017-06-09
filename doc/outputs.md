All cuIBM output files created during the course of the simulation are located in the simulation directory and its sub-folders. Here is a list of such files:

* `grid`: Binary file containing the points along a gridline of the structured Cartesian mesh in each direction. The first line stores the number of cells in each direction. Then, coordinates are stored in ascending order along each direction, starting with the x-direction, immediately followed by the y-direction.

* `forces`: this files stores the hydrodynamic forces acting on the immersed body at each time-step. The first column contains the time values. Following columns contains the force component in the x- and y-directions. This file is only generated when a body is immersed in the computational domain.

* `iterations`: this file consists of three columns - the first columns contains the time-step index; the second and third columns contain the number of iterations at a time-step required by the iterative solvers to converge (second column for the velocity solver, third one for the Poisson solver).

* Every given time-step interval (see `nsave` in the input file `simParams.yaml`), the numerical solution if saved in a sub-folder whose name is the time-step index. Each sub-folder contains the flux solution (binary file `q`); the pressure field--and the Lagrangian forces when using the method from Taira and Colonius (binary file `lambda`); and the coordinates the immersed boundaries (ascii file `bodies`).
