
### Re-designing cuIBM:
-removing io namespace, removing property class

Domain
-parse(domain.yaml)
-list(Mesh)

Mesh
-generate
x, y, z, spacings

Bodies
-parse(bodies.yaml)
-list(Body)

Body
- X, Y, Z, spacings
- kinematics, prescribed motion

FlowDescription
(just parameters)
- parse(flowDescription.yaml)
- list(BoundaryCondition)
xPlus, xMinus, same for y and z.
viscosity

SimulationParameters
- parse(simulationParameters.yaml)

Data
- input and output methods
- data structure
- time-step

* command-line parameters should overwrite database parameters (read from .yaml files).

NavierStokesSolver

