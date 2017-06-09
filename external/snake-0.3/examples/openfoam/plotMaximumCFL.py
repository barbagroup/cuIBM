"""
Plots the maximum CFL number over the course of the simulation.

The maximum CFL values are read from the log file (output) of OpenFOAM.
"""

import os
from snake.openfoam.simulation import OpenFOAMSimulation


simulation = OpenFOAMSimulation()
filepath = os.path.join('log.run', 'log.icoFoam')
simulation.read_maximum_cfl(filepath)
simulation.get_mean_maximum_cfl(limits=(60.0, 80.0))
simulation.plot_maximum_cfl(display_extrema=True, order=200,
                            limits=(0.0, 100.0, 0.5, 2.0),
                            save_name='maximumCFL')
