"""
Post-processes the force coefficients from a IBAMR simulation.

This script reads the forces, computes the mean forces within a given range,
computes the Strouhal number within a range, plots the force coefficients,
saves the figure, and prints a data-frame that contains the mean values.
"""

from snake.ibamr.simulation import IBAMRSimulation


simulation = IBAMRSimulation()
simulation.read_forces()
time_limits = (32.0, 64.0)
simulation.get_mean_forces(limits=time_limits)
simulation.get_strouhal(limits=time_limits, order=200)

simulation.plot_forces(display_coefficients=True,
                       coefficient=-2.0,
                       display_extrema=True, order=200,
                       limits=(0.0, 80.0, 0.0, 3.0),
                       style='mesnardo',
                       save_name='forceCoefficients')

dataframe = simulation.create_dataframe_forces(display_strouhal=True,
                                               display_coefficients=True,
                                               coefficient=-2.0)
print(dataframe)
