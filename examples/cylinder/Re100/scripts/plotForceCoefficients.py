"""
Post-processes the force coefficients from a cuIBM simulation.

This script reads the forces, computes the mean forces within a given range,
computes the Strouhal number within a range, plots the force coefficients,
saves the figure, and prints a data-frame that contains the mean values.
"""

from snake.cuibm.simulation import CuIBMSimulation


simulation = CuIBMSimulation()
simulation.read_forces()
time_limits = (60.0, 80.0)
simulation.get_mean_forces(limits=time_limits)
simulation.get_strouhal(limits=time_limits, order=200)

simulation.plot_forces(display_coefficients=True,
                       coefficient=2.0,
                       limits=(0.0, 80.0, -0.5, 1.5),
                       style='seaborn-dark',
                       save_name='forceCoefficients')

dataframe = simulation.create_dataframe_forces(display_strouhal=True,
                                               display_coefficients=True,
                                               coefficient=2.0)
print(dataframe)
