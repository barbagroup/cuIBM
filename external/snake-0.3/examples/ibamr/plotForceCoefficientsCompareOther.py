"""
Post-processes the force coefficients from a PetIBM simulation and compare them
to another simulation.

This script reads the forces, computes the mean forces within a given range,
computes the Strouhal number within a range, plots the force coefficients,
saves the figure, and prints a data-frame that contains the mean values.
"""

from snake.simulation import Simulation
from snake.ibamr.simulation import IBAMRSimulation


simulation = IBAMRSimulation(description='IBAMR')
simulation.read_forces()
time_limits = (32.0, 64.0)
simulation.get_mean_forces(limits=time_limits)
simulation.get_strouhal(limits=time_limits, order=200)

other = Simulation(description='',
                   directory='',
                   software='')
other.read_forces()
other.get_mean_forces(limits=time_limits)
other.get_strouhal(limits=time_limits, order=200)

simulation.plot_forces(display_coefficients=True,
                       coefficient=-2.0,
                       display_extrema=True, order=200,
                       limits=(0.0, 80.0, 0.0, 3.0),
                       other_simulations=other,
                       other_coefficients=2.0,
                       style='mesnardo',
                       save_name='forceCoefficientsCompareOther')

dataframe = simulation.create_dataframe_forces(display_strouhal=True,
                                               display_coefficients=True,
                                               coefficient=-2.0)
dataframe2 = other.create_dataframe_forces(display_strouhal=True,
                                           display_coefficients=True,
                                           coefficient=2.0)
print(dataframe.append(dataframe2))
