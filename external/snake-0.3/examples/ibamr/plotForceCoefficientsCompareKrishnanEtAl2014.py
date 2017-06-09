"""
Post-processes the force coefficients from a IBAMR simulation and compare them
to cuIBM ones (obtained by Anush Krishnan and
published into Krishnan et al., 2014).

This script reads the forces, computes the mean forces within a given range,
computes the Strouhal number within a range, plots the force coefficients,
saves the figure, and prints a data-frame that contains the mean values.
"""

import os

from snake.ibamr.simulation import IBAMRSimulation
from snake.cuibm.simulation import CuIBMSimulation


simulation = IBAMRSimulation(description='IBAMR')
simulation.read_forces()
time_limits = (32.0, 64.0)
simulation.get_mean_forces(limits=time_limits)
simulation.get_strouhal(limits=time_limits, order=200)

krishnan = CuIBMSimulation(description='Krishnan et al. (2014)')
filepath = os.path.join(os.environ['SNAKE'],
                        'resources',
                        'flyingSnake2d_cuibm_anush',
                        'flyingSnake2dRe2000AoA35',
                        'forces')
krishnan.read_forces(file_path=filepath)
krishnan.get_mean_forces(limits=time_limits)
krishnan.get_strouhal(limits=time_limits, order=200)

simulation.plot_forces(display_coefficients=True,
                       coefficient=-2.0,
                       display_extrema=True, order=200,
                       limits=(0.0, 80.0, 0.0, 3.0),
                       other_simulations=krishnan,
                       other_coefficients=2.0,
                       style='mesnardo',
                       save_name='forceCoefficientsCompareKrishnanEtAl2014')

dataframe = simulation.create_dataframe_forces(display_strouhal=True,
                                               display_coefficients=True,
                                               coefficient=-2.0)
dataframe2 = krishnan.create_dataframe_forces(display_strouhal=True,
                                              display_coefficients=True,
                                              coefficient=2.0)
print(dataframe.append(dataframe2))
