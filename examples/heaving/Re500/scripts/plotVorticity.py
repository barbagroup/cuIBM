"""
Plots the vorticity field of a 2D cuIBM simulation.
"""

import os
from snake.cuibm.simulation import CuIBMSimulation
from snake.body import Body


simulation = CuIBMSimulation()
simulation.read_grid()

for time_step in simulation.get_time_steps():
  body = Body(file_path=os.path.join('{:0>7}'.format(time_step), 'bodies'))
  simulation.read_fields('vorticity', time_step)
  simulation.plot_contour('vorticity',
                          field_range=(-10.0, 10.0, 10),
                          filled_contour=False,
                          bodies=body,
                          view=(-1.0, -1.0, 3.5, 1.0),
                          style='seaborn-dark',
                          colors='grey',
                          colorbar=False,
                          width=8.0)
