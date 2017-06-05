"""
Plots the vorticity field of a 2D cuIBM simulation.
"""

import os
from snake.cuibm.simulation import CuIBMSimulation
from snake.body import Body


simulation = CuIBMSimulation()
simulation.read_grid()

for time_step in simulation.get_time_steps():
  all_bodies = Body(file_path=os.path.join('{:0>7}'.format(time_step),
                                           'bodies'))
  n_total = all_bodies.x.size
  bodies = [Body(), Body()]
  bodies[0].x = all_bodies.x[:n_total // 2]
  bodies[0].y = all_bodies.y[:n_total // 2]
  bodies[1].x = all_bodies.x[n_total // 2:]
  bodies[1].y = all_bodies.y[n_total // 2:]
  simulation.read_fields('vorticity', time_step)
  simulation.plot_contour('vorticity',
                          field_range=[-2.0, 2.0, 40],
                          filled_contour=True,
                          bodies=bodies,
                          view=[-3.0, -5.0, 15.0, 5.0],
                          style='seaborn-dark',
                          cmap='viridis',
                          width=8.0)
