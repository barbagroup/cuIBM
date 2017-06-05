"""
Computes, plots, and saves the 2D vorticity field from a cuIBM simulation at
saved time-steps.
"""

from snake.cuibm.simulation import CuIBMSimulation


simulation = CuIBMSimulation()
simulation.read_grid()

for time_step in simulation.get_time_steps():
  simulation.read_fields('vorticity', time_step)
  simulation.plot_contour('vorticity',
                          field_range=(-5.0, 5.0, 101),
                          view=[-1.0, -5.0, 15.0, 5.0],
                          cmap='viridis',
                          style='seaborn-dark',
                          width=8.0)
