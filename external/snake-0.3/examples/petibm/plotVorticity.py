"""
Computes, plots, and saves the 2D vorticity field from a PetIBM simulation at
saved time-steps.
"""

from snake.petibm.simulation import PetIBMSimulation


simulation = PetIBMSimulation()
simulation.read_grid()

for time_step in simulation.get_time_steps():
  simulation.read_fields('vorticity', time_step)
  simulation.plot_contour('vorticity',
                          field_range=(-5.0, 5.0, 101),
                          filled_contour=True,
                          view=[-2.0, -5.0, 15.0, 5.0],
                          style='mesnardo',
                          width=8.0)
