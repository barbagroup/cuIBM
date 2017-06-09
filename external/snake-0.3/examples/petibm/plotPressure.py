"""
Plots and saves the 2D pressure field from a PetIBM simulation at saved
time-steps.
"""

from snake.petibm.simulation import PetIBMSimulation


simulation = PetIBMSimulation()
simulation.read_grid()

for time_step in simulation.get_time_steps():
  simulation.read_fields('pressure', time_step)
  simulation.plot_contour('pressure',
                          field_range=(-1.0, 0.5, 101),
                          filled_contour=True,
                          view=[-15.0, -15.0, 15.0, 15.0],
                          style='mesnardo',
                          cmap='viridis',
                          save_name='pressure',
                          width=8.0)
