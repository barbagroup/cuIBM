"""
Plots the 2D velocity magnitude from a IBAMR simulation at saved time-steps
using the visualization software VisIt.
"""

from snake.ibamr.simulation import IBAMRSimulation


simulation = IBAMRSimulation()

body_name = 'flyingSnake2dAoA35ds004filledInside'  # file name (no extension)
simulation.plot_field_contours_visit('velocity-magnitude', (0.0, 2.0),
                                     body=body_name,
                                     solution_folder='numericalSolution',
                                     view=(-0.75, -1.0, 1.50, 1.0),
                                     width=800)
