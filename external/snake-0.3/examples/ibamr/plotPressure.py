"""
Plots the 2D pressure field from a IBAMR simulation at saved time-steps using
the visualization software VisIt.
"""

from snake.ibamr.simulation import IBAMRSimulation


simulation = IBAMRSimulation()

body_name = 'flyingSnake2dAoA35ds004filledInside'  # file name (no extension)
simulation.plot_field_contours_visit('pressure', (-1.0, 0.5),
                                     body=body_name,
                                     solution_folder='numericalSolution',
                                     view=(-15.0, -15.0, 15.0, 15.0),
                                     width=800)
