"""
Plots the 2D vorticity field from a IBAMR simulation at saved time-steps using
the visualization software VisIt.
"""

from snake.ibamr.simulation import IBAMRSimulation


simulation = IBAMRSimulation()

body_name = 'flyingSnake2dAoA35ds004filledInside'  # file name (no extension)
simulation.plot_field_contours_visit('vorticity', (-5.0, 5.0),
                                     body=body_name,
                                     solution_folder='numericalSolution',
                                     view=(-2.0, -5.0, 15.0, 5.0),
                                     width=800)
