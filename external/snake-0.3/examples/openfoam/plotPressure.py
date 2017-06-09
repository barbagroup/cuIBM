"""
Plots the 2D pressure field from a OpenFOAM simulation at given times using the
visualization software ParaView.
"""

from snake.openfoam.simulation import OpenFOAMSimulation


simulation = OpenFOAMSimulation()
simulation.plot_field_contours_paraview('pressure',
                                        field_range=(-1.0, 0.5),
                                        view=(-2.0, -5.0, 15.0, 5.0),
                                        times=(0.0, 100.0, 2.0),
                                        width=800)
