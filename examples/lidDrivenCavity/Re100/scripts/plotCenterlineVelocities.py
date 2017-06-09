"""
Plots the velocities along the centerlines of the 2D cavity at Reynolds number
100 and compares with the numerical data reported in Ghia et al. (1982).

_References:_
* Ghia, U. K. N. G., Ghia, K. N., & Shin, C. T. (1982).
  High-Re solutions for incompressible flow using the Navier-Stokes equations
  and a multigrid method.
  Journal of computational physics, 48(3), 387-411.
"""

import os

from matplotlib import pyplot

from snake.cuibm.simulation import CuIBMSimulation
from snake.solutions.ghiaEtAl1982 import GhiaEtAl1982


directory = os.getcwd()

# Reads the velocity fields from files.
simulation = CuIBMSimulation(directory=directory,
                             description='cuIBM')
simulation.read_grid(file_path=os.path.join(directory, 'grid'))
time_step = simulation.get_time_steps()[-1]
simulation.read_fields(['x-velocity', 'y-velocity'], time_step)
# Grabs the mid-cavity velocity values.
y, u = simulation.fields['x-velocity'].get_vertical_gridline_values(0.5)
x, v = simulation.fields['y-velocity'].get_horizontal_gridline_values(0.5)

# Gets the centerline velocities at Re=100 reported in Ghia et al. (1982).
file_path = os.path.join(os.environ['CUIBM_DIR'],
                         'data',
                         'ghia_et_al_1982_lid_driven_cavity.dat')
ghia = GhiaEtAl1982(file_path=file_path, Re=100)

# Plots the instantaneous drag coefficients.
images_directory = os.path.join(directory, 'images')
if not os.path.isdir(images_directory):
  os.makedirs(images_directory)
pyplot.style.use('seaborn-dark')
simu_kwargs = {'label': 'cuIBM',
               'color': '#336699', 'linestyle': '-', 'linewidth': 3,
               'zorder': 10}
ghia_kwargs = {'label': 'Ghia et al. (1982)',
               'color': '#993333', 'linewidth': 0,
               'markeredgewidth': 2, 'markeredgecolor': '#993333',
               'markerfacecolor': 'none',
               'marker': 'o', 'markersize': 8,
               'zorder': 10}

# Plots the x-velocity along a vertical gridline a mid-cavity.
fig, ax = pyplot.subplots(figsize=(6, 6))
ax.grid(True, zorder=0)
ax.set_xlabel('y-coordinate', fontsize=16)
ax.set_ylabel('x-velocity', fontsize=16)
ax.plot(y, u, **simu_kwargs)
ax.plot(ghia.y, ghia.u, **ghia_kwargs)
ax.axis([0.0, 1.0, -0.75, 1.25])
ax.legend(prop={'size': 16})
pyplot.savefig(os.path.join(images_directory,
                            'uVerticalGridline{:0>7}.png'.format(time_step)))

# Plots the y-velocity along an horizontal gridline a mid-cavity.
fig, ax = pyplot.subplots(figsize=(6, 6))
ax.grid(True, zorder=0)
ax.set_xlabel('x-coordinate', fontsize=16)
ax.set_ylabel('y-velocity', fontsize=16)
ax.plot(x, v, **simu_kwargs)
ax.plot(ghia.x, ghia.v, **ghia_kwargs)
ax.axis([0.0, 1.0, -0.75, 1.25])
ax.legend(prop={'size': 16})
pyplot.savefig(os.path.join(images_directory,
                            'vHorizontalGridline{:0>7}.png'.format(time_step)))
