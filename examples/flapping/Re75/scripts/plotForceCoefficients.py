"""
Plots the force coefficients and compares the solution with results from
Li et al. (2015), Wang et al. (2004), and Eldredge (2007).

_References:_
* Wang, Z. J., Birch, J. M., & Dickinson, M. H. (2004).
  Unsteady forces and flows in low Reynolds number hovering flight:
  two-dimensional computations vs robotic wing experiments.
  Journal of Experimental Biology, 207(3), 449-460.
* Eldredge, J. D. (2007).
  Numerical simulation of the fluid dynamics of 2D rigid body motion with the
  vortex particle method.
  Journal of Computational Physics, 221(2), 626-648.
* Li, C., Dong, H., & Liu, G. (2015).
  Effects of a dynamic trailing-edge flap on the aerodynamic performance and
  flow structures in hovering flight.
  Journal of Fluids and Structures, 58, 49-65.
"""

import os
import math
import numpy
from matplotlib import pyplot


# Read forces from file
filepath = 'forces'
with open(filepath, 'r') as infile:
  data = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)

# Simulation parameters
c = 1.0  # chord-length
A0 = 2.8 * c  # stroke amplitude
f = 0.25  # flapping frequency
rho = 1.0  # density
Umax = A0 * numpy.pi * f  # maximum translational velocity
Re = 75.0
nu = Umax * c / Re
print('Kinematic viscosity: {}'.format(nu))

# Scale the time by the frequency to get t/T
times = f * data[0]
# Scale the forces to get the force coefficients.
# Forces are normalized by the maxima of the quasi-steady forces
# as described in Wang et al. (2004).
cd = data[1] / (0.5 * rho * Umax**2 * c) / 1.4
cl = data[2] / (0.5 * rho * Umax**2 * c) / 1.2
# Handle stroke reversal
for i in range(cd.size):
  if (times[i] - math.floor(times[i]) > 0.5):
    cd[i] *= -1

# Apply smoothing to reduce noise
# N = 10  # length of the window
# window = numpy.ones(N) / N
# cl = numpy.convolve(cl, window, 'same')
# cd = numpy.convolve(cd, window, 'same')

directory = os.path.join(os.environ['CUIBM_DIR'], 'data')
# Load lift coefficient from Li et al. (2015)
filepath = os.path.join(directory,
                        'li_et_al_2015_flapping_liftCoefficientRe75.dat')
with open(filepath, 'r') as infile:
  cl_li = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)
# Load drag coefficient from Li et al. (2015)
filepath = os.path.join(directory,
                        'li_et_al_2015_flapping_dragCoefficientRe75.dat')
with open(filepath, 'r') as infile:
  cd_li = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)
# Load lift coefficient from Wang et al. (2004) digitized by Li et al. (2015)
filepath = os.path.join(directory,
                        'wang_et_al_2004_flapping_liftCoefficientRe75.dat')
with open(filepath, 'r') as infile:
  cl_wang = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)
# Load drag coefficient from Wang et al. (2004) digitized by Li et al. (2015)
filepath = os.path.join(directory,
                        'wang_et_al_2004_flapping_dragCoefficientRe75.dat')
with open(filepath, 'r') as infile:
  cd_wang = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)
# Load lift coefficient from Eldredge (2007) digitized by Li et al. (2015)
filepath = os.path.join(directory,
                        'eldredge_2007_flapping_liftCoefficientRe75.dat')
with open(filepath, 'r') as infile:
  cl_eldredge = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)
# Load drag coefficient from Eldredge (2007) digitized by Li et al. (2015)
filepath = os.path.join(directory,
                        'eldredge_2007_flapping_dragCoefficientRe75.dat')
with open(filepath, 'r') as infile:
  cd_eldredge = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)

# Plot the force coefficients and compare to results from the literature
pyplot.style.use('seaborn-dark')
fig, ax = pyplot.subplots(2, figsize=(10.0, 6.0), sharex=True)
ax[0].grid()
ax[0].set_ylabel('$C_L$', fontsize=16)
ax[0].plot(times, cl,
           label='cuIBM',
           color='#f03b20', linewidth=1)
ax[0].plot(cl_li[0], cl_li[1],
           label='Li et al. (2015)',
           color='#004529', linestyle='-', linewidth=2, zorder=10)
ax[0].plot(cl_eldredge[0], cl_eldredge[1],
           label='Eldredge (2007)',
           color='#252525', linestyle='--', linewidth=2, zorder=11)
ax[0].scatter(cl_wang[0], cl_wang[1],
              label='Wang et al. (2004)',
              marker='o', facecolors='none', edgecolors='#2c7fb8', linewidth=2,
              zorder=12)
ax[0].set_xlim(0.0, 4.0)
ax[0].set_ylim(-0.25, 1.5)
ax[1].grid()
ax[1].set_xlabel('t/T', fontsize=16)
ax[1].set_ylabel('$C_D$', fontsize=16)
ax[1].plot(times, cd,
           label='cuIBM',
           color='#f03b20', linewidth=1)
ax[1].plot(cd_li[0], cd_li[1],
           label='Li et al. (2015)',
           color='#004529', linestyle='-', linewidth=2, zorder=10)
ax[1].plot(cd_eldredge[0], cd_eldredge[1],
           label='Eldredge (2007)',
           color='#252525', linestyle='--', linewidth=2, zorder=11)
ax[1].scatter(cd_wang[0], cd_wang[1],
              label='Wang et al. (2004)',
              marker='o', facecolors='none', edgecolors='#2c7fb8', linewidth=2,
              zorder=12)
ax[1].set_xlim(0.0, 4.0)
ax[1].set_ylim(-2.0, 2.0)

handles, labels = ax[1].get_legend_handles_labels()
fig.legend(handles, labels,
           ncol=2, loc='center', prop={'size': 12}, frameon=False,
           bbox_to_anchor=(0.5, 0.95))

# Save figure as PNG
directory = os.path.join(os.getcwd(), 'images')
if not os.path.isdir(directory):
  os.makedirs(directory)
filepath = os.path.join(directory, 'forceCoefficients.png')
pyplot.savefig(filepath, dpi=100)
pyplot.show()
