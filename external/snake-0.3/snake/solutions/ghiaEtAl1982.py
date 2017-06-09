"""
Implementation of the class `GhiaEtAl1982` that reads the centerline velocities
reported in Ghia et al. (1982).

_References:_
* Ghia, U. K. N. G., Ghia, K. N., & Shin, C. T. (1982).
  High-Re solutions for incompressible flow using the Navier-Stokes equations
  and a multigrid method.
  Journal of computational physics, 48(3), 387-411.
"""

import os

import numpy


class GhiaEtAl1982(object):
  """
  Container to store results from Ghia et al. (1982).
  """

  def __init__(self, Re=None, file_path=None):
    """
    Initialization.

    Parameters
    ----------
    Re: float, optional
      Desired Reynolds number;
      default: None.
    file_path: string, optional
      Path of the file containing the validation data;
      default: None.
    """
    self.y, self.u = None, None
    self.x, self.v = None, None
    if Re:
      self.read_centerline_velocities(Re, file_path=file_path)

  def read_centerline_velocities(self, Re, file_path=None):
    """
    Reads the centerline velocities from file and for a given Reynolds number.

    Parameters
    ----------
    Re: float
      Desired Reynolds number.
    file_path: string, optional
      Path of the file containing the validation data;
      default: None (will be read the file located in `resources` directory of
      the `snake` package).
    """
    if not file_path:
      file_path = os.path.join(os.environ['SNAKE'],
                               'resources',
                               'results',
                               'ghia_et_al_1982_lid_driven_cavity.dat')
    Re = str(int(round(Re)))
    # column indices in file with experimental results
    cols = {'100': {'u': 1, 'v': 7},
            '1000': {'u': 2, 'v': 8},
            '3200': {'u': 3, 'v': 9},
            '5000': {'u': 4, 'v': 10},
            '10000': {'u': 5, 'v': 11}}

    with open(file_path, 'r') as infile:
      y, u, x, v = numpy.loadtxt(infile,
                                 dtype=float,
                                 usecols=(0, cols[Re]['u'], 6, cols[Re]['v']),
                                 unpack=True)
    self.y, self.u = y, u
    self.x, self.v = x, v
