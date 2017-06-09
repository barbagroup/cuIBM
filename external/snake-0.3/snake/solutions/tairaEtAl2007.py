"""
Implementation of the class `TairaEtAl2007` that contains the
steady-state force coefficients of an inclined flat-plate with aspect-ratio 2
at Reynolds number 100 for angles of attack between 0 and 90 degrees.

The force coefficient data files are located in the folder
`resources/results` of the snake package.

_References:_
* Taira, K., Dickson, W. B., Colonius,
  T., Dickinson, M. H., & Rowley, C. W. (2007).
  Unsteadiness in flow over a flat plate at angle-of-attack at low Reynolds
  numbers.
  AIAA Paper, 710, 2007.
"""

import numpy


class TairaEtAl2007(object):
  """
  Contains experimental force coefficients for various angles-of-attack
  of a flat-plate with aspect-ratio 2.

  The experimental results are reported in Taira et al. (2007)
  and Taira (2008).
  """

  def __init__(self):
    self.cd = None
    self.cl = None

  def read_drag_coefficients(self, file_path):
    """
    Reads the angles-of-attack and the drag coefficients from file.

    Parameters
    ----------
    file_path: string
      Path of the file containing the angles and drag coefficients.
    """
    print('[info] reading drag coefficients '
          'from file {} ...'.format(file_path))
    with open(file_path, 'r') as infile:
      angles, values = numpy.loadtxt(infile,
                                     dtype=float,
                                     comments='#',
                                     unpack=True)
    self.cd = [angles, values]

  def read_lift_coefficients(self, file_path):
    """
    Reads the angles-of-attack and the lift coefficients from file.

    Parameters
    ----------
    file_path: string
      Path of the file containing the angles and lift coefficients.
    """
    print('[info] reading lift coefficients '
          'from file {} ...'.format(file_path))
    with open(file_path, 'r') as infile:
      angles, values = numpy.loadtxt(infile,
                                     dtype=float,
                                     comments='#',
                                     unpack=True)
    self.cl = [angles, values]
