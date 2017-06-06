"""
Implementation of the class `Body`.
"""

import numpy


class Body(object):
  """
  Contains information about an immersed body.
  """

  def __init__(self, file_path=None, label=None):
    """
    Reads the body coordinates from given file.

    Only works with two-dimensional bodies for now.

    Parameters
    ----------
    file_path: string, optional
      Path of the file containing the body coordinates.
      default: None.
    label: string, optional
      Label of the body;
      default: None.
    """
    self.label = label
    self.file_path = file_path
    if file_path:
      self.x, self.y = self.read_coordinates(file_path)

  def read_coordinates(self, file_path):
    """
    Reads the coordinates from file.

    Parameters
    ----------
    file_path: string, optional
      Path of the file with the body coordinates.

    Returns
    -------
    coords: numpy ndarray of floats
      The coordinates of the body.
    """
    print('[info] reading body coordinates from file '
          '{} ...'.format(self.file_path))
    with open(self.file_path, 'r') as infile:
      return numpy.loadtxt(infile, dtype=float, skiprows=1, unpack=True)
