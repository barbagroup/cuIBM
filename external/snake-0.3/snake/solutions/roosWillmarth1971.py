"""
Implementation of the class `RoosWillmarth1971` to read the steady-state drag
coefficient from a file.

The drag coefficients have been reported in Roos and Willmarth (1971).

_References:_
* Roos, F. W., & Willmarth, W. W. (1971).
  Some experimental results on sphere and disk drag.
  AIAA Journal, 9(2), 285-291.
"""

import numpy


class RoosWillmarth1971(object):
  """
  Contains info about the experimental data from Roos and Willmarth (1971).
  """

  def __init__(self, file_path=None):
    """
    If provided, read the drag coefficient versus the Reynolds number from
    file.

    Parameters
    ----------
    file_path: string, optional
      Path of the file to read;
      default: None.
    """
    self.description = 'Roos and Willmarth (1971)'
    if file_path:
      self.read_drag_coefficient(file_path)

  def read_drag_coefficient(self, file_path):
    """
    Reads the drag coefficient versus the Reynolds number from given file path.

    Parameters
    ----------
    file_path: string
      Path of the file to read.
    """
    print('[info] reading drag coefficient '
          'from file {} ...'.format(file_path)),
    with open(file_path, 'r') as infile:
      self.re, self.cd = numpy.loadtxt(infile,
                                       dtype=float,
                                       comments='#',
                                       unpack=True)
    print('done')
