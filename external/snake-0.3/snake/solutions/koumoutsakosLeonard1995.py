"""
Implementation of the class `KoumoutsakosLeonard1995` that contains the
instantaneous drag coefficients of an impulsively-started 2D cylinder for
Reynolds numbers 40, 550, and 3000.

The drag coefficient data files are located in the folder
`resources/results` of the snake package.

_References:_
* Koumoutsakos, P., & Leonard, A. (1995).
  High-resolution simulations of the flow around an impulsively started
  cylinder using vortex methods.
  Journal of Fluid Mechanics, 296, 1-38.
"""

import os

import numpy

from snake.force import Force


class KoumoutsakosLeonard1995(object):
  """
  Container to store results from Koumoutsakos and Leonard (1995).
  """

  def __init__(self, file_path=None, Re=None):
    """
    Initializes.

    Parameters
    ----------
    file_path: string, optional
      Path of the file containing the instantaneous drag coefficients;
      default: None.
    Re: float, optional
      Reynolds number;
      default: None.
    """
    self.description = 'Koumoutsakos and Leonard (1995)'
    self.cd = None
    if file_path or Re:
      self.read_drag(file_path=file_path, Re=Re)

  def read_drag(self, file_path=None, Re=None):
    """
    Reads the instantaneous drag coefficients from file.

    Parameters
    ----------
    file_path: string, optional
      Path of the file containing the instantaneous drag coefficients;
      default: None.
    Re: float, optional
      Reynolds number;
      default: None.
    """
    if not (file_path or Re):
      print('[error] please provide path of file of Reynolds number')
      return
    print('[info] reading drag coefficients ...'),
    if not file_path:
      file_name = ('koumoutsakos_leonard_1995_'
                   'cylinder_dragCoefficientRe{}.dat'.format(Re))
      file_path = os.path.join(os.environ['SNAKE'],
                               'resources',
                               'results',
                               file_name)
    with open(file_path, 'r') as infile:
      times, drag = numpy.loadtxt(infile,
                                  dtype=float,
                                  comments='#',
                                  unpack=True)
    self.cd = Force(0.5 * times, drag)
    print('done')
