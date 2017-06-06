"""
Implementation of the class `Force`.
A `Force` object contains information about an instantaneous force
(e.g., fx or fy).
"""

import numpy
from scipy import signal


class Force(object):
  """
  Contains info about an instantaneous force.
  """

  def __init__(self, times=None, values=None, label=None):
    """
    Initializes the force.

    Parameters
    ----------
    times: numpy 1D array of floats, optional
      Discrete time values;
      default: None.
    values: numpy 1D array of floats, optional
      Instantaneous values of the force;
      default: None.
    label: string, optional
      Description of the force;
      default: None.
    """
    self.label = None
    self.times = None
    self.values = None
    self.mean = None
    self.deviations = None
    self.strouhal = None
    if numpy.any(times) and numpy.any(values):
      self.set(times, values, label=label)

  def set(self, times, values, label=None):
    """
    Sets the discrete time and values of the instantaneous force, as well as
    its label.

    Parameters
    ----------
    times: numpy 1D array of floats
      Discrete time values.
    values: numpy 1D array of floats
      Instantaneous values of the force.
    label: string, optional
      Description of the force;
      default: None.
    """
    assert times.size == values.size
    self.times, self.values = times, values
    self.label = label

  def get_mean(self, limits=(0.0, float('inf')), last_period=False, order=5):
    """
    Computes the mean force.

    Parameters
    ----------
    limits: 2-tuple of floats, optional
      Time-limits to compute the mean value;
      default: (0.0, float('inf')).
    last_period: boolean, optional
      If 'True': computes the mean value over the last period;
      default: False.
    order: integer, optional
      If `last_period=True`: number of neighbors used to define an extremum;
      default: 5.

    Returns
    -------
    mean: dictionary of (string, float) items
      Keys: - 'value': the mean value;
            - 'start': the actual starting time value used;
            - 'end': the actual end time value used.
    """
    if last_period:
      minima, maxima = self.get_extrema(order=order)
      mask = (minima if minima[-1] > maxima[-1] else maxima)[-2:]
    else:
      mask = numpy.where(numpy.logical_and(self.times >= limits[0],
                                           self.times <= limits[1]))[0]
    self.mean = {'value': numpy.mean(self.values[mask]),
                 'start': self.times[mask[0]],
                 'end': self.times[mask[-1]]}
    return self.mean

  def get_deviations(self, limits=(0.0, float('inf')), order=5):
    """
    Computes the deviations around the mean value.

    Parameters
    ----------
    limits: 2-tuple of floats, optional
      Time-limits to compute the mean value;
      default: [0.0, float('inf')].
    order: integer, optional
      Number of neighboring points used to define an extreme;
      default: 5.

    Returns
    -------
    deviations: dictionary of 2  (string, 1D array of floats) items
      Absolute deviations of the minima and maxima with respect to the mean
      value.
    """
    minima, maxima = self.get_extrema(limits=limits, order=order)
    mean = self.get_mean(limits=limits)['value']
    self.deviations = {'min': numpy.absolute(self.values[minima] - mean),
                       'max': numpy.absolute(self.values[maxima] - mean)}
    return self.deviations

  def get_extrema(self, limits=(0.0, float('inf')), order=5):
    """
    Computes masks (i.e. arrays of indices) of the extrema of the force.

    Parameters
    ----------
    order: integer, optional
      Number of neighboring points used to define an extreme;
      default: 5.

    Returns
    -------
    minima: 1D array of integers
      Index of all minima.
    maxima: 1D array of integers
      Index of all maxima.
    """
    minima = signal.argrelextrema(self.values, numpy.less_equal,
                                  order=order)[0][:-1]
    maxima = signal.argrelextrema(self.values, numpy.greater_equal,
                                  order=order)[0][:-1]
    mask = numpy.where(numpy.logical_and(self.times >= limits[0],
                                         self.times <= limits[1]))[0]
    minima = numpy.intersect1d(minima, mask, assume_unique=True)
    maxima = numpy.intersect1d(maxima, mask, assume_unique=True)
    # remove indices that are too close
    minima = minima[numpy.append(True, minima[1:] - minima[:-1] > order)]
    maxima = maxima[numpy.append(True, maxima[1:] - maxima[:-1] > order)]
    return minima, maxima

  def get_strouhal(self, L=1.0, U=1.0, limits=(0.0, float('inf')), order=5):
    """
    Computes the Strouhal number based on the frequency of the signal.

    The frequency is computed using the minima of the signal.

    Parameters
    ----------
    L: float, optional
      Characteristics length of the body;
      default: 1.0.
    U: float, optional
      Characteristics velocity of the body;
      default: 1.0.
    n_periods: integer, optional
      Number of periods (starting from end) to average the Strouhal number;
      default: 1.
    limits: 2-tuple of floats, optional
      Time-limits used to compute the Strouhal number;
      default: (0.0, inf).
    order: integer, optional
      Number of neighbors used on each side to define an extremum;
      default: 5.

    Returns
    -------
    strouhal: dictionary
      Returns the mean Strouhal number, the number of periods used for
      averaging, the actual time-limits used, and the value over each period
      used.
    """
    minima, _ = self.get_extrema(limits=limits, order=order)
    strouhals = L / U / (self.times[minima[1:]] - self.times[minima[:-1]])
    self.strouhal = {'n-periods': minima.size - 1,
                     'time-limits': (self.times[minima[0]],
                                     self.times[minima[-1]]),
                     'values': strouhals,
                     'mean': strouhals.mean()}
    return self.strouhal
