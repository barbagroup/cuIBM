"""
Implementation of the class `Simulation`.
"""

import os
import sys

import numpy
from matplotlib import pyplot
import pandas


class Simulation(object):
  """
  Simulation manager.
  """

  def __init__(self,
               software,
               description=None,
               directory=os.getcwd(),
               **kwargs):
    """
    Registers a simulation.

    Parameters
    ----------
    software: string
      Name of the software used for the simulation;
      choices: 'cuibm', 'petibm', 'ibamr', 'openfoam'.
    description: string, optional
      Description of the simulation;
      default: None.
    directory: string, optional
      Directory of the simulation;
      default: <current working directory>.
    **kwargs: dictionary
      Other attributes to create.
    """
    if description:
      description = description.replace('_', ' ')
    self.description = description
    self.directory = directory
    self.software = software.lower()
    # set extra arguments
    for key, value in kwargs.items():
      setattr(self, key, value)
    self._print_registration()
    self.fields = {}
    self.forces = []
    self._derive_class()

  def _print_registration(self):
    """
    Prints global info of the simulation.
    """
    print('\n[info] registering simulation ...')
    print('\t- description: {}'.format(self.description))
    print('\t- software: {}'.format(self.software))
    print('\t- directory: {}\n'.format(self.directory))

  def _derive_class(self):
    """
    Finds the appropriate child class based on the software used.
    """
    if self.software == 'cuibm':
      from .cuibm.simulation import CuIBMSimulation
      self.__class__ = CuIBMSimulation
    elif self.software == 'petibm':
      from .petibm.simulation import PetIBMSimulation
      self.__class__ = PetIBMSimulation
    elif self.software == 'openfoam':
      from .openfoam.simulation import OpenFOAMSimulation
      self.__class__ = OpenFOAMSimulation
    elif self.software == 'ibamr':
      from .ibamr.simulation import IBAMRSimulation
      self.__class__ = IBAMRSimulation
    else:
      print('[error] software indicated: {}'.format(self.software))
      print('[error] should be one of the followings: '
            'cuibm, petibm, openfoam, or ibamr')
      sys.exit(0)

  def get_mean_forces(self,
                      limits=(0.0, float('inf')),
                      last_period=False,
                      order=5):
    """
    Computes the time-averaged forces (or force-coefficients).

    Parameters
    ----------
    limits: 2-list of floats, optional
      Time-limits used to compute the mean;
      default: (0.0, +inf).
    last_period: boolean, optional
      Set 'True' to compute the mean over the last period of the signal;
      default: False.
    order: integer, optional
      Number of neighboring points used to define an extremum;
      default: 5.
    """
    for index, force in enumerate(self.forces):
      self.forces[index].get_mean(limits=limits,
                                  last_period=last_period,
                                  order=order)

  def get_strouhal(self,
                   L=1.0, U=1.0,
                   limits=(0.0, float('inf')),
                   order=5,
                   index=1):
    """
    Computes the Strouhal number based on the frequency of the force signal.

    Parameters
    ----------
    L: float, optional
      Characteristics length of the body;
      default: 1.0.
    U: float, optional
      Characteristics velocity of the body;
      default: 1.0.
    limits: 2-tuple of floats, optional
      Time-limits used as reference to compute the Strouhal number;
      default: (0.0, inf).
    order: integer, optional
      Number of neighbors used on each side to define an extremum;
      default: 5.
    index: integer, optional
      Index of the list of forces to use to compute the Strouhal number;
      default: 1 (most of the time, 1 corresponds to the lift force).
    """
    return self.forces[index].get_strouhal(L=L, U=U,
                                           limits=limits, order=order)

  def plot_forces(self,
                  indices=None, labels=None,
                  display_coefficients=False, coefficient=1.0,
                  limits=(0.0, float('inf'), 0.0, float('inf')),
                  style=None,
                  save_directory=None, save_name=None, fmt='png', dpi=100,
                  display_extrema=False, order=5,
                  display_guides=False, fill_between=False,
                  other_simulations=[], other_coefficients=[],
                  show=False):
    """
    Displays the forces into a figure.

    Parameters
    ----------
    indices: list of integers, optional
      List of the index of each force to display;
      default: None (all forces).
    labels: list of strings, optional
      Labels for each force to display;
      default: None (default labels are set).
    display_coefficients: boolean, optional
      Set 'True' if plotting force coefficients instead of forces;
      default: False.
    coefficient: float, optional
      scale coefficient to convert a force in a force coefficient;
      default: 1.0.
    limits: 4-tuple of floats, optional
      Limits of the axes (xmin, xmax, ymin, ymax);
      default: (0.0, +inf, 0.0, +inf).
    style: string, optional
      Path of the style-sheet to use for the figure;
      default: None.
    save_directory: string, optional
      Directory where to save the figure;
      default: None.
    save_name: string, optional
      Name of the .PNG file to save;
      default: None (does not save).
    fmt: string, optional
      Format of the file to save;
      default: 'png'.
    dpi: integer, optional
      Dots per inch (resolution);
      default: 100.
    display_extrema: boolean, optional
      Set 'True' to emphasize the extrema of the curves;
      default: False.
    order: integer, optional
      Number of neighbors used on each side to define an extreme;
      default: 5.
    display_guides: boolean, optional
      Set 'True' to display guides to judge steady regime;
      default: False.
    fill_between: boolean, optional
      Set 'True' to fill between lines defined by the extrema;
      default: False.
    other_simulations: list of Simulation objects or single Simulation object,
    optional
      List of other simulations to add to plot;
      default: [].
    other_coefficients: list of floats, optional
      Scale coefficients for each other simulation;
      default: [].
    show: boolean, optional
      Set 'True' to display the figure;
      default: False.
    """
    if not (save_name or show):
      return
    # convert other_simulations in list if single Simulation provided
    try:
      assert isinstance(other_simulations, (list, tuple))
    except:
      other_simulations = [other_simulations]
    # convert other_coefficients in list if single float provided
    try:
      assert isinstance(other_coefficients, (list, tuple))
    except:
      other_coefficients = [other_coefficients]
    print('\n[info] plotting forces ...')
    if style:
      try:
        pyplot.style.use(style)
      except:
        try:
          pyplot.style.use(os.path.join(os.environ['SNAKE'],
                                        'snake',
                                        'styles',
                                        style + '.mplstyle'))
        except:
          pass
    fig, ax = pyplot.subplots(figsize=(8, 6))
    color_cycle = ax._get_lines.prop_cycler
    ax.grid(True, zorder=0)
    ax.set_xlabel('time-unit')
    ax.set_ylabel('force coefficients' if display_coefficients else 'forces')
    if not labels:
      labels = (['$C_d$', '$C_l$', '$C_m$'] if display_coefficients
                else ['$F_x$', '$F_y$', '$F_z$'])
    if not indices:
      indices = numpy.arange(0, len(self.forces) + 1, 1)
    for index, force in enumerate(self.forces):
      if index not in indices:
        continue
      color = next(color_cycle)['color']
      line, = ax.plot(force.times, coefficient * force.values,
                      label=' - '.join(filter(None,
                                              [self.description,
                                               labels[index]])),
                      color=color, linestyle='-', zorder=9)
      if display_extrema:
        minima, maxima = force.get_extrema(order=order)
        ax.scatter(force.times[minima], coefficient * force.values[minima],
                   c=color, marker='o', zorder=10)
        ax.scatter(force.times[maxima], coefficient * force.values[maxima],
                   c=color, marker='o', zorder=10)
        if fill_between:
          line.remove()
          ax.plot(force.times[minima], coefficient * force.values[minima],
                  color='white', linestyle='-', zorder=9)
          ax.plot(force.times[maxima], coefficient * force.values[maxima],
                  color='white', linestyle='-', zorder=9)
          times = numpy.concatenate((force.times[minima],
                                     force.times[maxima][::-1]))
          vals = coefficient * numpy.concatenate((force.values[minima],
                                                  force.values[maxima][::-1]))
          ax.fill(times, vals,
                  label=' - '.join(filter(None,
                                          [self.description,
                                           labels[index]])),
                  facecolor=color, alpha=0.8, zorder=8)
        if display_guides:
          ax.axhline(coefficient * force.values[minima[-1]],
                     color=color, linestyle=':', zorder=10)
          ax.axhline(coefficient * force.values[maxima[-1]],
                     color=color, linestyle=':', zorder=10)

    for i, simulation in enumerate(other_simulations):
      for index, force in enumerate(simulation.forces):
        if index not in indices:
          continue
        color = next(color_cycle)['color']
        line, = ax.plot(force.times, other_coefficients[i] * force.values,
                        label=' - '.join(filter(None,
                                                [simulation.description,
                                                 labels[index]])),
                        color=color, linestyle='--', zorder=9)
        if fill_between:
          line.remove()
          minima, maxima = force.get_extrema(order=order)
          ax.scatter(force.times[minima],
                     other_coefficients[i] * force.values[minima],
                     c=color, marker='o', zorder=10)
          ax.scatter(force.times[maxima],
                     other_coefficients[i] * force.values[maxima],
                     c=color, marker='o', zorder=10)
          ax.plot(force.times[minima],
                  other_coefficients[i] * force.values[minima],
                  color='white', linestyle='-', zorder=9)
          ax.plot(force.times[maxima],
                  other_coefficients[i] * force.values[maxima],
                  color='white', linestyle='-', zorder=9)
          times = numpy.concatenate((force.times[minima],
                                     force.times[maxima][::-1]))
          values = (other_coefficients[index]
                    * numpy.concatenate((force.values[minima],
                                         force.values[maxima][::-1])))
          ax.fill(times, values,
                  label=' - '.join(filter(None,
                                          [simulation.description,
                                           labels[index]])),
                  facecolor=color, alpha=0.5, zorder=7)
    legend = ax.legend()
    legend.set_zorder(20)  # put legend on top
    ax.axis(limits)
    if save_name:
      if not save_directory:
        save_directory = os.path.join(self.directory, 'images')
      if not os.path.isdir(save_directory):
        os.makedirs(save_directory)
      print('[info] saving figure {}.{} in directory {} ...'
            ''.format(save_name, fmt, save_directory))
      pyplot.savefig(os.path.join(save_directory, '.'.join([save_name, fmt])),
                     dpi=dpi,
                     bbox_inches='tight',
                     format=fmt)
    if show:
      print('[info] displaying figure ...')
      pyplot.show()
    pyplot.close()

  def create_dataframe_forces(self,
                              indices=[],
                              display_coefficients=False,
                              coefficient=1.0,
                              display_strouhal=False,
                              labels=[]):
    """
    Creates a data frame with Pandas to display time-averaged forces
    (or force coefficients).

    Parameters
    ----------
    indices: list of integers, optional
      List of the index of each force to display in the dataframe;
      default: None (all forces);
    display_coefficients: boolean, optional
      Set 'True' if force coefficients are to be displayed;
      default: False.
    coefficient: float, optional
      Scale factor to convert a force into a force-coefficient;
      default: 1.0.
    display_strouhal: boolean, optional
      Set 'True' to display the mean Strouhal number in the dataframe;
      default: False.

    Returns
    -------
    dataframe: Pandas dataframe
      The dataframe of the simulation.
    """
    print('[info] instantaneous signals are averaged between '
          '{} and {} time-units.'.format(self.forces[0].mean['start'],
                                         self.forces[0].mean['end']))
    descriptions = ['<no description>' if not self.description
                    else self.description]
    # set default indices
    if not any(indices):
      indices = numpy.arange(0, len(self.forces), 1)
    # set default labels
    if not any(labels):
      if display_coefficients:
        labels = numpy.array(['<Cd>', '<Cl>', '<Cm>'])
      else:
        labels = numpy.array(['<Fx>', '<Fy>', '<Fz>'])
    # grab mean values of interest
    values = []
    for index, force in enumerate(self.forces):
      if index not in indices:
        continue
      values.append('{0:.4f}'.format(coefficient * force.mean['value']))
    dataframe = pandas.DataFrame([values],
                                 index=descriptions,
                                 columns=labels[indices])
    if display_strouhal:
      # assuming Strouhal number based on lift curve signal
      strouhal = self.forces[1].strouhal
      print('[info] Strouhal number is averaged '
            'over the {} periods of the lift curve '
            'between {} and {} time-units.'
            ''.format(strouhal['n-periods'],
                      strouhal['time-limits'][0], strouhal['time-limits'][1]))
      dataframe['<St>'] = '{0:.4f}'.format(strouhal['mean'])
    return dataframe
