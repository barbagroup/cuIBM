"""
Implementation of the class `IBAMRSimulation`.
"""

import os

import numpy

from ..simulation import Simulation
from ..force import Force


class IBAMRSimulation(Simulation):
  """
  Contains info about a IBAMR simulation.
  Inherits from class Simulation.
  """

  def __init__(self, description=None, directory=os.getcwd(), **kwargs):
    """
    Initializes by calling the parent constructor.

    Parameters
    ----------
    description: string, optional
      Description of the simulation;
      default: None.
    directory: string, optional
      Directory of the simulation;
      default: <current working directory>.
    """
    super(IBAMRSimulation, self).__init__(software='ibamr',
                                          description=description,
                                          directory=directory,
                                          **kwargs)

  def read_forces(self, file_path=None, labels=None):
    """
    Reads forces from files.

    Parameters
    ----------
    file_path: string, optional
      Path of the file containing the forces;
      default: None.
    labels: list of strings, optional
      Label of each force to read;
      default: None.
    """
    if not file_path:
      file_path = os.path.join(self.directory,
                               'dataIB',
                               'ib_Drag_force_struct_no_0')
    print('[info] reading forces from {} ...'.format(file_path)),
    with open(file_path, 'r') as infile:
      times, force_x, force_y = numpy.loadtxt(infile,
                                              dtype=float,
                                              usecols=(0, 4, 5),
                                              unpack=True)
    self.forces = []
    self.forces.append(Force(times, force_x, '$F_x$'))
    self.forces.append(Force(times, force_y, '$F_y$'))
    print('done')

  def write_visit_summary_files(self, time_steps):
    """
    Writes summary files for Visit with list of sub-directories to look into.

    Parameters
    ----------
    time_steps: 3-tuple of integers
      Staring and and ending time_steps followed by the time-step increment.
    """
    print('[info] writing summary files for VisIt ...'),
    time_steps = numpy.arange(time_steps[0], time_steps[1] + 1, time_steps[2])
    # list of SAMRAI files to VisIt
    dumps_visit_list = [os.path.join('visit_dump.{0:05}'.format(time_step),
                                     'summary.samrai')
                        for time_step in time_steps]
    dumps_visit = numpy.array(dumps_visit_list)
    # list of SILO file to VisIt
    lag_visit_list = [os.path.join('lag_data.cycle_{:06d}'.format(time_step),
                                   'lag_data.cycle_{:06d}.summary.silo'
                                   ''.format(time_step))
                      for time_step in time_steps]
    lag_data_visit = numpy.array(lag_visit_list)
    # write files
    with open(os.path.join(self.directory, 'dumps.visit'), 'w') as outfile:
      numpy.savetxt(outfile, dumps_visit, fmt='%s')
    with open(os.path.join(self.directory, 'lag_data.visit'), 'w') as outfile:
      numpy.savetxt(outfile, lag_data_visit, fmt='%s')
    print('done')

  def plot_field_contours_visit(self, field_name,
                                field_range,
                                body=None,
                                solution_folder='numericalSolution',
                                states=(0, 20000, 1),
                                view=(-2.0, -2.0, 2.0, 2.0),
                                width=800):
    """
    Plots the contour of a given field using VisIt.

    Parameters
    ----------
    field_name: string
      Name of field to plot;
      choices: vorticity, pressure, velocity-magnitude, x-velocity, y-velocity.
    field_range: 2-tuple of floats
      Range of the field to plot (min, max).
    body: string, optional
      Name of the immersed body;
      default: None.
    solution_folder: string, optional
      Relative path of the folder containing the numerical solution;
      default: 'numericalSolution'.
    states: 3-tuple of integers, optional
      Limits of index of the states to plot followed by the increment;
      default: (0, 20000, 1).
    view: 4-tuple of floats, optional
      Bottom-left and top-right coordinates of the view to display;
      default: (-2.0, -2.0, 2.0, 2.0).
    width: integer, optional
      Width (in pixels) of the figure;
      default: 800.
    """
    args = {}
    args['--directory'] = self.directory
    args['--field'] = field_name
    args['--range'] = '{} {}'.format(*field_range)
    if body:
      args['--body'] = body
    args['--solution-folder'] = solution_folder
    args['--states'] = '{} {} {}'.format(*states)
    args['--view'] = '{} {} {} {}'.format(*view)
    args['--width'] = str(width)
    script = os.path.join(os.environ['SNAKE'],
                          'snake',
                          'ibamr',
                          'plotField2dVisIt.py')
    arguments = ' '.join([key + ' ' + value
                          for key, value in args.items()])
    os.system('visit -nowin -cli -s {} {}'.format(script, arguments))

  def compute_mean_number_cells_visit(self,
                                      solution_folder='numericalSolution',
                                      states=(0, 20000, 1),
                                      time_limits=(0.0, float('inf'))):
    """
    Computes the number of a cells, on average, in the AMR mesh at requested
    states.

    Parameters
    ----------
    solution_folder: string, optional
      Relative path of the folder containing the numerical solution;
      default: 'numericalSolution'.
    states: 3-tuple of integers, optional
      Starting and ending states followed by the increment;
      default: (0, 20000, 1).
    time_limits: 2-tuple of floats, optional
      Time-limits within which the mean value is calculated;
      default: (0.0, inf).
    """
    args = {}
    args['--directory'] = self.directory
    args['--solution-folder'] = solution_folder
    args['--states'] = '{} {} {}'.format(*states)
    args['--time-limits'] = '{} {}'.format(*time_limits)
    script = os.path.join(os.environ['SNAKE'],
                          'snake',
                          'ibamr',
                          'getNumberCellsVisIt.py')
    arguments = ' '.join([key + ' ' + value
                          for key, value in args.items()])
    os.system('visit -nowin -cli -s {} {}'.format(script, arguments))
