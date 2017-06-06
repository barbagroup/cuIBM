"""
Implementation of the class `OpenFOAMSimulation`.
"""

import os

import numpy
from scipy import signal
from matplotlib import pyplot

from ..simulation import Simulation
from ..force import Force


class OpenFOAMSimulation(Simulation):
  """
  Contains info about a OpenFOAM simulation.
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
    super(OpenFOAMSimulation, self).__init__(software='openfoam',
                                             description=description,
                                             directory=directory,
                                             **kwargs)

  def read_forces(self,
                  display_coefficients=False,
                  labels=None,
                  forces_folder=os.path.join('postProcessing',
                                             'forces'),
                  force_coefficients_folder=os.path.join('postProcessing',
                                                         'forceCoeffs'),
                  usecols=(0, 2, 3)):
    """
    Reads forces from files.

    Parameters
    ----------
    display_coefficients: boolean, optional
      Set to 'True' if force coefficients are required;
      default: False (i.e. forces).
    labels: list of strings, optional
      Label of each force to read;
      default: None.
    forces_folder: string, optional
      Relative path from the simulation directory to the folder containing
      the forces;
      default: 'postProcessing/forces'.
    force_coefficients_folder: string, optional
      Relative path from the simulation directory to the folder containing
      the force coefficients;
      default: 'postProcessing/forceCoeffs'.
    usecols: tuple of integers, optional
      Index of columns to read from file, including the time-column index;
      default: (0, 2, 3).
    """
    if display_coefficients:
      info = {'directory': os.path.join(self.directory,
                                        force_coefficients_folder),
              'file-name': 'forceCoeffs.dat',
              'description': 'force-coefficients'}
      if not labels:
        labels = ['$C_d$', '$C_l$']
    else:
      info = {'directory': os.path.join(self.directory,
                                        forces_folder),
              'file-name': 'forces.dat',
              'description': 'forces'}
      if not labels:
        labels = ['$F_x$', '$F_y$']
    info['usecols'] = usecols
    info['labels'] = labels
    # backward compatibility from 2.2.2 to 2.0.1
    if not os.path.isdir(info['directory']):
      info['directory'] = '{}/forces'.format(self.directory)
      info['usecols'] = (0, 1, 2)
    # end of backward compatibility
    print('[info] reading {} in {} ...'.format(info['description'],
                                               info['directory']))
    subdirectories = sorted(os.listdir(info['directory']))
    times = numpy.empty(0)
    force_x, force_y = numpy.empty(0), numpy.empty(0)
    for subdirectory in subdirectories:
      forces_path = os.path.join(info['directory'],
                                 subdirectory,
                                 info['file-name'])
      with open(forces_path, 'r') as infile:
        t, fx, fy = numpy.loadtxt(infile,
                                  dtype=float,
                                  comments='#',
                                  usecols=info['usecols'],
                                  unpack=True)
      times = numpy.append(times, t)
      force_x, force_y = numpy.append(force_x, fx), numpy.append(force_y, fy)
    # set Force objects
    self.forces = []
    self.forces.append(Force(times, force_x, label=labels[0]))
    self.forces.append(Force(times, force_y, label=labels[1]))

  def read_maximum_cfl(self, file_path):
    """
    Reads the instantaneous maximum CFL number from a given log file.

    Parameters
    ----------
    file_path: string
      Path of the logging file containing the instantaneous maximum CFL number.

    Returns
    -------
    cfl: dictionary of (string, 1D array of floats) items
      Contains the discrete time and cfl values.
    """
    print('[info] reading CFL from {} ...'.format(file_path)),
    with open(file_path, 'r') as infile:
      times = numpy.array([float(line.split()[-1])
                           for line in infile if line.startswith('Time = ')])
    with open(file_path, 'r') as infile:
      cfl = numpy.array([float(line.split()[-1])
                         for line in infile
                         if line.startswith('Courant Number mean')])
    assert(times.shape == cfl.shape)
    self.cfl = {'times': times, 'values': cfl}
    print('done')
    return self.cfl

  def get_mean_maximum_cfl(self, limits=(0.0, float('inf'))):
    """
    Computes the mean CFL number.

    Parameters
    ----------
    limits: list of floats, optional
      Time-limits to compute the mean value;
      default: (0.0, float('inf')).

    Returns
    -------
    mean: dictionary of (string, float) items
      The mean value and the actual time-limits used to average the CFL.
    """
    print('[info] computing the mean CFL number ...')
    mask = numpy.where(numpy.logical_and(self.cfl['times'] >= limits[0],
                                         self.cfl['times'] <= limits[1]))[0]
    self.cfl['mean'] = {'start': self.cfl['times'][mask[0]],
                        'end': self.cfl['times'][mask[-1]],
                        'value': self.cfl['values'].mean()}
    print('[info] averaging the maximum CFL number '
          'between {} and {} time-units:'.format(self.cfl['mean']['start'],
                                                 self.cfl['mean']['end']))
    print('\t<max(CFL)> = {}'.format(self.cfl['mean']['value']))
    return self.cfl['mean']

  def plot_maximum_cfl(self,
                       display_extrema=False, order=5,
                       limits=(0.0, float('inf'), 0.0, float('inf')),
                       save_directory=None, save_name=None, fmt='png',
                       style='mesnardo',
                       show=False):
    """
    Plots the instantaneous maximum CFL number.

    Parameters
    ----------
    time: 1d array of floats
      Discrete time values.
    cfl: 1d array of floats
      Maximum CFL values.
    display_extrema: boolean, optional
      Set 'True' to emphasize the extrema of the curves;
      default: False.
    order: integer, optional
      Number of neighbors used on each side to define an extremum;
      default: 5.
    limits: list of floats, optional
      Limits of the axes [xmin, xmax, ymin, ymax];
      default: [0.0, +inf, 0.0, +inf].
    directory: string, optional
      Directory of the simulation;
      default: <current directory>.
    save_directory: string, optional
      Directory where to save the figure;
      default: None (will be '<simulation directory>/images').
    save_name: string, optional
      Name of the file to save;
      default: None (does not save).
    fmt: string, optional
      Format to save the figure;
      default: 'png'.
    style: string, optional
      Name of the .mplstyle file that contains to the style.
      The file should be located in the folder 'snake/styles';
      default: 'mesnardo'.
    show: boolean, optional
      Set 'True' to display the figure;
      default: False.
    """
    print('[info] plotting cfl ...')
    pyplot.style.use(os.path.join(os.environ['SNAKE'],
                                  'snake',
                                  'styles',
                                  style + '.mplstyle'))
    fig, ax = pyplot.subplots(figsize=(8, 6))
    color_cycle = ax._get_lines.prop_cycler
    color = next(color_cycle)['color']
    ax.grid(True, zorder=0)
    ax.set_xlabel('time', fontsize=18)
    ax.set_ylabel('maximum CFL', fontsize=18)
    ax.plot(self.cfl['times'], self.cfl['values'], color=color, zorder=10)
    if display_extrema:
      minima = signal.argrelextrema(self.cfl['values'], numpy.less_equal,
                                    order=order)[0][:-1]
      maxima = signal.argrelextrema(self.cfl['values'], numpy.greater_equal,
                                    order=order)[0][:-1]
      # remove indices that are too close
      minima = minima[numpy.append(True, minima[1:] - minima[:-1] > order)]
      maxima = maxima[numpy.append(True, maxima[1:] - maxima[:-1] > order)]
      ax.scatter(self.cfl['times'][minima], self.cfl['values'][minima],
                 c=color, marker='o', zorder=10)
      ax.scatter(self.cfl['times'][maxima], self.cfl['values'][maxima],
                 c=color, marker='o', zorder=10)
    ax.axis(limits)
    if save_name:
      if not save_directory:
        save_directory = os.path.join(self.directory, 'images')
      print('[info] saving figure in directory {} ...'.format(save_directory))
      if not os.path.isdir(save_directory):
        os.makedirs(save_directory)
      pyplot.savefig(os.path.join(save_directory, save_name + '.' + fmt),
                     bbox_inches='tight',
                     format=fmt)
    if show:
      print('[info] displaying figure ...')
      pyplot.show()
    pyplot.close()

  def create_matplotlib_colormap(self, colormap_name, file_path=None):
    """
    Writes the values of a Matplotlib colormap into a temporary file located in
    the current working directory.
    The list of Matplotlib colormaps is available
    [here](http://matplotlib.org/examples/color/colormaps_reference.html))

    Parameters
    ----------
    colormap_name: string
      Name of the Matplotlib colormap to write into a file.
    file_path: string, optional
      Path of the file to write;
      default: None (will be '<colormap_name>_tmp.dat').

    Returns
    -------
    file_path: string
      Path of the file created.
    """
    from matplotlib import cm
    file_path = os.path.join(os.getcwd(),
                             colormap_name + '_tmp.dat')
    print('[info] write colormap {} from Matplotlib into file {} ...'
          ''.format(colormap_name, file_path))
    with open(file_path, 'w') as outfile:
      colormap_object = getattr(cm, colormap_name)
      try:
        colors = colormap_object.colors
      except:
        colors = []
        for i in range(colormap_object.N):
          colors.append(colormap_object(i)[:-1])
      for color in colors:
        outfile.write('{}, {}, {}\n'.format(*color))
    return file_path

  def plot_field_contours_paraview(self, field_name,
                                   field_range=(-1.0, 1.0),
                                   view=(-2.0, -2.0, 2.0, 2.0),
                                   times=(0, 0, 0),
                                   width=800,
                                   colormap=None,
                                   display_scalar_bar=True,
                                   display_time_text=True,
                                   display_mesh=False):
    """
    Plots the contour of a given field using ParaView.

    Parameters
    ----------
    field_name: string
      Name of field to plot;
      choices: vorticity, pressure, x-velocity, y-velocity.
    field_range: 2-tuple of floats, optional
      Range of the field to plot (min, max);
      default: (-1.0, 1.0).
    view: 4-tuple of floats, optional
      Bottom-left and top-right coordinates of the view to display;
      default: (-2.0, -2.0, 2.0, 2.0).
    times: 3-tuple of floats, optional
      Time-limits followed by the time-increment to consider;
      default: (0, 0, 0).
    width: integer, optional
      Width (in pixels) of the figure;
      default: 800.
    colormap: string, optional
      Name of the Matplotlib colormap to use;
      default: None.
    display_scalar_bar: boolean, optional
      Displays the scalar bar;
      default: True.
    display_time_text: boolean, optional
      Displays the time-unit in the top-left corner;
      default: True.
    display_mesh: boolean, optional
      Displays the mesh (Surface with Edges);
      default: False
    """
    # create the command-line parameters
    arguments = []
    arguments.append('--directory ' + self.directory)
    arguments.append('--field ' + field_name)
    arguments.append('--range {} {}'.format(*field_range))
    arguments.append('--times {} {} {}'.format(*times))
    arguments.append('--view {} {} {} {}'.format(*view))
    arguments.append('--width {}'.format(width))
    if display_mesh:
      arguments.append('--mesh')
    if not display_scalar_bar:
      arguments.append('--no-scalar-bar')
    if not display_time_text:
      arguments.append('--no-time-text')
    if colormap:
      colormap_path = self.create_matplotlib_colormap(colormap_name=colormap)
      arguments.append('--colormap ' + colormap_path)
    # execute the Python script with pvbatch
    script = os.path.join(os.environ['SNAKE'],
                          'snake',
                          'openfoam',
                          'plotField2dParaView.py')
    os.system('pvbatch {} {}'.format(script, ' '.join(arguments)))
    if colormap:
      os.remove(colormap_path)

  def plot_mesh_paraview(self,
                         view=(-2.0, -2.0, 2.0, 2.0),
                         width=800):
    """
    Plots the mesh (black lines on white background) using ParaView.

    Parameters
    ----------.
    view: 4-tuple of floats, optional
      Bottom-left and top-right coordinates of the view to display;
      default: (-2.0, -2.0, 2.0, 2.0).
    width: integer, optional
      Width (in pixels) of the figure;
      default: 800.
    """
    arguments = []
    arguments.append('--directory ' + self.directory)
    arguments.append('--view {} {} {} {}'.format(*view))
    arguments.append('--width {}'.format(width))
    # execute the Python script with pvbatch
    script = os.path.join(os.environ['SNAKE'], 'snake', 'openfoam',
                          'plotMesh2dParaView.py')
    os.system('pvbatch {} {}'.format(script, ' '.join(arguments)))
