"""
Implementation of the class `Field`.
"""

import os

import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


class Field(object):
  """
  Contains information about a 2D field (pressure for example) on a structured
  Cartesian grid.
  """

  def __init__(self, x=None, y=None, values=None, time_step=None, label=None):
    """
    Initializes the field by its grid and its values.

    Parameters
    ----------
    x: numpy 1D array of floats, optional
      Stations along a gridline in the x-direction;
      default: None.
    y: numpy 1D array of floats, optional
      Stations along a gridline in the y-direction;
      default: None.
    values: numpy 2D array of floats, optional
      Discrete field values;
      default: None.
    time_step: integer, optional
      Time-step;
      default: None.
    label: string, optional
      Description of the field;
      default: None.
    """
    self.label = None
    self.x, self.y = None, None
    self.values = None
    self.time_step = None
    if numpy.any(x) and numpy.any(y) and values.shape == (y.size, x.size):
      self.set(x, y, values, time_step=time_step, label=label)

  def set(self, x, y, values, time_step=None, label=None):
    """
    Sets the stations along a gridline in each direction, the field values,
    the time-step, and the label.

    Parameters
    ----------
    x: numpy 1D array of floats
      Stations along a gridline in the x-direction.
    y: numpy 1D array of floats
      Stations along a gridline in the y-direction.
    values: numpy 2D array of floats
      Discrete field values.
    time_step: integer, optional
      Time-step;
      default: None.
    label: string, optional
      Description of the field;
      default: None.
    """
    assert values.shape == (y.size, x.size)
    self.x, self.y = x, y
    self.values = values
    self.time_step = time_step
    self.label = label

  def subtract(self, other, label=None, atol=1.0E-12):
    """
    Subtracts a given field to the current one (returns 'self' - 'other').

    Note: the other field must be defined on the same grid.

    Parameters
    ----------
    other: Field object
      The field that is subtracted.
    label: string, optional
      Label of the Field object to create;
      default: None (will be '<current label>-subtracted').
    atol: float, optional
      Absolute-tolerance to define if two grid nodes have the same location;
      default: 1.0E-12.

    Returns
    -------
    subtracted_field: Field object
      The subtracted field.
    """
    # check the two solutions share the same grid
    assert numpy.allclose(self.x, other.x, atol=atol)
    assert numpy.allclose(self.y, other.y, atol=atol)
    assert self.values.shape == other.values.shape
    if not label:
      label = self.label + '-subtracted'
    return Field(label=label,
                 time_step=self.time_step,
                 x=self.x, y=self.y,
                 values=self.values - other.values)

  def restrict(self, x, y, label=None, atol=1.0E-12):
    """
    Restricts the field solution onto a coarser grid.
    Note: all nodes on the coarser grid are present in the actual grid.

    Parameters
    ----------
    x: numpy 1D array of floats
      Stations along a gridline in the x-direction.
    y: numpy 1D array of floats
      Stations along a gridline in the y-direction.
    label: string, optional
      Label of the restricted field;
      default: None (will be '<current label>-restricted').
    atol: float, optional
      Absolute tolerance used to define shared nodes between two grids;
      default: 1.0E-06.

    Returns
    -------
    restricted_field: Field object
      Field restricted onto the coarser grid.
    """
    def intersection(a, b, atol=atol):
      return numpy.any(numpy.abs(a - b[:, numpy.newaxis]) <= atol, axis=0)
    mask_x = intersection(self.x, x, atol=atol)
    mask_y = intersection(self.y, y, atol=atol)
    if not label:
      label = self.label + '-restricted'
    return Field(x=self.x[mask_x],
                 y=self.y[mask_y],
                 values=numpy.array([self.values[j][mask_x]
                                     for j in range(self.y.size)
                                     if mask_y[j]]),
                 time_step=self.time_step,
                 label=label)

  def get_difference(self, other, x, y, norm='L2'):
    """
    Returns the difference between two fields in a given norm and on a given
    grid.
    The grid is defined by its stations along a gridline in the x-direction and
    its stations along a gridline in the y-direction.

    Parameters
    ----------
    other: Field object
      The other field used for comparison.
    x: numpy 1D array of floats
      Stations along a gridline in the x-direction.
    y: numpy 1D array of floats
      Stations along a gridline in the y-direction.
    norms: string, optional
      Norm to use;
      choices: 'L2', 'Linf';
      default: 'L2'.

    Returns
    -------
    norm: float
      The difference in the given norm.
    """
    norms = {'L2': None, 'Linf': numpy.inf}
    field = self.restrict(x, y)
    other = other.restrict(x, y)
    subtracted = field.subtract(other)
    return numpy.linalg.norm(subtracted.values, ord=norms[norm])

  def get_gridline_values(self, x=None, y=None):
    """
    Returns the field values along either a vertical or an horizontal
    gridline.

    The vertical gridline is defined by its x-position.
    The horizontal gridline is defined by its y-position.

    Parameters
    ----------
    x: float, optional
      x-position of the vertical gridline;
      default: None.
    y: float, optional
      y-position of the horizontal gridline;
      default: None.

    Returns
    -------
    stations, values: two numpy 1D arrays of floats
      Stations and values along the gridline.
    """
    if (x and y) or not (x or y):
      print('[error] use either x or y keyword arguments '
            'to define the gridline position')
      return
    elif x:
      return self.get_vertical_gridline_values(x)
    elif y:
      return self.get_horizontal_gridline_values(y)

  def get_vertical_gridline_values(self, x):
    """
    Returns field values along a vertical gridline defined by its x-position.

    If the x-position of the gridline does not match any gridline of the
    Cartesian grid, we interpolate the values.

    Parameters
    ----------
    x: float
      x-position of the vertical gridline.

    Returns
    -------
    y, u: two numpy 1D arrays of floats
      Stations and values along the vertical gridline.
    """
    indices = numpy.where(numpy.abs(self.x - x) <= 1.0E-06)[0]
    # if no station matches the given value, we interpolate
    if indices.size == 0:
      i = numpy.where(self.x > x)[0][0]
      return (self.y, (abs(self.x[i] - x) * self.values[:, i - 1]
                       + abs(self.x[i - 1] - x) * self.values[:, i])
              / abs(self.x[i] - self.x[i - 1]))
    else:
      i = indices[0]
      return self.y, self.values[:, i]

  def get_horizontal_gridline_values(self, y):
    """
    Returns field values along an horizontal gridline defined by its
    y-position.

    If the y-position of the gridline does not match any gridline of the
    Cartesian grid, we interpolate the values.

    Parameters
    ----------
    y: float
      y-position of the horizontal gridline.

    Returns
    -------
    x, u: two numpy 1D arrays of floats
      Stations and values along the horizontal gridline.
    """
    indices = numpy.where(numpy.abs(self.y - y) <= 1.0E-06)[0]
    # if no station matches the given value, we interpolate
    if indices.size == 0:
      j = numpy.where(self.y > y)[0][0]
      return (self.x, (abs(self.y[j] - y) * self.values[j - 1, :]
                       + abs(self.y[j - 1] - y) * self.values[j, :])
              / abs(self.y[j] - self.y[j - 1]))
    else:
      j = indices[0]
      return self.x, self.values[j, :]

  def plot_vertical_gridline_values(self, x,
                                    boundaries=(None, None),
                                    plot_settings={},
                                    plot_limits=(None, None, None, None),
                                    save_directory=None,
                                    show=False,
                                    other_data=None,
                                    other_plot_settings={},
                                    style=None):
    """
    Plots the field values along a group of vertical gridlines.

    Parameters
    ----------
    x: list of floats
      Group of vertical gridlines defined by their x-position.
    boundaries: 2-tuple of floats, optional
      Gridline limits to consider;
      default: (None, None) (the entire gridline)
    plot_settings: dictionary of (string, object) items, optional
      Contains optional arguments to call pyplot.plot function for the gridline
      data;
      default: empty dictionary.
    plot_limits: 4-tuple of floats, optional
      Limits of the plot (x-start, x-end, y-start, y-end);
      default: (None, None, None, None)
    save_directory: string, optional
      Directory where to save the figure;
      default: None (does not save).
    show: boolean, optional
      Set 'True' if you want to display the figure;
      default: False.
    other_data: 2-tuple of 1d arrays of floats, optional
      Other data to add to the figure (1st array contains the y-stations,
      2nd array contains the values at the stations);
      default: None.
    other_plot_settings: dictionary of (string, object) items, optional
      Contains optional arguments to call pyplot.plot function
      for the other data;
      default: empty dictionary.
    style: string, optional
      Path of a Matplotlib style-sheet;
      default: None.
    """
    print('[info] plotting field values along vertical gridline(s) ...'),
    if style:
      pyplot.style.use(style)
    fig, ax = pyplot.subplots(figsize=(6, 6))
    ax.grid(True, zorder=0)
    ax.set_xlabel('y-coordinate', fontsize=16)
    ax.set_ylabel('{} along vertical gridline'.format(self.label), fontsize=16)
    if not isinstance(x, (list, tuple)):
      x = [x]
    for x_target in x:
      y, u = self.get_vertical_gridline_values(x_target)
      if all(boundaries):
        mask = numpy.where(numpy.logical_and(y >= boundaries[0],
                                             y <= boundaries[1]))[0]
        y, u = y[mask], u[mask]
      ax.plot(y, u, **plot_settings)
    if other_data:
      y, u = other_data
      if all(boundaries):
        mask = numpy.where(numpy.logical_and(y >= boundaries[0],
                                             y <= boundaries[1]))[0]
        y, u = y[mask], u[mask]
      ax.plot(y, u, **other_plot_settings)
    ax.axis(plot_limits)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax.legend(prop={'size': 16})
    if save_directory and os.path.isdir(save_directory):
      file_name = '{}VerticalGridline{:0>7}.png'.format(self.label,
                                                        self.time_step)
      pyplot.savefig(os.path.join(save_directory, file_name))
    if show:
      pyplot.show()
    print('done')

  def plot_horizontal_gridline_values(self, y,
                                      boundaries=(None, None),
                                      plot_settings={},
                                      plot_limits=(None, None, None, None),
                                      save_directory=None,
                                      show=False,
                                      other_data=None,
                                      other_plot_settings={},
                                      style=None):
    """
    Plots the field values along a group of horizontal gridlines.

    Parameters
    ----------
    y: list of floats
      Group of horizontal gridlines defined by their y-position.
    boundaries: 2-tuple of floats, optional
      Gridline limits to consider;
      default: (None, None) (the entire gridline)
    plot_settings: dictionary of (string, object) items, optional
      Contains optional arguments to call pyplot.plot function for the gridline
      data;
      default: empty dictionary.
    plot_limits: 4-tuple of floats, optional
      Limits of the plot (x-start, x-end, y-start, y-end);
      default: (None, None, None, None)
    save_directory: string, optional
      Directory where to save the figure;
      default: None (does not save).
    show: boolean, optional
      Set 'True' if you want to display the figure;
      default: False.
    other_data: 2-tuple of 1d arrays of floats, optional
      Other data to add to the figure (1st array contains the y-stations,
      2nd array contains the values at the stations);
      default: None.
    other_plot_settings: dictionary of (string, object) items, optional
      Contains optional arguments to call pyplot.plot function for the other
      data;
      default: empty dictionary.
    style: string, optional
      Path of a Matplotlib style-sheet;
      default: None.
    """
    print('[info] plotting field values along horizontal gridline(s) ...'),
    if style:
      pyplot.style.use(style)
    fig, ax = pyplot.subplots(figsize=(6, 6))
    ax.grid(True, zorder=0)
    ax.set_xlabel('x-coordinate', fontsize=16)
    ax.set_ylabel('{} along horizontal gridline'.format(self.label),
                  fontsize=16)
    if not isinstance(y, (list, tuple)):
      y = [y]
    for y_target in y:
      x, u = self.get_horizontal_gridline_values(y_target)
      if all(boundaries):
        mask = numpy.where(numpy.logical_and(x >= boundaries[0],
                                             x <= boundaries[1]))[0]
        x, u = x[mask], u[mask]
      ax.plot(x, u, **plot_settings)
    if other_data:
      x, u = other_data
      if all(boundaries):
        mask = numpy.where(numpy.logical_and(x >= boundaries[0],
                                             x <= boundaries[1]))[0]
        x, u = x[mask], u[mask]
      ax.plot(x, u, **other_plot_settings)
    ax.axis(plot_limits)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax.legend(prop={'size': 16})
    if save_directory and os.path.isdir(save_directory):
      file_name = '{}HorizontalGridline{:0>7}.png'.format(self.label,
                                                          self.time_step)
      pyplot.savefig(os.path.join(save_directory, file_name))
    if show:
      pyplot.show()
    print('done')

  def plot_contour(self,
                   field_range=None,
                   filled_contour=True,
                   view=[float('-inf'), float('-inf'),
                         float('inf'), float('inf')],
                   bodies=[],
                   time_increment=None,
                   save_name=None,
                   save_directory=os.getcwd(),
                   fmt='png',
                   colorbar=True,
                   cmap=None,
                   colors=None,
                   width=8.0,
                   dpi=100):
    """
    Plots and saves the field.

    Parameters
    ----------
    field_range: 3-list of floats, optional
      Min, max and number of contours to plot;
      default: None.
    filled_contour: boolean, optional
      Set 'True' to create a filled contour;
      default: True.
    view: 4-list of floats, optional
      Bottom-left and top-right coordinates of the rectangular view to plot;
      default: the whole domain.
    bodies: list of Body objects or single Body object, optional
      The immersed bodies to add to the figure;
      default: [] (no immersed body).
    time_increment: float, optional
      Time-increment used to advance the simulation;
      default: None.
    save_name: string, optional
      Prefix used to create the images directory and to save the .png files;
      default: None.
    save_directory: string, optional
      Directory where to save the image;
      default: '<current directory>'.
    fmt: string, optional
      Format of the file to save;
      default: 'png'.
    colorbar: boolean, optional
      Set 'True' to display an horizontal colorbar at the bottom-left of the
      figure;
      default: True.
    cmap: string, optional
      The Matplotlib colormap to use;
      default: None.
    colors: string, optional
      The Matplotlib colors to use;
      default: None.
    width: float, optional
      Width of the figure (in inches);
      default: 8.
    dpi: integer, optional
      Dots per inch (resolution);
      default: 100
    """
    if abs(self.values.min() - self.values.max()) <= 1.0E-06:
      print('[warning] uniform field; plot contour skipped!')
      return
    # convert bodies in list if single body provided
    if not isinstance(bodies, (list, tuple)):
      bodies = [bodies]
    print('[time-step {}] plotting the {} contour ...'.format(self.time_step,
                                                              self.label))
    height = width * (view[3] - view[1]) / (view[2] - view[0])
    fig, ax = pyplot.subplots(figsize=(width, height), dpi=dpi)
    ax.tick_params(axis='x', labelbottom='off')
    ax.tick_params(axis='y', labelleft='off')
    # create filled contour
    if field_range:
      levels = numpy.linspace(*field_range)
      print('\tmin={}, max={}'.format(self.values.min(), self.values.max()))
      colorbar_ticks = numpy.linspace(field_range[0], field_range[1], 5)
      colorbar_format = '%.01f'
    else:
      levels = numpy.linspace(self.values.min(), self.values.max(), 101)
      print('\tmin={}, max={}, steps={}'.format(levels[0],
                                                levels[-1],
                                                levels.size))
      colorbar_ticks = numpy.linspace(self.values.min(), self.values.max(), 3)
      colorbar_format = '%.04f'
    color_map = {'pressure': cm.jet, 'vorticity': cm.RdBu_r,
                 'x-velocity': cm.RdBu_r, 'y-velocity': cm.RdBu_r}
    X, Y = numpy.meshgrid(self.x, self.y)
    contour_type = ax.contourf if filled_contour else ax.contour
    if not colors:
      if not cmap:
        cmap = (cm.RdBu_r if self.label not in color_map.keys()
                else color_map[self.label])
    cont = contour_type(X, Y, self.values,
                        levels=levels, extend='both',
                        colors=colors,
                        cmap=cmap)
    if colorbar:
      ains = inset_axes(pyplot.gca(), width='30%', height='2%', loc=3)
      cont_bar = fig.colorbar(cont,
                              cax=ains, orientation='horizontal',
                              ticks=colorbar_ticks, format=colorbar_format)
      cont_bar.ax.tick_params(labelsize=10)
      cont_bar.ax.xaxis.set_ticks_position('top')
    if time_increment:
      ax.text(0.05, 0.85,
              '{} time-units'.format(time_increment * self.time_step),
              transform=ax.transAxes, fontsize=10)
    # draw body
    for body in bodies:
      ax.plot(body.x, body.y,
              color='black', linewidth=1, linestyle='-')
    # set limits
    ax.set_xlim(view[::2])
    ax.set_ylim(view[1::2])
    ax.set_aspect('equal')
    # save image
    save_name = (self.label if not save_name else save_name)
    file_path = os.path.join(save_directory,
                             '{}{:0>7}.{}'.format(save_name,
                                                  self.time_step,
                                                  fmt))
    pyplot.savefig(file_path,
                   dpi=dpi, bbox_inches='tight', pad_inches=0,
                   format=fmt)
    pyplot.close()
