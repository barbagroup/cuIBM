"""
Implementation of the class `MovingVortices`.
"""

import os
import sys

import numpy

from ..field import Field


class MovingVortices(object):
  """
  Analytical plug-in for the moving vortices case.
  """

  def __init__(self, x, y, time):
    """
    Computes the velocity and pressure fields on a given grid.

    Parameters
    ----------
    x, y: 2D numpy arrays of floats
      Contains the stations along the gridline in each direction.
    time: float
      Time at which the analytical solution will be computed.
    """
    self.bottom_left, self.top_right = [x[0], y[0]], [x[-1], y[-1]]
    self.fields = {}
    x_u, y_u = x[1:-1], 0.5 * (y[:-1] + y[1:])
    self.fields['x-velocity'], _ = self.get_velocity(x_u, y_u, time)
    x_v, y_v = 0.5 * (x[:-1] + x[1:]), y[1:-1]
    _, self.fields['y-velocity'] = self.get_velocity(x_v, y_v, time)
    x_p, y_p = x_v, y_u
    self.fields['pressure'] = self.get_pressure(x_p, y_p, time)

  def mapped_meshgrid(self, x, y):
    """
    Maps the grid to a $[0, 2\pi]x[0, 2\pi]$ domain and returns the mesh-grid.

    Parameters
    ----------
    x, y: 1D numpy arrays of floats
      Stations along a grid-line.

    Returns
    -------
    X, Y: numpy meshgrid
      The mesh-grid.
    """
    X1, X2 = 0.0, 2.0 * numpy.pi
    x = (X1
         + (X2 - X1)
         * (x - self.bottom_left[0])
         / (self.top_right[0] - self.bottom_left[0]))
    y = (X1
         + (X2 - X1)
         * (y - self.bottom_left[1])
         / (self.top_right[1] - self.bottom_left[1]))
    return numpy.meshgrid(x, y)

  def get_velocity(self, x, y, time):
    """
    Computes the analytical solution of the velocity field.

    Parameters
    ----------
    x, y: 1D numpy arrays of floats
      Nodal stations along each direction.
    time: float
      The time.

    Returns
    -------
    ux, uy: Field objects
      The velocity components.
    """
    X, Y = self.mapped_meshgrid(x, y)
    return (Field(label='x-velocity',
                  x=x, y=y,
                  values=(1.0 - 2.0 * numpy.cos(X - 2.0 * numpy.pi * time)
                                    * numpy.sin(Y - 2.0 * numpy.pi * time))),
            Field(label='y-velocity',
                  x=x, y=y,
                  values=(1.0 + 2.0 * numpy.sin(X - 2.0 * numpy.pi * time)
                                    * numpy.cos(Y - 2.0 * numpy.pi * time))))

  def get_pressure(self, x, y, time):
    """Computes the analytical solution of the pressure field.

    Parameters
    ----------
    x, y: 1D numpy arrays of floats
      Nodal stations along each direction.
    time: float
      The time.

    Returns
    -------
    p: Field object
      The pressure field.
    """
    X, Y = self.mapped_meshgrid(x, y)
    return Field(label='pressure',
                 x=x, y=y,
                 values=(- numpy.cos(2.0 * (X - 2.0 * numpy.pi * time))
                         - numpy.cos(2.0 * (Y - 2.0 * numpy.pi * time))))

  def plot_fields(self, time_step,
                  view=[float('-inf'), float('-inf'),
                        float('inf'), float('inf')],
                  save_directory=os.path.join(os.getcwd(), 'images'),
                  save_name='analytical',
                  fmt='png',
                  dpi=100):
    """
    Plots the velocity and pressure fields.

    Parameters
    ----------
    time_step: integer
      Index used to as a suffix for the file names.
    view: 4-list of floats, optional
      Bottom-left and top-right coordinates of the view to plot;
      default: entire domain.
    save_directory: string, optional
      Directory where to save the figures;
      default: '<current directory>/images'.
    save_name: string, optional
      Prefix of the folder name that will contain the files;
      default: 'analytical'.
    fmt: string, optional
      Format of the files to save;
      default: 'png'.
    dpi: integer, optional
      Dots per inch (resolution);
      default: 100.
    """
    for name, field in self.fields.items():
      self.fields[name].time_step = time_step
      field.plot_contour(view=view,
                         save_directory=save_directory,
                         save_name=save_name,
                         fmt=fmt,
                         spi=100)

  def write_fields_petsc_format(self, x, y, time,
                                periodic_directions=None,
                                save_directory=None):
    """
    Computes and writes velocity and pressure fields into PETSc-readable files.
    The files are saved in the sub-folder 0000000.

    Parameters
    ----------
    x, y: 1D numpy arrays of floats
      Nodal stations along each direction.
    time: float
      The time at which the solution is computed.
    periodic_directions: list of strings
      Directions with periodic condition at the ends;
      default: None
    save_directory: string
      Directory of the simulation;
      default: None.
    """
    # create flux fields on staggered grid
    n_xu = x.size - (1 if 'x' in periodic_directions else 2)
    xu, yu = x[1: n_xu + 1], 0.5 * (y[:-1] + y[1:])
    qx = (self.get_velocity(xu, yu, time)[0].values
          * numpy.outer(y[1:] - y[:-1], numpy.ones(n_xu)))
    n_yv = y.size - (1 if 'y' in periodic_directions else 2)
    xv, yv = 0.5 * (x[:-1] + x[1:]), y[1: n_yv + 1]
    qy = (self.get_velocity(xv, yv, time)[1].values
          * numpy.outer(numpy.ones(n_yv), x[1:] - x[:-1]))
    # create directory where to save files
    if not save_directory:
      save_directory = os.path.join(os.getcwd(), '0000000')
    if not os.path.isdir(save_directory):
      os.makedirs(save_directory)
    sys.path.append(os.path.join(os.environ['PETSC_DIR'],
                                 'bin'))  # compatible petsc-3.7
    sys.path.append(os.path.join(os.environ['PETSC_DIR'],
                                 'bin',
                                 'pythonscripts'))  # compatible petsc-3.5
    import PetscBinaryIO
    # write fluxes
    vec = qx.flatten().view(PetscBinaryIO.Vec)
    file_path = os.path.join(save_directory, 'qx.dat')
    print('[info] writing fluxes in x-direction in file ...')
    PetscBinaryIO.PetscBinaryIO().writeBinaryFile(file_path, [vec, ])
    vec = qy.flatten().view(PetscBinaryIO.Vec)
    file_path = os.path.join(save_directory, 'qy.dat')
    print('[info] writing fluxes in y-direction in file ...')
    PetscBinaryIO.PetscBinaryIO().writeBinaryFile(file_path, [vec, ])
    # write pressure -- pressure field set to zero everywhere
    sizes = (y.size - 1, x.size - 1)
    vec = numpy.zeros(sizes).flatten().view(PetscBinaryIO.Vec)
    file_path = os.path.join(save_directory, 'phi.dat')
    print('[info] writing pressure in file ...')
    PetscBinaryIO.PetscBinaryIO().writeBinaryFile(file_path, [vec, ])
