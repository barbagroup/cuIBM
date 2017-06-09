"""
Implementation of the class `DecayingVortices`.
"""

import os
import sys
import math

import numpy

from ..field import Field


class DecayingVortices(object):
  """
  Analytical plug-in for the decaying vortices case.
  """

  def __init__(self, x, y, time, Re, amplitude):
    """
    Computes the velocities and pressure fields on a given grid.

    Parameters
    ----------
    x, y: 1D numpy arrays of floats
      Contains the stations along a gridline in each direction.
    time: float
      Time at which the analytical solution will be computed.
    Re: float
      Reynolds number.
    amplitude: float
      Amplitude of the Taylor-Green vortex.
    """
    self.bottom_left, self.top_right = [x[0], y[0]], [x[-1], y[-1]]
    self.fields = {}
    x_u, y_u = x[1:-1], 0.5 * (y[:-1] + y[1:])
    self.fields['x-velocity'], _ = self.get_velocity(x_u, y_u,
                                                     time,
                                                     Re,
                                                     amplitude)
    x_v, y_v = 0.5 * (x[:-1] + x[1:]), y[1:-1]
    _, self.fields['y-velocity'] = self.get_velocity(x_v, y_v,
                                                     time,
                                                     Re,
                                                     amplitude)
    x_p, y_p = x_v, y_u
    self.fields['pressure'] = self.get_pressure(x_p, y_p,
                                                time,
                                                Re)
    fluxes = self.get_flux_from_velocity(x, y)
    self.fields['x-flux'], self.fields['y-flux'] = fluxes

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

  def get_velocity(self, x, y, time, Re, amplitude):
    """
    Computes the analytical solution of the velocity field.

    Parameters
    ----------
    x, y: 1D numpy arrays of floats
      Nodal stations along each direction.
    time: float
      The time.
    Re: float
      The Reynolds number.
    amplitude: float
      amplitude of the vortices.

    Returns
    -------
    ux, uy: Field objects
      The velocity components.
    """
    X, Y = self.mapped_meshgrid(x, y)
    return (Field(label='x-velocity',
                  x=x, y=y,
                  values=(- amplitude * numpy.cos(X) * numpy.sin(Y)
                          * math.exp(-2.0 * (2.0 * numpy.pi)**2 * time / Re))),
            Field(label='y-velocity',
                  x=x, y=y,
                  values=(amplitude * numpy.sin(X) * numpy.cos(Y)
                          * math.exp(-2.0 * (2.0 * numpy.pi)**2 * time / Re))))

  def get_flux_from_velocity(self, x, y):
    """
    Computes the flux fields from the velocity fields given the stations
    along a gridline in each direction of a 2D structured Cartesian grid.

    Parameters
    ----------
    x, y: 1D numpy arrays of floats
      Nodal stations along each direction.

    Returns
    -------
    qx, qy: Field objects
      The flux fields.
    """
    dx, dy = x[1:] - x[:-1], y[1:] - y[:-1]
    return (Field(label='x-flux',
                  x=self.fields['x-velocity'].x,
                  y=self.fields['x-velocity'].y,
                  values=self.fields['x-velocity'].values * dy[:, None]),
            Field(label='y-flux',
                  x=self.fields['y-velocity'].x,
                  y=self.fields['y-velocity'].y,
                  values=self.fields['y-velocity'].values * dx[None, :]))

  def get_pressure(self, x, y, time, Re):
    """Computes the analytical solution of the pressure field.

    Parameters
    ----------
    x, y: 1D numpy arrays of floats
      Nodal stations along each direction.
    time: float
      The time.
    Re: float
      The Reynolds number.

    Returns
    -------
    p: Field object
      The pressure field.
    """
    X, Y = self.mapped_meshgrid(x, y)
    return Field(label='pressure',
                 x=x, y=y,
                 values=(- 0.25 * (numpy.cos(2.0 * X) + numpy.cos(2.0 * Y))
                         * math.exp(-4.0 * (2.0 * numpy.pi)**2 * time / Re)))

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
      Directory where to save the images;
      default: '<current directory>/images'.
    save_name: string, optional
      Prefix of the folder name that will contain the .png files;
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
                         dpi=dpi)

  def write_fields_petsc_format(self, x, y, time, Re, amplitude,
                                periodic_directions=None,
                                save_directory=None):
    """
    Computes and writes velocity and pressure fields into PETSc-readable files.

    Parameters
    ----------
    x, y: 1D numpy arrays of floats
      Nodal stations along each direction.
    time: float
      The time at which the solution is computed.
    Re: float
      The Reynolds number.
    amplitude: float
      The amplitude of the vortices.
    periodic_directions: list of strings, optional
      Directions with periodic condition at the ends;
      default: None
    save_directory: string, optional
      Directory of the simulation;
      default: None.
    """
    # create flux fields on staggered grid
    n_xu = x.size - (1 if 'x' in periodic_directions else 2)
    xu, yu = x[1: n_xu + 1], 0.5 * (y[:-1] + y[1:])
    qx = (self.get_velocity(xu, yu, time, Re, amplitude)[0].values
          * numpy.outer(y[1:] - y[:-1], numpy.ones(n_xu)))
    n_yv = y.size - (1 if 'y' in periodic_directions else 2)
    xv, yv = 0.5 * (x[:-1] + x[1:]), y[1: n_yv + 1]
    qy = (self.get_velocity(xv, yv, time, Re, amplitude)[1].values
          * numpy.outer(numpy.ones(n_yv), x[1:] - x[:-1]))
    # create directory where to save files
    if not save_directory:
      save_directory = os.path.join(os.getcwd(), '0000000')
    if not os.path.isdir(save_directory):
      os.makedirs(save_directory)
    sys.path.append(os.path.join(os.environ['PETSC_DIR'],
                                 'bin'))  # compatible petsc-3.5
    sys.path.append(os.path.join(os.environ['PETSC_DIR'],
                                 'bin',
                                 'pythonscripts'))  # compatible petsc-3.7
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
