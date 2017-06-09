"""
Implementation of the class `CuIBMSimulation`.
"""

import os
import struct

import numpy

from ..barbaGroupSimulation import BarbaGroupSimulation
from ..field import Field
from ..force import Force


class CuIBMSimulation(BarbaGroupSimulation):
  """
  Contains info about a cuIBM simulation.
  Inherits from the class BarbaGroupSimulation.
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
    **kwargs: dictionary
      Extra keyword-arguments.
    """
    super(CuIBMSimulation, self).__init__(software='cuibm',
                                          description=description,
                                          directory=directory,
                                          **kwargs)

  def read_grid(self, file_path=None):
    """
    Reads the computational grid from file.

    Parameters
    ----------
    file_path: string, optional
      Path of the file containing grid stations along each direction;
      default: None.
    """
    print('[info] reading grid ...')
    if not file_path:
      file_path = os.path.join(self.directory, 'grid')
    # test if file written in binary format
    textchars = (bytearray({7, 8, 9, 10, 12, 13, 27}
                           | set(range(0x20, 0x100)) - {0x7f}))
    is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))
    infile = open(file_path, 'rb')
    binary_format = is_binary_string(infile.read(1024))
    infile.close()
    if binary_format:
      with open(file_path, 'rb') as infile:
        # x-direction
        nx = struct.unpack('i', infile.read(4))[0]
        x = numpy.array(struct.unpack('d' * (nx + 1),
                                      infile.read(8 * (nx + 1))))
        # y-direction
        ny = struct.unpack('i', infile.read(4))[0]
        y = numpy.array(struct.unpack('d' * (ny + 1),
                                      infile.read(8 * (ny + 1))))
    else:
      with open(file_path, 'r') as infile:
        data = numpy.loadtxt(infile, dtype=numpy.float64)
        # x-direction
        nx = int(data[0])
        x, data = data[1:nx + 2], data[nx + 2:]
        # y-direction
        ny = int(data[0])
        y = data[1:]
    self.grid = x, y
    print('\tgrid-size: {}x{}'.format(x.size - 1, y.size - 1))

  def read_forces(self, file_path=None, labels=None, usecols=(0, 1, 2)):
    """
    Reads forces from files.

    Parameters
    ----------
    file_path: string, optional
      Path of the file containing the forces data;
      default: None.
    labels: list of strings, optional
      Label to give to each force that will be read from file;
      default: None
    usecols: tuple of integers, optional
      Index of each column to read in the forces file
      (including the time column);
      default: (0, 1, 2)
    """
    if not file_path:
      file_path = os.path.join(self.directory, 'forces')
    print('[info] reading forces ...')
    with open(file_path, 'r') as infile:
      data = numpy.loadtxt(infile,
                           dtype=numpy.float64,
                           usecols=usecols,
                           unpack=True)
    times = data[0]
    if not labels:
      labels = ['f_x', 'f_y']  # default labels
    self.forces = []  # reset forces if already present
    for index, values in enumerate(data[1:]):
      self.forces.append(Force(times, values, label=labels[index]))

  def read_fluxes(self, time_step, directory=None, **kwargs):
    """
    Reads the flux fields from file at a given time-step.

    Parameters
    ----------
    time_step: integer
      Time-step at which to read the fluxes.
    directory: string, optional
      Directory containing the saved time-step folders;
      default: None (will be <simulation-directory>/<time-step>).
    **kwargs: dictionary
      Extra keyword-arguments.

    Returns
    -------
    qx, qy: Field objects
      Fluxes in the x- and y-directions.
    """
    print('[time-step {}] reading fluxes from file ...'.format(time_step))
    # get grid-stations and number of cells along each direction
    x, y = self.grid
    nx, ny = x.size - 1, y.size - 1
    # read fluxes from file
    if not directory:
      directory = os.path.join(self.directory, '{:0>7}'.format(time_step))
    file_path = os.path.join(directory, 'q')
    # test if file written in binary format
    textchars = bytearray({7, 8, 9, 10, 12, 13, 27}
                          | set(range(0x20, 0x100)) - {0x7f})
    is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))
    infile = open(file_path, 'rb')
    binary_format = is_binary_string(infile.read(1024))
    infile.close()
    if binary_format:
      with open(file_path, 'rb') as infile:
        nq = struct.unpack('i', infile.read(4))[0]
        q = numpy.array(struct.unpack('d' * nq, infile.read(8 * nq)))
    else:
      with open(file_path, 'r') as infile:
        nq = int(infile.readline())
        q = numpy.loadtxt(infile, dtype=numpy.float64)
    # set flux Field objects
    qx = Field(label='x-flux',
               time_step=time_step,
               x=x[1:-1],
               y=0.5 * (y[:-1] + y[1:]),
               values=q[:(nx - 1) * ny].reshape(ny, nx - 1))
    offset = qx.values.size
    qy = Field(label='y-flux',
               time_step=time_step,
               x=0.5 * (x[:-1] + x[1:]),
               y=y[1:-1],
               values=q[offset:].reshape(ny - 1, nx))
    return qx, qy

  def read_pressure(self, time_step, directory=None, **kwargs):
    """
    Reads pressure field from solution file at given time-step.

    Parameters
    ----------
    time_step: integer
      Time-step at which to read the pressure field.
    directory: string, optional
      Directory containing the saved time-step folders;
      default: None (will be <simulation-directory>/<time-step>).
    **kwargs: dictionary
      Extra keyword-arguments.

    Returns
    -------
    p: Field object
      The pressure field.
    """
    print('[time-step {}] reading pressure from file ...'.format(time_step))
    # get info about mesh-grid
    x, y = self.grid
    nx, ny = x.size - 1, y.size - 1
    # read pressure from file
    if not directory:
      directory = os.path.join(self.directory, '{:0>7}'.format(time_step))
    file_path = os.path.join(directory, 'lambda')
    # test if file written in binary format
    textchars = bytearray({7, 8, 9, 10, 12, 13, 27}
                          | set(range(0x20, 0x100)) - {0x7f})
    is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))
    infile = open(file_path, 'rb')
    binary_format = is_binary_string(infile.read(1024))
    infile.close()
    if binary_format:
      with open(file_path, 'rb') as infile:
        nlambda = struct.unpack('i', infile.read(4))[0]
        p = numpy.array(struct.unpack('d' * nlambda,
                                      infile.read(8 * nlambda)))[:nx * ny]
    else:
      with open(file_path, 'r') as infile:
        nlambda = int(infile.readline())
        p = numpy.loadtxt(infile, dtype=numpy.float64)[:nx * ny]
    # set pressure Field object
    p = Field(label='pressure',
              time_step=time_step,
              x=0.5 * (x[:-1] + x[1:]),
              y=0.5 * (y[:-1] + y[1:]),
              values=p.reshape(ny, nx))
    return p
