"""
Implementation of the class `PetIBMSimulation`.
"""

import os
import sys
import struct

import numpy

try:
  sys.path.append(os.path.join(os.environ['PETSC_DIR'], 'bin'))
  import PetscBinaryIO
except:
  pass

from ..barbaGroupSimulation import BarbaGroupSimulation
from ..field import Field
from ..force import Force


class PetIBMSimulation(BarbaGroupSimulation):
  """
  Contains info about a PetIBM simulation.
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
      default: present working directory.
    """
    super(PetIBMSimulation, self).__init__(software='petibm',
                                           description=description,
                                           directory=directory,
                                           **kwargs)

  def read_grid(self, file_path=None):
    """
    Reads the grid from the file containing the grid nodes.

    Parameters
    ----------
    file_path: string, optional
      Path of the file containing grid-node stations along each direction;
      default: None.
    """
    print('[info] reading the grid ...')
    if not file_path:
      file_path = os.path.join(self.directory, 'grid.dat')
      if not os.path.exists(file_path):
        file_path = os.path.join(self.directory, 'grid.txt')
    # test if file written in binary format
    textchars = bytearray({7, 8, 9, 10, 12, 13, 27}
                          | set(range(0x20, 0x100)) - {0x7f})
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
        self.grid = numpy.array([x, y])
    else:
      with open(file_path, 'r') as infile:
        n_cells = numpy.array([int(n)
                               for n in infile.readline().strip().split()])
        coords = numpy.loadtxt(infile, dtype=numpy.float64)
      self.grid = numpy.array(numpy.split(coords,
                                          numpy.cumsum(n_cells[:-1] + 1)))
    if self.grid.size == 2:
      print('\tgrid-size: {}x{}'.format(self.grid[0].size - 1,
                                        self.grid[1].size - 1))
    elif self.grid.size == 3:
      print('\tgrid-size: {}x{}x{}'.format(self.grid[0].size - 1,
                                           self.grid[1].size - 1,
                                           self.grid[2].size - 1))

  def write_grid(self, file_path, fmt='%0.16g'):
    """
    Writes the stations along a gridline in each direction into a file.

    Parameters
    ----------
    file_path: string
      Path of the file to write into.
    fmt: string, optional
      Format to use for the stations;
      default: '%0.16g'.
    """
    with open(file_path, 'w') as outfile:
      if self.grid.size == 3:
        outfile.write('{}\t{}\t{}\n'.format(self.grid[0].size - 1,
                                            self.grid[1].size - 1,
                                            self.grid[2].size - 1))
      else:
        outfile.write('{}\t{}\n'.format(self.grid[0].size - 1,
                                        self.grid[1].size - 1))
    with open(file_path, 'ab') as outfile:
        numpy.savetxt(outfile, numpy.c_[self.grid[0]], fmt=fmt)
        numpy.savetxt(outfile, numpy.c_[self.grid[1]], fmt=fmt)
        if self.grid.size == 3:
          numpy.savetxt(outfile, numpy.c_[self.grid[2]], fmt=fmt)

  def read_forces(self, file_path=None, labels=None):
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
    """
    if not file_path:
      file_path = os.path.join(self.directory, 'forces.txt')
    print('[info] reading forces ...'),
    with open(file_path, 'r') as infile:
      data = numpy.loadtxt(infile, dtype=numpy.float64, unpack=True)
    times = data[0]
    if not labels:
      labels = ['f_x', 'f_z', 'f_z']  # default labels
    self.forces = []
    for index, values in enumerate(data[1:]):
      self.forces.append(Force(times, values, label=labels[index]))
    print('done')

  def read_fluxes(self, time_step,
                  periodic_directions=[],
                  directory=None,
                  **kwargs):
    """
    Reads the flux fields at a given time-step.

    Parameters
    ----------
    time_step: integer
      Time-step at which the field will be read.
    periodic_directions: list of strings, optional
      Directions that have periodic boundary conditions;
      default: [].
    directory: string, optional
      Directory where are saved the flux fields;
      default: None (defined as <simulation-directory>/<time-step>).

    Returns
    -------
    qx, qy, qz: Field objects
      Fluxes in the x-, y-, and z-directions.
    """
    print('[time-step {}] reading fluxes from files ...'.format(time_step)),
    dim3 = (len(self.grid) == 3)
    # directory with numerical solution
    if not directory:
      directory = os.path.join(self.directory, '{:0>7}'.format(time_step))
    # read grid-stations and fluxes
    x, y = self.grid[:2]
    nx, ny = x.size - 1, y.size - 1
    qx_file_path = os.path.join(directory, 'qx.dat')
    qx = PetscBinaryIO.PetscBinaryIO().readBinaryFile(qx_file_path)[0]
    qy_file_path = os.path.join(directory, 'qy.dat')
    qy = PetscBinaryIO.PetscBinaryIO().readBinaryFile(qy_file_path)[0]
    if dim3:
      z = self.grid[2]
      nz = z.size - 1
      qz_file_path = os.path.join(directory, 'qz.dat')
      qz = PetscBinaryIO.PetscBinaryIO().readBinaryFile(qz_file_path)[0]
    # create flux Field objects in staggered arrangement
    # reshape fluxes in multi-dimensional arrays
    if dim3:
      qx = qx.reshape((nz, ny, (nx if 'x' in periodic_directions else nx - 1)))
      qx = qx[:, :, :(-1 if 'x' in periodic_directions else None)]
      qx = Field(label='x-flux',
                 time_step=time_step,
                 x=x[1:-1],
                 y=0.5 * (y[:-1] + y[1:]),
                 z=0.5 * (z[:-1] + z[1:]),
                 values=qx)
      qy = qy.reshape((nz, (ny if 'y' in periodic_directions else ny - 1), nx))
      qy = qy[:, :(-1 if 'y' in periodic_directions else None), :]
      qy = Field(label='y-flux',
                 time_step=time_step,
                 x=0.5 * (x[:-1] + x[1:]),
                 y=y[1:-1],
                 z=0.5 * (z[:-1] + z[1:]),
                 values=qy)
      qz = qz.reshape(((nz if 'z' in periodic_directions else nz - 1), ny, nx))
      qz = qz[:(-1 if 'z' in periodic_directions else None), :, :]
      qz = Field(label='z-flux',
                 time_step=time_step,
                 x=0.5 * (x[:-1] + x[1:]),
                 y=0.5 * (y[:-1] + y[1:]),
                 z=z[1:-1],
                 values=qz)
      print('done')
      return qx, qy, qz
    else:
      qx = qx.reshape((ny, (nx if 'x' in periodic_directions else nx - 1)))
      qx = qx[:, :(-1 if 'x' in periodic_directions else None)]
      qx = Field(label='x-flux',
                 time_step=time_step,
                 x=x[1:-1],
                 y=0.5 * (y[:-1] + y[1:]),
                 values=qx)
      qy = qy.reshape(((ny if 'y' in periodic_directions else ny - 1), nx))
      qy = qy[:(-1 if 'y' in periodic_directions else None), :]
      qy = Field(label='y-flux',
                 time_step=time_step,
                 x=0.5 * (x[:-1] + x[1:]),
                 y=y[1:-1],
                 values=qy)
      print('done')
      return qx, qy

  def read_pressure(self, time_step, directory=None, **kwargs):
    """
    Reads the pressure field from file given the time-step.

    Parameters
    ----------
    time_step: integer
      Time-step at which the field will be read.
    directory: string, optional
      Directory where is saved the pressure field;
      default: None (set to <simulation-directory>/<time-step>).

    Returns
    -------
    p: Field object
      The pressure field.
    """
    print('[time-step {}] reading pressure field ...'.format(time_step)),
    dim3 = (len(self.grid) == 3)
    # get grid stations and number of cells along each direction
    x, y = self.grid[:2]
    nx, ny = x.size - 1, y.size - 1
    if dim3:
      z = self.grid[2]
      nz = z.size - 1
    # directory with numerical solution
    if not directory:
      directory = os.path.join(self.directory, '{:0>7}'.format(time_step))
    # read pressure
    phi_file_path = os.path.join(directory, 'phi.dat')
    p = PetscBinaryIO.PetscBinaryIO().readBinaryFile(phi_file_path)[0]
    # set pressure Field object
    if dim3:
      p = Field(label='pressure',
                time_step=time_step,
                x=0.5 * (x[:-1] + x[1:]),
                y=0.5 * (y[:-1] + y[1:]),
                z=0.5 * (z[:-1] + z[1:]),
                values=p.reshape((nz, ny, nx)))
    else:
      p = Field(label='pressure',
                time_step=time_step,
                x=0.5 * (x[:-1] + x[1:]),
                y=0.5 * (y[:-1] + y[1:]),
                values=p.reshape((ny, nx)))
    print('done')
    return p
