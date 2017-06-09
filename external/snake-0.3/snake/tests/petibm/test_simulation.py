"""
Tests for the class `PetIBMSimulation`.
"""

import os
import sys
import shutil
import unittest
import numpy

if sys.version_info < (2, 8):
  sys.path.append(os.path.join(os.environ['PETSC_DIR'], 'bin'))
  import PetscBinaryIO

from snake.petibm.simulation import PetIBMSimulation


atol = 1.0E-12


class PetIBMSimulationTest(unittest.TestCase, PetIBMSimulation):
  def __init__(self, *args, **kwargs):
    super(PetIBMSimulationTest, self).__init__(*args, **kwargs)
    self.generate_stubs()

  def generate_stubs(self):
    self.directory = 'data'
    self.grid = numpy.array([numpy.linspace(0.0, 10.0, 11),
                             numpy.linspace(-1.0, 1.0, 101)])

  def test_read_grid(self):
    x, y = self.grid
    file_name = 'grid_test.txt'
    with open(file_name, 'w') as outfile:
      # write number of cells
      outfile.write('{}\t{}\n'.format(x.size - 1, y.size - 1))
    with open(file_name, 'ab') as outfile:
      # write cell-boundaries in x-direction
      numpy.savetxt(outfile, x, fmt='%.18f')
      # write cell-boundaries in y-direction
      numpy.savetxt(outfile, y, fmt='%.18f')
    self.read_grid(file_path=file_name)
    assert numpy.allclose(x, self.grid[0], atol=atol)
    assert numpy.allclose(y, self.grid[1], atol=atol)
    os.remove(file_name)

  def test_write_grid(self):
    x, y = self.grid
    file_name = 'grid_test.txt'
    self.write_grid(file_name, fmt='%0.18f')
    self.read_grid(file_path=file_name)
    assert numpy.allclose(x, self.grid[0], atol=atol)
    assert numpy.allclose(y, self.grid[1], atol=atol)
    os.remove(file_name)

  def test_read_forces(self):
    self.read_forces()
    assert len(self.forces) == 2
    assert numpy.allclose(self.forces[0].times, self.forces[1].times,
                          atol=atol)

  @unittest.skipIf(sys.version_info >= (3, 0), 'Not supported with Python 3')
  def test_read_fluxes(self):
    x, y = self.grid
    # create flux fields on staggered grid
    xu, yu = x[1: -1], 0.5 * (y[:-1] + y[1:])
    qx_ref = numpy.random.rand(yu.size, xu.size)
    xv, yv = 0.5 * (x[:-1] + x[1:]), y[1: -1]
    qy_ref = numpy.random.rand(yv.size, xv.size)
    # create solution folder
    directory = os.path.join(self.directory, '0000000')
    if not os.path.isdir(directory):
      os.makedirs(directory)
    # write fluxes
    vec = qx_ref.flatten().view(PetscBinaryIO.Vec)
    file_path = os.path.join(directory, 'qx.dat')
    PetscBinaryIO.PetscBinaryIO().writeBinaryFile(file_path, [vec, ])
    vec = qy_ref.flatten().view(PetscBinaryIO.Vec)
    file_path = os.path.join(directory, 'qy.dat')
    PetscBinaryIO.PetscBinaryIO().writeBinaryFile(file_path, [vec, ])
    # read fluxes
    qx, qy = self.read_fluxes(0)
    assert numpy.allclose(qx_ref, qx.values, atol=atol)
    assert numpy.allclose(xu, qx.x, atol=atol)
    assert numpy.allclose(yu, qx.y, atol=atol)
    assert numpy.allclose(qy_ref, qy.values, atol=atol)
    assert numpy.allclose(xv, qy.x, atol=atol)
    assert numpy.allclose(yv, qy.y, atol=atol)
    shutil.rmtree(directory)

  @unittest.skipIf(sys.version_info >= (3, 0), 'Not supported with Python 3')
  def test_read_pressure(self):
    x, y = self.grid
    # create pressure field
    xp, yp = 0.5 * (x[:-1] + x[1:]), 0.5 * (y[:-1] + y[1:])
    p_ref = numpy.random.rand(yp.size, xp.size)
    # create solution folder
    directory = os.path.join(self.directory, '0000000')
    if not os.path.isdir(directory):
      os.makedirs(directory)
    # write pressure
    vec = p_ref.flatten().view(PetscBinaryIO.Vec)
    file_path = os.path.join(directory, 'phi.dat')
    PetscBinaryIO.PetscBinaryIO().writeBinaryFile(file_path, [vec, ])
    # read pressure
    p = self.read_pressure(0)
    assert numpy.allclose(p_ref, p.values, atol=atol)
    assert numpy.allclose(xp, p.x, atol=atol)
    assert numpy.allclose(yp, p.y, atol=atol)
    shutil.rmtree(directory)


if __name__ == '__main__':
  unittest.main()
