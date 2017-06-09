"""
Tests for the class `CuIBMSimulation`.
"""

import unittest
import numpy

from snake.cuibm.simulation import CuIBMSimulation


atol = 1.0E-12

# reference values
dim = 2
nx, ny = 10, 13
xmin, xmax = -1.0, 1.0
ymin, ymax = -1.1, 1.1
nt = 100
last_time = 0.04
last_fx, last_fy = 1.3524, -9.76969E-09


class CuIBMSimulationTest(unittest.TestCase, CuIBMSimulation):
  def __init__(self, *args, **kwargs):
    super(CuIBMSimulationTest, self).__init__(*args, **kwargs)
    self.generate_stubs()

  def generate_stubs(self):
    self.directory = 'data'

  def test_read_grid(self):
    self.read_grid()
    assert len(self.grid) == dim
    assert self.grid[0].size == nx + 1
    assert self.grid[1].size == ny + 1
    assert abs(self.grid[0].min() - xmin) <= atol
    assert abs(self.grid[0].max() - xmax) <= atol
    assert abs(self.grid[1].min() - ymin) <= atol
    assert abs(self.grid[1].max() - ymax) <= atol

  def test_read_forces(self):
    self.read_forces()
    assert len(self.forces) == dim
    for force in self.forces:
      assert force.times.size == force.values.size == nt
    assert numpy.allclose(self.forces[0].times, self.forces[1].times, atol)
    assert abs(self.forces[0].times[-1] - last_time) <= atol
    assert abs(self.forces[0].values[-1] - last_fx) <= atol
    assert abs(self.forces[1].values[-1] - last_fy) <= atol

  def test_read_fluxes(self):
    self.read_grid()
    qx, qy = self.read_fluxes(nt)
    x, y = self.grid[0][1:-1], 0.5 * (self.grid[1][1:] + self.grid[1][:-1])
    assert qx.label == 'x-flux'
    assert qx.time_step == nt
    assert qx.x.size == nx - 1
    assert numpy.allclose(qx.x, x, atol)
    assert qx.y.size == ny
    assert numpy.allclose(qx.y, y, atol)
    assert qx.values.shape == (qx.y.size, qx.x.size)
    x, y = 0.5 * (self.grid[0][1:] + self.grid[0][:-1]), self.grid[1][1:-1]
    assert qy.label == 'y-flux'
    assert qy.time_step == nt
    assert qy.x.size == nx
    assert numpy.allclose(qy.x, x, atol)
    assert qy.y.size == ny - 1
    assert numpy.allclose(qy.y, y, atol)
    assert qy.values.shape == (qy.y.size, qy.x.size)

  def test_read_pressure(self):
    self.read_grid()
    p = self.read_pressure(nt)
    x = 0.5 * (self.grid[0][1:] + self.grid[0][:-1])
    y = 0.5 * (self.grid[1][1:] + self.grid[1][:-1])
    assert p.label == 'pressure'
    assert p.time_step == nt
    assert p.x.size == nx
    assert numpy.allclose(p.x, x, atol)
    assert p.y.size == ny
    assert numpy.allclose(p.y, y, atol)
    assert p.values.shape == (p.y.size, p.x.size)


if __name__ == '__main__':
  unittest.main()
