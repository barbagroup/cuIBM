"""
Tests functions of the module `convergence`.
"""

import unittest
import numpy

from snake.field import Field
from snake import convergence


atol = 1.0E-12


class ConvergenceTest(unittest.TestCase):
  def __init__(self, *args, **kwargs):
    super(ConvergenceTest, self).__init__(*args, **kwargs)
    self.generate_stubs()

  def generate_stubs(self):
    nx, ny = 11, 21
    ratio = 3
    # grid to be used as a mask
    self.grid = [numpy.random.random(nx), numpy.random.random(ny)]
    # create fields
    self.fine = Field(x=numpy.ones(nx * ratio**2),
                      y=numpy.ones(ny * ratio**2),
                      values=numpy.random.rand(ny * ratio**2, nx * ratio**2),
                      label='fine')
    self.fine.x[::ratio**2] = self.grid[0][:]
    self.fine.y[::ratio**2] = self.grid[1][:]
    self.medium = Field(x=numpy.ones(nx * ratio),
                        y=numpy.ones(ny * ratio),
                        values=self.fine.values[::ratio, ::ratio] + 1.0,
                        label='medium')
    self.medium.x[::ratio] = self.grid[0][:]
    self.medium.y[::ratio] = self.grid[1][:]
    self.coarse = Field(x=self.grid[0],
                        y=self.grid[1],
                        values=(self.fine.values[::ratio**2, ::ratio**2] +
                                (1.0 + ratio)),
                        label='coarse')
    self.coarse2 = Field(x=self.grid[0],
                         y=self.grid[1],
                         values=numpy.random.rand(*self.coarse.values.shape),
                         label='coarse2')
    self.ratio = ratio

  def test_get_observed_order(self):
    p = convergence.get_observed_order(self.coarse, self.coarse2, self.coarse,
                                       self.ratio, self.grid, order=None)
    assert abs(p - 0.0) <= atol
    p = convergence.get_observed_order(self.coarse, self.coarse2, self.coarse,
                                       self.ratio, self.grid, order=numpy.inf)
    assert abs(p - 0.0) <= atol
    p = convergence.get_observed_order(self.coarse, self.medium, self.fine,
                                       self.ratio, self.grid, order=None)
    assert abs(p - 1.0) <= atol
    p = convergence.get_observed_order(self.coarse, self.medium, self.fine,
                                       self.ratio, self.grid, order=numpy.inf)
    assert abs(p - 1.0) <= atol


if __name__ == '__main__':
  unittest.main()
