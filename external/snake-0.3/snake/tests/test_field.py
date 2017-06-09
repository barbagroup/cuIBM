"""
Tests for the class `Field`.
"""

import unittest
import numpy

from snake.field import Field


atol = 1.0E-12


class FieldTest(unittest.TestCase, Field):
  def __init__(self, *args, **kwargs):
    super(FieldTest, self).__init__(*args, **kwargs)
    self.generate_stubs()

  def generate_stubs(self):
    self.label = 'test'
    self.time_step = 0
    self.x = numpy.linspace(0.0, 10.0, 9 * 4)
    self.y = numpy.linspace(-1.0, 1.0, 9 * 5)
    self.values = numpy.random.rand(self.y.size, self.x.size)

  def test_restrict(self):
    field1 = self.restrict(self.x, self.y, atol=atol)
    field2 = self.restrict(self.x[::3], self.y[::3], atol=atol)
    field3 = self.restrict(self.x[::9], self.y[::9], atol=atol)
    assert numpy.allclose(field1.x, self.x, atol=atol)
    assert numpy.allclose(field1.y, self.y, atol=atol)
    assert numpy.allclose(field1.values, self.values, atol=atol)
    assert numpy.allclose(field2.x, self.x[::3], atol=atol)
    assert numpy.allclose(field2.y, self.y[::3], atol=atol)
    assert numpy.allclose(field2.values, self.values[::3, ::3], atol=atol)
    assert numpy.allclose(field3.x, self.x[::9], atol=atol)
    assert numpy.allclose(field3.y, self.y[::9], atol=atol)
    assert numpy.allclose(field3.values, self.values[::9, ::9], atol=atol)

  def test_get_difference(self):
    assert (self.get_difference(self, self.x, self.y, norm='L2') == 0.0)
    assert (self.get_difference(self, self.x, self.y, norm='Linf') == 0.0)

  def test_subtract(self):
    subtracted_field = self.subtract(self)
    assert numpy.count_nonzero(subtracted_field.values) == 0

  def test_get_horizontal_gridline_values(self):
    index = 4
    x, u = self.get_horizontal_gridline_values(self.y[index])
    assert numpy.allclose(x, self.x, atol=atol)
    assert numpy.allclose(u, self.values[index, :], atol=atol)
    y = 0.5 * (self.y[index] + self.y[index + 1])
    x, u = self.get_horizontal_gridline_values(y)
    assert numpy.allclose(x, self.x, atol=atol)
    assert numpy.allclose(u, 0.5 * (self.values[index, :] +
                                    self.values[index + 1, :]),
                          atol=atol)

  def test_get_vertical_gridline_values(self):
    index = 4
    y, v = self.get_vertical_gridline_values(self.x[index])
    assert numpy.allclose(y, self.y, atol=atol)
    assert numpy.allclose(v, self.values[:, index], atol=atol)
    x = 0.5 * (self.x[index] + self.x[index + 1])
    y, v = self.get_vertical_gridline_values(x)
    assert numpy.allclose(y, self.y, atol=atol)
    assert numpy.allclose(v, 0.5 * (self.values[:, index] +
                                    self.values[:, index + 1]),
                          atol=atol)

  def test_get_gridline_values(self):
    index = 4
    y = 0.5 * (self.y[index] + self.y[index + 1])
    x, u = self.get_gridline_values(y=y)
    assert numpy.allclose(x, self.x, atol=atol)
    assert numpy.allclose(u, 0.5 * (self.values[index, :] +
                                    self.values[index + 1, :]),
                          atol=atol)
    x = 0.5 * (self.x[index] + self.x[index + 1])
    y, v = self.get_gridline_values(x=x)
    assert numpy.allclose(y, self.y, atol=atol)
    assert numpy.allclose(v, 0.5 * (self.values[:, index] +
                                    self.values[:, index + 1]),
                          atol=atol)


if __name__ == '__main__':
  unittest.main()
