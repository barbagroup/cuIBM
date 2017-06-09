"""
Implementation of the class `GEOFile` used to generate a .geo file.
"""

import os

import numpy


class GEOFile(object):
  def __init__(self, name='mesh'):
    self.name = name
    self.domain = None
    self.body = None
    self.boxes = []

  def define_domain(self,
                    bottom_left=[-1.0, -1.0],
                    top_right=[1.0, 1.0],
                    n_inlet=10):
    self.domain = Domain(bottom_left, top_right, n_inlet)

  def define_body(self, coordinates_path, levels):
    self.body = Body(coordinates_path, self.domain, levels)

  def define_refinement_boxes(self, boxes_info=[]):
    n_boxes = len(boxes_info) / 5
    for i in range(n_boxes):
      bottom_left = boxes_info[i], boxes_info[i + 1]
      top_right = boxes_info[i + 2], boxes_info[i + 3]
      levels = boxes_info[i + 4]
      self.boxes = Box(i, bottom_left, top_right,
                       levels, self.domain.ref_length)

  def write(self, save_directory):
    file_path = os.path.join(save_directory, self.name + '.geo')
    with open(file_path, 'r') as outfile:
      outfile.write('# {}\n'.format(self.name))
    self.body.write_geo(file_path)
    self.domain.write_geo(file_path, counter=self.body.n)
    with open(file_path, 'a') as outfile:
      # plane surface
      outfile.write('// plane surface\n')
      outfile.write('Plane Surface(1) = {1, 2};\n')
      # physical volume
      outfile.write('// physical volume\n')
      outfile.write('Physical Volume(1) = {1};\n')
    for i, box in enumerate(self.boxes):
      box.write_geo(file_path, counter=i)


class Domain(object):
  def __init__(self,
               bottom_left=[-1.0, -1.0],
               top_right=[1.0, 1.0],
               n_inlet=10):
    self.bottom_left, self.top_right = bottom_left, top_right
    self.n_inlet = n_inlet
    self.ref_length = (top_right[1] - bottom_left[1]) / self.n_inlet

  def write_geo(self, outfile_path, mode='a', counter=0):
    with open(outfile_path, mode) as outfile:
      outfile.write('cl_exterior = {};\n'.format(self.ref_length))
      # points
      outfile.write('// domain points\n')
      outfile.write('Point({}) = {{}, {}, 0.0, {}};\n'
                    ''.format(counter + 1,
                              self.bottom_left[0],
                              self.bottom_left[1],
                              self.ref_length))
      outfile.write('Point({}) = {{}, {}, 0.0, {}};\n'
                    ''.format(counter + 2,
                              self.top_right[0],
                              self.bottom_left[1],
                              self.ref_length))
      outfile.write('Point({}) = {{}, {}, 0.0, {}};\n'
                    ''.format(counter + 3,
                              self.top_right[0],
                              self.top_right[1],
                              self.ref_length))
      outfile.write('Point({}) = {{}, {}, 0.0, {}};\n'
                    ''.format(counter + 4,
                              self.bottom_left[0],
                              self.top_right[1],
                              self.ref_length))
      # lines
      outfile.write('// domain lines\n')
      outfile.write('Line({0}) = {{0}, {1}};\n'.format(counter + 1,
                                                       counter + 2))
      outfile.write('Line({0}) = {{0}, {1}};\n'.format(counter + 2,
                                                       counter + 3))
      outfile.write('Line({0}) = {{0}, {1}};\n'.format(counter + 3,
                                                       counter + 4))
      outfile.write('Line({0}) = {{0}, {1}};\n'.format(counter + 4,
                                                       counter + 1))
      # line-loop
      outfile.write('// domain line-loop\n')
      outfile.write('Line Loop(2) = {{}};\n'
                    ''.format(', '.join([str(counter + i)
                                         for i in [1, 2, 3, 4]])))


class Body(object):
  def __init__(self, domain, levels=0, name='body'):
    self.name = name
    self.ref_length = domain.ref_length / 2.0**levels

  def read_coordinates_from_file(self, file_path):
    """
    Reads two-dimensional coordinates from input file.

    Parameters
    ----------
    file_path: string
      Path of the file containing the 2d coordinates.
    """
    with open(file_path, 'r') as infile:
      self.x, self.y = numpy.loadtxt(infile,
                                     dtype=float,
                                     delimiter='\t',
                                     skiprows=1,
                                     unpack=True)
    self.n = self.x.size
    # assume closed loop
    x, y = self.x, self.y
    self.lengths = numpy.append(numpy.sqrt((x[:-1] - x[1:])**2
                                           + (y[:-1] - y[1:])**2),
                                numpy.sqrt((x[0] - x[-1])**2
                                           + (y[0] - y[-1])**2))

  def write_geo(self, outfile_path, mode='a', counter=0):
    with open(outfile_path, mode) as outfile:
      outfile.write('cl_body = {};\n'.format(self.ref_length))
      # points
      outfile.write('// body points\n')
      for i in range(self.n):
        outfile.write('Point({}) = {{}, {}, 0.0};\n'.format(counter + i + 1,
                                                            self.x[i],
                                                            self.y[i]))
      # lines
      outfile.write('// body lines\n')
      for i in range(self.n - 1):
        outfile.write('Line({0}) = {{0}, {1}};\n'.format(counter + i + 1,
                                                         counter + i + 2))
      outfile.write('Line({0}) = {{0}, {1}};\n'.format(counter + self.n,
                                                       counter + 1))
      # line-loop
      outfile.write('// body line-loop\n')
      outfile.write('Line Loop(1) = {{}};\n'
                    ''.format(', '.join([str(counter + i + 1)
                                         for i in range(self.n)])))


class Box(object):
  def __init__(self, index,
               bottom_left=[-1.0, -1.0],
               top_right=[1.0, 1.0],
               levels=0,
               domain_ref_length=1.0):
    self.index = index
    self.bottom_left, self.top_right = bottom_left, top_right
    self.ref_length_in = domain_ref_length
    self.ref_left_out = domain_ref_length

  def write_geo(self, outfile_path, mode='a'):
    """
    Writes the box into a file.

    Parameters
    ---------_
    outfile_path: string
      Path of the file to write.
    mode: string, optional
      Writing-mode;
      default: 'a' (append).
    """
    with open(outfile_path, mode) as outfile:
      outfile.write('// box field {}\n'.format(self.index))
      outfile.write('Field[{}] = Box;\n'.format(self.index))
      outfile.write('Field[{}].VIn = {};\n'.format(self.index,
                                                   self.ref_length_in))
      outfile.write('Field[{}].VOut = {};\n'.format(self.index,
                                                    self.ref_length_out))
      outfile.write('Field[{}].XMin = {};\n'.format(self.index,
                                                    self.bottom_left[0]))
      outfile.write('Field[{}].XMax = {};\n'.format(self.index,
                                                    self.top_right[0]))
      outfile.write('Field[{}].YMin = {};\n'.format(self.index,
                                                    self.bottom_left[1]))
      outfile.write('Field[{}].YMax = {};\n'.format(self.index,
                                                    self.top_right[1]))
