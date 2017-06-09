"""
Implementation of the class `CartesianStructuredMesh` and its sub-classes.
"""

import os
import sys
import math
from decimal import Decimal

import numpy
import yaml


class Segment(object):
  """
  Contains information about a segment.
  """

  def __init__(self, data=None, vertices=None, mode=None):
    """
    Creates the segment vertices.

    Parameters
    ----------
    data: dictionary, optional
      YAML data;
      default: None.
    vertices: 1D array of floats, optional
      vertices along the segment;
      default: None.
    mode: string, optional
      Mode to use to create the segment;
      choices: ['cuibm'];
      default: None.
    """
    if data:
      self.create_from_yaml_data(data, mode=mode)
    elif numpy.all(vertices):
      self.create_from_vertices(vertices)
    self.nb_divisions = self.vertices.size - 1

  def create_from_vertices(self, vertices):
    """
    Stores vertices.
    """
    self.vertices = vertices

  def create_from_yaml_data(self, data, mode=None):
    """
    Creates vertices from provided YAML data.

    Parameters
    ----------
    data: dictionary
      YAML data.
    mode: string, optional
      Mode to use to create the gridline;
      choices: ['cuibm'];
      default: None.
    """
    self.start, self.end = data['start'], data['end']
    self.width = data['width']
    if mode == 'cuibm':
      self.stretch_ratio = data['stretchRatio']
    else:
      ratio = self.get_stretch_ratio(self.width,
                                     stretch_ratio=data['stretchRatio'],
                                     aspect_ratio=data['aspectRatio'],
                                     precision=data['precision'])
      self.stretch_ratio = ratio
    self.vertices = self.get_vertices(reverse=data['reverse'], mode=mode)

  def print_parameters(self):
    """
    Prints the parameters of the segment.
    """
    print('\n[info] printing parameters of segment ...')
    print('\tlimits: [{}, {}]'.format(self.start, self.end))
    print('\treal limits: [{}, {}]'.format(self.vertices[0],
                                           self.vertices[-1]))
    print('\tstretching ratio: {}'.format(self.stretch_ratio))
    print('\taspect ratio: {}'.format(self.aspect_ratio))
    print('\tnumber of divisions: {}\n'.format(self.nb_divisions))

  def get_vertices(self, reverse=False, mode=None):
    """
    Computes the vertices of the segment.

    Parameters
    ----------
    reverse: boolean, optional
      Set 'True' if you want to reverse the stretching order;
      default: False.
    mode: string, optional
      Mode to use to create the gridline;
      choices: ['cuibm'];
      default: None.

    Returns
    -------
    vertices: 1D array of floats
      The vertices of the segment.
    """
    length = abs(self.end - self.start)
    width = self.width
    ratio = self.stretch_ratio
    if abs(ratio - 1.0) < 1.0E-06:
      # uniform discretization
      self.aspect_ratio = 1.0
      if abs(int(round(length / width)) - length / width) > 1.0E-12:
        print('[error] uniform discretization: '
              'length of the segment should be a multiple of the width '
              'provided')
        sys.exit(-1)
      return numpy.arange(self.start, self.end + width / 2.0, width)
    # stretched discretization a la cuibm
    if mode == 'cuibm':
      n = 2
      h1 = length * (ratio - 1.0) / (ratio**n - 1.0)
      while h1 > width:
        n += 1
        h1 = length * (ratio - 1.0) / (ratio**n - 1.0)
      n -= 1
      h1 = length * (ratio - 1.0) / (ratio**n - 1.0)
      h2 = h1 * ratio**(n - 1)
      width0 = (h2 if reverse else h1)
      self.stretch_ratio = (1.0 / ratio if reverse else ratio)
      widths = numpy.empty(n, dtype=numpy.float64)
      widths[0], widths[1:] = width0, self.stretch_ratio
      widths = numpy.cumprod(widths)
      self.aspect_ratio = (widths[0] / widths[-1] if reverse
                           else widths[-1] / widths[0])
      return numpy.insert(self.start + numpy.cumsum(widths), 0, self.start)
    # stretched discretization
    n = int(round(math.log(1.0 - length / width * (1.0 - ratio))
                  / math.log(ratio)))
    widths = numpy.empty(n, dtype=numpy.float64)
    widths[0], widths[1:] = width, ratio
    widths = numpy.cumprod(widths)
    # compute the aspect ratio
    self.aspect_ratio = widths[-1] / widths[0]
    # return the vertices
    if reverse:
      # inverse the stretching ratio and reverse the widths
      self.stretch_ratio = 1.0 / ratio
      return numpy.insert(self.end - numpy.cumsum(widths), 0, self.end)[::-1]
    else:
      return numpy.insert(self.start + numpy.cumsum(widths), 0, self.start)

  def get_stretch_ratio(self, width,
                        stretch_ratio=1.0, aspect_ratio=1.0, precision=6):
    """
    Computes the optimal stretching ratio given a targeted stretching ratio
    or a targeted aspect ratio.

    Parameters
    ----------
    width: float
      Width of the first division.
    stretch_ratio: float, optional
      Targeted stretching ratio;
      default: 1.0.
    aspect_ratio: float, optional
      Targeted aspect ratio between the first and last divisions;
      default: 1.0.
    precision: integer, optional
      Precision of the optimal stretching ratio to compute;
      default: 6.

    Returns
    -------
    ratio: float
      The optimal stretching ratio.
    """
    # if stretching ratio provided
    if abs(stretch_ratio - 1.0) > 1.0E-06:
      return self.compute_optimal_stretch_ratio(width, stretch_ratio,
                                                precision=precision)
    # if aspect ratio provided
    elif abs(aspect_ratio - 1.0) > 1.0E-06:
      ratio = self.compute_stretch_ratio(width, aspect_ratio,
                                         precision=precision)
      return self.compute_optimal_stretch_ratio(width, ratio,
                                                precision=precision)
    # uniform discretization
    else:
      return 1.0

  def compute_stretch_ratio(self, width, aspect_ratio,
                            precision=6):
    """
    Computes the stretching ratio provided the a targeted aspect ratio
    between the first and last divisions of the segment.

    Parameters
    ----------
    width: float
      Width of the first division.
    aspect_ratio: float
      Targeted aspect ratio.
    precision: integer, optional
      Precision of the stretching ratio to compute;
      default: 6.

    Returns
    -------
    ratio: float
      The stretching ratio.
    """
    length = abs(self.end - self.start)
    current_precision = 1
    ratio = 2.0
    while current_precision < precision:
        n = int(round(math.log(1.0 - length / width * (1.0 - ratio))
                      / math.log(ratio)))
        candidate_aspect_ratio = ratio**(n - 1)
        if candidate_aspect_ratio < aspect_ratio:
            ratio += (0.1)**current_precision
            current_precision += 1
        else:
            ratio -= (0.1)**current_precision
    return ratio

  def compute_optimal_stretch_ratio(self, width, ratio,
                                    precision=6):
    """
    Computes the optimal stretching ratio provided a targeted one.

    Parameters
    ----------
    width: float
      Width of the first division of the segment.
    ratio: float
      Targeted stretching ratio.
    precision: integer, optional
      Precision of the stretching ratio to compute;
      default: 6.

    Returns
    -------
    ratio: float
      The optimal stretching ratio.
    """
    def geometric_sum(a, r, n):
      """
      Computes the sum of the geometric progression.
      """
      return a * (1.0 - r**n) / (1.0 - r)
    length = abs(self.end - self.start)
    precision_ratio = abs(Decimal(str(ratio)).as_tuple().exponent)
    while precision_ratio < precision:
      n = int(math.log(1.0 - (1.0 - ratio) * length / width)
              / math.log(ratio))
      deviation_inf = abs(length - geometric_sum(width, ratio, n))
      deviation_sup = abs(length - geometric_sum(width, ratio, n + 1))
      precision_ratio += 1
      if deviation_inf < deviation_sup:
        ratio += 0.1**precision_ratio
      else:
        ratio -= 0.1**precision_ratio
    return ratio

  def generate_yaml_info(self):
    """
    Generates a dictionary with segment's information ready for YAML.

    The dictionary contains the end, the number of divisions, and the
    stretching ratio of the segment.
    """
    info = {}
    info['end'] = self.end
    info['cells'] = self.nb_divisions
    info['stretchRatio'] = self.stretch_ratio
    return info


class GridLine(object):
  """
  Contains information about a gridline.
  """

  def __init__(self, data=None, vertices=None, label=None, mode=None):
    """
    Creates a gridline from provided YAML data or vertices.

    Parameters
    ----------
    data: dictionary, optional
      YAML data about the gridline;
      default: None.
    vertices: 1D array of floats, optional
      Vertices along the gridline;
      default: None.
    label: string, optional
      Label of the direction;
      default: None.
    mode: string, optional
      Mode to use to create the gridline;
      choices: ['cuibm'];
      default: None.
    """
    self.label = label
    self.segments = []
    if data:
      self.create_from_yaml_data(data, mode=mode)
    elif numpy.all(vertices):
      self.create_from_vertices(vertices)
    self.nb_divisions = sum(segment.nb_divisions for segment in self.segments)

  def create_from_vertices(self, vertices):
    """
    Defines the gridline from provided vertices as a single segment.

    Parameters
    ----------
    vertices: 1D array of floats
      The vertices along the gridline.
    """
    self.start, self.end = vertices[0], vertices[-1]
    self.segments.append(Segment(vertices=vertices))

  def create_from_yaml_data(self, data, mode=None):
    """
    Initializes the gridline parameters and computes its vertices.

    A gridline is defined as a sequence of uniform and/or stretched segments.

    Parameters
    ----------
    data: dictionary
      Parameters of the gridline in a YAML format.
    mode: string, optional
      Mode to use to create the gridline;
      choices: ['cuibm'];
      default: None.
    """
    self.label = data['direction']
    self.start = data['start']
    self.end = data['subDomains'][-1]['end']
    for index, node in enumerate(data['subDomains']):
      # store starting point
      start = (data['start'] if index == 0
               else data['subDomains'][index - 1]['end'])
      data['subDomains'][index]['start'] = start
      # set default parameters if not present
      if 'reverse' not in node.keys():
        data['subDomains'][index]['reverse'] = False
      if 'precision' not in node.keys():
        data['subDomains'][index]['precision'] = 6
      if 'aspectRatio' not in node.keys():
        data['subDomains'][index]['aspectRatio'] = 1.0
      if 'stretchRatio' not in node.keys():
        data['subDomains'][index]['stretchRatio'] = 1.0
      # create a segment
      self.segments.append(Segment(data=node, mode=mode))

  def get_vertices(self, precision=6):
    """
    Gets the vertices removing the repeated values at boundaries
    between consecutive segments.

    Parameters
    ----------
    precision: integer, optional
      Precision used to round the vertices so we can remove repeated values;
      default: 6.

    Returns
    -------
    vertices: 1D array of floats
      The vertices along the gridline.
    """
    return numpy.unique(numpy.concatenate(([numpy.round(segment.vertices,
                                                        precision)
                                            for segment in self.segments])))

  def print_parameters(self):
    """
    Prints parameters of the gridline.
    """
    print('[info] printing parameters of the gridline '
          '{} ...'.format(self.label))
    vertices = self.get_vertices()
    print('\tlimits: [{}, {}]'.format(self.start, self.end))
    print('\treal limits: [{}, {}]'.format(vertices[0], vertices[-1]))
    print('\tnumber of divisions: {}'.format(self.nb_divisions))
    for segment in self.segments:
      segment.print_parameters()

  def generate_yaml_info(self):
    """
    Generates a dictionary with gridline's information ready for YAML.

    The dictionary contains the direction, the start, and segments information
    of the gridline.
    """
    info = {}
    info['direction'] = self.label
    info['start'] = self.start
    info['subDomains'] = []
    for segment in self.segments:
      info['subDomains'].append(segment.generate_yaml_info())
    return info


class CartesianStructuredMesh(object):
  """
  Contains info related to a Cartesian structured mesh (stretched or uniform).
  """

  def __init__(self):
    """
    Instantiates an empty mesh.
    """
    self.gridlines = []

  def get_number_cells(self):
    """
    Gets the number of divisions along each gridline and the total number of
    cells.
    """
    nb_divisions = []
    for gridline in self.gridlines:
      nb_divisions.append(gridline.nb_divisions)
    total = 1
    for nb_division in nb_divisions:
      total *= nb_division
    return total, nb_divisions

  def create(self, data, mode=None):
    """
    Creates the gridlines.

    Parameters
    ----------
    data: list of dictionaries
      Contains YAML information about the gridlines.
    mode: string, optional
      Mode to use to generate the grid;
      choices: ['cuibm'];
      default: None.
    """
    print('[info] generating Cartesian grid ...')
    for node in data:
      self.gridlines.append(GridLine(data=node, mode=mode))

  def print_parameters(self):
    """Prints parameters of the Cartesian structured mesh."""
    print('[info] printing parameters of the Cartesian structured mesh ...')
    nb_cells, _ = self.get_number_cells()
    print('\tnumber of cells: {}'.format(nb_cells))
    for gridline in self.gridlines:
      gridline.print_parameters()

  def write(self, file_path, precision=6, direction='all'):
    """
    Writes the gridlines into a file.

    The first line of the file contains the number of divisions
    along each gridline.
    Then the vertices along each gridline are written in single column
    starting with the first gridline.

    The last argument allows to output the gridline in one given direction.
    In that case, the stations are listed in a single column.

    Parameters
    ----------
    file_path: string
      Path of the file to write in.
    precision: integer, optional
      Precision at which the vertices will be written;
      default: 6.
    direction: string, optional
      Which gridline to write into file;
      default: 'all' (all directional gridlines).
    """
    print('[info] writing gridlines into {} ...'.format(file_path))
    with open(file_path, 'w') as outfile:
      if direction == 'all':
        _, nb_cells_directions = self.get_number_cells()
        outfile.write('\t'.join(str(nb) for nb in nb_cells_directions) + '\n')
    with open(file_path, 'ab') as outfile:
      if direction in ['x', 'all']:
        numpy.savetxt(outfile,
                      self.gridlines[0].get_vertices(precision=precision))
      if direction in ['y', 'all']:
        numpy.savetxt(outfile,
                      self.gridlines[1].get_vertices(precision=precision))
      if len(self.gridlines) == 3 and direction in ['z', 'all']:
        numpy.savetxt(outfile,
                      self.gridlines[2].get_vertices(precision=precision))

  def read(self, file_path):
    """
    Reads the coordinates from the file.

    Parameters
    ----------
    file_name: string
      Name of file containing grid-node stations along each direction.
    """
    print('[info] reading vertices from {} ...'.format(file_path))
    with open(file_path, 'r') as infile:
      nb_divisions = numpy.array([int(n)
                                  for n in infile.readline().strip().split()])
      vertices = numpy.loadtxt(infile, dtype=numpy.float64)
    vertices = numpy.array(numpy.split(vertices,
                                       numpy.cumsum(nb_divisions[:-1] + 1)))
    labels = ['x', 'y', 'z']
    for index, vertices_gridline in enumerate(vertices):
      self.gridlines.append(GridLine(vertices=vertices_gridline,
                                     label=labels[index]))

  def read_yaml_file(self, file_path):
    """
    Parses the YAML file.

    Parameters
    ----------
    file_path: string
      Path the YAML file.

    Returns
    -------
    data: list of dictionaries
      Parsed YAML information.
    """
    print('[info] reading grid parameters from {} ...'.format(file_path))
    with open(file_path, 'r') as infile:
      return yaml.load(infile)

  def write_yaml_file(self, file_path):
    """
    Writes a YAML readable file with information about the Cartesian structured
    mesh.

    Parameters
    ----------
    file_path: string
      Path of the file to write.
    """
    print('[info] writing grid parameters into {} ...'.format(file_path))
    data = []
    for gridline in self.gridlines:
      data.append(gridline.generate_yaml_info())
    nb_cells, nb_cells_directions = self.get_number_cells()
    with open(file_path, 'w') as outfile:
      outfile.write('# {}\n'.format(os.path.basename(file_path)))
      n_cells_string = 'x'.join(str(nb) for nb in nb_cells_directions)
      outfile.write('# {} = {}\n\n'.format(n_cells_string, nb_cells))
      outfile.write(yaml.dump(data, default_flow_style=False))
