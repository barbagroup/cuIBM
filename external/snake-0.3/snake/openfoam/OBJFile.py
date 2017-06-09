"""
Implementation of the classes to generate OBJ files.
"""

import os
import math

import numpy


class Vertex(object):
  """
  Contains info about a vertex.
  """

  def __init__(self, index, x, y, z=0.0):
    """
    Initializes a vertex by its coordinates and its index in the mesh.

    Parameters
    ----------
    index: integer
      Index of the vertex in the mesh.
    x, y: floats
      x- and y-coordinates of the vertex.
    z: float, optional
      z-coordinate of the vertex;
      default: 0.0.
    """
    self.index = index
    self.x, self.y, self.z = x, y, z


class Face(object):
  """
  Contains info about a face.
  """

  def __init__(self, vertex1, vertex2, vertex3):
    """
    Defines a face by its vertices.

    Parameters
    ----------
    vertex1, vertex2, vertex3: Vertex objects
      Vertex objects that compose the face.
    """
    self.vertex1 = vertex1
    self.vertex2 = vertex2
    self.vertex3 = vertex3


class OBJFile(object):
  """
  Contains information about an OBJ file.
  """
  def __init__(self, name):
    """
    Stores name of the OBJ file.

    Parameters
    ----------
    name: string
      Name of the OBJ file.
    """
    self.name = name


class Box2d(OBJFile):
  """
  Contains information about a 2d box OBJ file.
  """
  def __init__(self, name,
               bottom_left=[-1.0, -1.0],
               top_right=[1.0, 1.0],
               n=[10, 10],
               z=0.0):
    """
    Creates the 2d box.

    Parameters
    ----------
    name: string
      Name of the box.
    bottom_left: 2-list of floats, optional
      Bottom-left corner of the box;
      default: [-1.0, -1.0].
    top_right: 2-list of floats, optional
      Top-right corner of the box;
      default: [1.0, 1.0].
    n: 2-list of integers, optional
      Number of points for the discretization in each direction;
      default: [10, 10].
    z: float, optional
      z-location of the 2d box;
      default: 0.0.
    """
    OBJFile.__init__(self, name)
    self.x, self.y = self.create_coordinates(bottom_left, top_right, n)
    self.z = z
    self.vertices = self.create_vertices()
    self.faces = self.create_faces()

  def create_coordinates(self,
                         bottom_left=[-1.0, -1.0],
                         top_right=[1.0, 1.0],
                         n=[10, 10]):
    """
    Creates the coordinates of the box along each direction.

    Parameters
    ----------
    bottom_left: 2-list offloats, optional
      Bottom-left corner of the box;
      default: [-1.0, -1.0].
    top_right: 2-list of floats, optional
      Top-right corner of the box;
      default: [1.0, 1.0].
    n: 2-list of integers, optional
      Number of points for the discretization in each direction;
      default: [10, 10].

    Returns
    -------
    x, y: 1D arrays of floats
      Coordinates along each direction of the box.
    """
    return (numpy.linspace(bottom_left[0], top_right[0], n[0]),
            numpy.linspace(bottom_left[1], top_right[1], n[1]))

  def create_vertices(self):
    """
    Creates the vertices of the OBJ file.

    Returns
    -------
    vertices: 1D array of Vertex objects
      Array containing the vertices.
    """
    nx, ny = self.x.size, self.y.size
    vertices = numpy.empty(nx * ny, dtype=object)
    for j, yj in enumerate(self.y):
      for i, xi in enumerate(self.x):
        vertices[j * nx + i] = Vertex(j * nx + i, xi, yj, self.z)
    return vertices

  def create_faces(self):
    """
    Creates the faces of the OBJ file.

    Returns
    -------
    faces: 1D array of Face objects
      Array containing the faces.
    """
    nx, ny = self.x.size, self.y.size
    lower_faces = numpy.empty((nx - 1) * (ny - 1), dtype=object)
    upper_faces = numpy.empty((nx - 1) * (ny - 1), dtype=object)
    for j in range(ny - 1):
      for i in range(nx - 1):
        lower_faces[j * (nx - 1) + i] = Face(self.vertices[j * nx + i],
                                             self.vertices[j * nx + (i + 1)],
                                             self.vertices[(j + 1) * nx + i])
        upper_faces[j * (nx - 1) + i] = Face(self.vertices[(j + 1) * nx + (i + 1)],
                                             self.vertices[(j + 1) * nx + i],
                                             self.vertices[j * nx + (i + 1)])
    return numpy.insert(upper_faces,
                        numpy.arange(len(lower_faces)), lower_faces)

  def write(self,
            save_directory=os.getcwd()):
    """
    Writes object into a OBJ file.

    Parameters
    ----------
    save_directory: string, optional
      Directory where to save the OBJ file;
      default: <current directory>.
    """
    print('[info] writing OBJ file ...'),
    nx, ny = self.x.size, self.y.size
    header = ('# Wavefront OBJ file\n'
              '# points: {}\n'
              '# faces: {}\n'
              '# zones: 1\n'
              '# regions: 0 {}\n'.format(nx * ny, 2 * (nx - 1) * (ny - 1),
                                         self.name))
    obj_path = '{}/{}.obj'.format(save_directory, self.name)
    with open(obj_path, 'w') as outfile:
      outfile.write(header)
      for vertex in self.vertices:
        outfile.write('v {} {} {}\n'.format(vertex.x, vertex.y, vertex.z))
      outfile.write('g {}\n'.format(self.name))
      for face in self.faces:
        outfile.write('f {} {} {}\n'.format(face.vertex1.index + 1,
                                            face.vertex2.index + 1,
                                            face.vertex3.index + 1))
    print('done')
    print('path: {}'.format(obj_path))


class Body2d(OBJFile):
  """
  Contains information about the body OBJ file.
  """

  def __init__(self, file_path,
               name=None,
               extrusion_limits=[0.0, 1.0]):
    """
    Reads the coordinates of the 2d geometry.

    Parameters
    ----------
    file_path: string
      Path of the coordinates file.
    name: string, optional
      Name of the body;
      default: None.
    extrusion_limits: list of floats, optional
      Limits of the extrusion;
      default: [0.0, 1.0].
    """
    if not name:
      name = os.path.splitext(os.path.basename(file_path))[0]
    OBJFile.__init__(self, name)
    self.x, self.y = self.read_coordinates(file_path)
    self.extrusion_limits = extrusion_limits

  def read_coordinates(self, file_path):
    """
    Reads the x- and y- coordinates from input file.

    Parameters
    ----------
    file_path: string
      Path of the coordinates file.

    Returns
    -------
    x, y: 1d arrays of floats
      x- and y-coordinates in closed loop form.
    """
    print('[info] reading input coordinates from {} ...'.format(file_path)),
    # read the coordinates file
    with open(file_path, 'r') as infile:
      x, y = numpy.loadtxt(infile,
                           dtype=float,
                           delimiter='\t',
                           skiprows=1,
                           unpack=True)
    # append first element to the end of the array if different
    if math.sqrt((x[0] - x[-1])**2 + (y[0] - y[-1])**2) > 1.0E-06:
      x, y = numpy.append(x, x[0]), numpy.append(y, y[0])
    print('done')
    return x, y

  def write(self,
            save_directory=os.getcwd()):
    """
    Writes the coordinates in a .obj format.

    Parameters
    ----------
    save_directory: string, optional
      Directory where to save the .obj file;
      default: <current directory>.
    """
    print('[info] writing OBJ file ...'),
    outfile_path = '{}/{}.obj'.format(save_directory, self.name)
    header = ('# Wavefront OBJ file\n'
              '# points: {}\n'
              '# faces: {}\n'
              '# zones: 1\n'
              '# Regions: 0 {}\n'.format(2 * self.x.size,
                                         2 * self.x.size,
                                         self.name))
    with open(outfile_path, 'w') as outfile:
      outfile.write(header)
      for i in range(self.x.size):
        for z in self.extrusion_limits[::-1]:
          outfile.write('v {} {} {}\n'.format(self.x[i], self.y[i], z))
          # rotated geometry required by Nikos to run their code
          # outfile.write('v {} {} {}\n'.format(self.y[i], z, self.x[i]))
      outfile.write('g {}\n'.format(self.name))
      for i in range(1, self.x.size):
        outfile.write('f {} {} {}\n'.format(2 * i, 2 * i - 1, 2 * i + 1))
        outfile.write('f {} {} {}\n'.format(2 * i + 1, 2 * (i + 1), 2 * i))
      outfile.write('f {} {} {}\n'.format(2 * self.x.size,
                                          2 * self.x.size - 1,
                                          1))
      outfile.write('f {} {} {}\n'.format(1, 2, 2 * self.x.size))
    print('done')
    print('path: {}'.format(outfile_path))
