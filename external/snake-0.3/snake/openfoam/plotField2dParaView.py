"""
Plots the 2D solution from an OpenFOAM simulation using ParaView.
cli: pvbatch plotField2dParaView.py <arguments>
"""

import os
import sys
import argparse

import numpy
import paraview
from paraview.simple import *

sys.path.append(os.environ['SNAKE'])
from snake import miscellaneous


def parse_command_line():
  """
  Parses the command-line.
  """
  print('[info] parsing command-line...'),
  # create the parser
  formatter_class = argparse.ArgumentDefaultsHelpFormatter
  parser = argparse.ArgumentParser(description='Plots the vorticity field '
                                               'with ParaFOAM',
                                   formatter_class=formatter_class)
  # fill the parser with arguments
  parser.add_argument('--directory',
                      dest='directory',
                      type=str,
                      default=os.getcwd(),
                      help='directory of the OpenFOAM simulation')
  parser.add_argument('--field',
                      dest='field_name',
                      type=str,
                      choices=['vorticity', 'pressure',
                               'x-velocity', 'y-velocity'],
                      help='name of the field to plot')
  parser.add_argument('--range',
                      dest='field_range',
                      type=float, nargs=2,
                      default=(-1.0, 1.0),
                      metavar=('min', 'max'),
                      help='range of the field to plot')
  parser.add_argument('--times',
                      dest='times',
                      type=float, nargs=3,
                      default=(0, 0, 0),
                      metavar=('min', 'max', 'increment'),
                      help='times to plot')
  parser.add_argument('--view',
                      dest='view',
                      type=float, nargs=4,
                      default=(-2.0, -2.0, 2.0, 2.0),
                      metavar=('x-bl', 'y-bl', 'x-tr', 'y-tr'),
                      help='bottom-left coordinates followed by top-right '
                           'coordinates of the view to display')
  parser.add_argument('--width',
                      dest='width',
                      type=int,
                      default=800,
                      help='figure width in pixels')
  parser.add_argument('--colormap', dest='colormap_path',
                      type=str,
                      default=None,
                      help='path of the colormap file to use')
  parser.add_argument('--mesh', dest='display_mesh',
                      action='store_true',
                      help='displays the mesh (Surface with Edges)')
  parser.add_argument('--no-scalar-bar', dest='display_scalar_bar',
                      action='store_false',
                      help='does not display the scalar bar')
  parser.add_argument('--no-time-text', dest='display_time_text',
                      action='store_false',
                      help='does not display the time-unit')
  parser.set_defaults(display_scalar_bar=True, display_time_text=True)
  # parse given options file
  parser.add_argument('--options',
                      type=open,
                      action=miscellaneous.ReadOptionsFromFile,
                      help='path of the file with options to parse')
  print('done')
  return parser.parse_args()


def plot_field_contours(field_name,
                        field_range=(-1.0, 1.0),
                        directory=os.getcwd(),
                        view=(-2.0, -2.0, 2.0, 2.0),
                        times=(0, 0, 0),
                        width=800,
                        colormap_path=None,
                        display_scalar_bar=True,
                        display_time_text=True,
                        display_mesh=False):
  """
    Plots the contour of a given field using ParaView.

    Parameters
    ----------
    field_name: string
      Name of field to plot;
      choices: vorticity, pressure, x-velocity, y-velocity.
    field_range: 2-tuple of floats, optional
      Range of the field to plot (min, max);
      default: (-1.0, 1.0).
    directory: string, optional
      Directory where to save the .png files;
      default: <current working directory>.
    view: 4-tuple of floats, optional
      Bottom-left and top-right coordinates of the view to display;
      default: (-2.0, -2.0, 2.0, 2.0).
    times: 3-tuple of floats, optional
      Time-limits followed by the time-increment to consider;
      default: (0, 0, 0).
    width: integer, optional
      Width (in pixels) of the figure;
      default: 800.
    colormap: string, optional
      Name of the Matplotlib colormap to use;
      default: None.
    display_scalar_bar: boolean, optional
      Displays the scalar bar;
      default: True.
    display_time_text: boolean, optional
      Displays the time-unit in the top-left corner;
      default: True.
    display_mesh: boolean, optional
      Displays the mesh (Surface with Edges);
      default: False
    """
  print('Paraview: \n{}\n'.format(paraview.__path__))
  name = os.path.basename(os.path.normpath(directory))
  openfoam_file_name = '{}.OpenFOAM'.format(name)
  reader = PV4FoamReader(FileName=os.path.join(directory, openfoam_file_name))
  print('[info] plotting {} field ...'.format(field_name))
  variable_names = {'vorticity': 'vorticity',
                    'pressure': 'p',
                    'x-velocity': 'U',
                    'y-velocity': 'U'}
  reader.VolumeFields = [variable_names[field_name]]
  reader.MeshParts = ['front - patch']
  # set view
  render_view = create_render_view(view=view, width=width)
  # grab times available
  if all(times) == 0.0:
    times = numpy.array(reader.TimestepValues)
  else:
    start, end, increment = times
    times = numpy.arange(start, end + increment / 2.0, increment)
  # create images directory
  view_str = '{:.2f}_{:.2f}_{:.2f}_{:.2f}'.format(*view)
  images_directory = os.path.join(directory, 'images',
                                  field_name + '_' + view_str)
  if not os.path.isdir(images_directory):
    os.makedirs(images_directory)
  print('[info] .png files will be saved in: {}'.format(images_directory))
  # edit colormap
  PVLookupTable = edit_colormap(field_name, field_range,
                                colormap_path=colormap_path)
  # add a scalar bar
  if display_scalar_bar:
    scalar_bar = add_scalar_bar(field_name, PVLookupTable)
    # update view
    render_view.Representations.append(scalar_bar)
  # show field
  data_representation = Show()
  if display_mesh:
    data_representation.Representation = 'Surface With Edges'
    data_representation.AmbientColor = [0.0, 0.0, 0.0]
    data_representation.EdgeColor = [0.0, 0.0, 0.0]
  data_representation.ColorArrayName = variable_names[field_name]
  data_representation.LookupTable = PVLookupTable
  data_representation.ColorAttributeType = 'CELL_DATA'
  # add text to view
  if display_time_text:
    text = Text()
    data_representation_3 = Show()
    data_representation_3.FontSize = 12
    data_representation_3.TextScaleMode = 2
    data_representation_3.Position = [0.02, 0.9]  # 0.0, 0.0: bottom-left
    data_representation_3.Color = [0.0, 0.0, 0.0]
  # plot and save requested contours
  for time in times:
    print('[info] creating view at {} time-units ...'.format(time))
    render_view.ViewTime = time
    if display_time_text:
      text.Text = 'time = {}'.format(time)
    WriteImage(os.path.join(images_directory,
                            '{}{:06.2f}.png'.format(field_name, time)))


def create_render_view(view=(-2.0, -2.0, 2.0, 2.0), width=800):
  center = [0.5 * (view[0] + view[2]), 0.5 * (view[1] + view[3])]
  h = 1.0 + 0.5 * abs(view[3] - view[1])
  height = int(width * abs(view[3] - view[1]) / abs(view[2] - view[0]))
  render_view = GetRenderView()
  render_view.ViewSize = [width, height]
  Render()  # needed
  render_view.CenterAxesVisibility = 0
  render_view.OrientationAxesVisibility = 0
  render_view.CameraPosition = [center[0], center[1], h]
  render_view.CameraFocalPoint = [center[0], center[1], 0.0]
  render_view.CameraViewUp = [0.0, 1.0, 0.0]
  render_view.CenterOfRotation = [0.0, 0.0, 1.0]
  render_view.CameraViewAngle = 90.0
  render_view.Background = [0.34, 0.34, 0.34]
  Render()
  return render_view


def edit_colormap(field_name, field_range,
                  colormap_path=None):
  mode = {'vorticity': {'variable': 'vorticity',
                        'vectormode': 'Component',
                        'vectorcomponent': 2,
                        'colorspace': 'Diverging'},
          'x-velocity': {'variable': 'U',
                         'vectormode': 'Component',
                         'vectorcomponent': 0,
                         'colorspace': 'Diverging'},
          'y-velocity': {'variable': 'U',
                         'vectormode': 'Component',
                         'vectorcomponent': 1,
                         'colorspace': 'Diverging'},
          'pressure': {'variable': 'p',
                       'vectormode': 'Magnitude',
                       'vectorcomponent': 0,
                       'colorspace': 'HSV'}, }
  min_range, max_range = round(field_range[0], 2), round(field_range[1], 2)
  if colormap_path:
    RGBPoints = []
    with open(colormap_path, 'r') as infile:
      data = infile.readlines()
    increment = abs(max_range - min_range) / (len(data) - 1)
    for i, d in enumerate(data):
      r, g, b = d.strip().split(',')
      RGBPoints.append(min_range + i * increment)
      RGBPoints.append(float(r))
      RGBPoints.append(float(g))
      RGBPoints.append(float(b))
  else:
    RGBPoints = [min_range, 0.0, 0.0, 1.0, max_range, 1.0, 0.0, 0.0]
  return GetLookupTableForArray(mode[field_name]['variable'],
                                mode[field_name]['vectorcomponent'] + 1,
                                RGBPoints=RGBPoints,
                                VectorMode=mode[field_name]['vectormode'],
                                VectorComponent=mode[field_name]['vectorcomponent'],
                                NanColor=[0.0, 0.0, 0.0],
                                ColorSpace=mode[field_name]['colorspace'],
                                ScalarRangeInitialized=1.0,
                                LockScalarRange=1)


def add_scalar_bar(field_name, PVLookupTable):
  return CreateScalarBar(ComponentTitle='',
                         Title='',
                         Enabled=1,
                         LabelFontSize=10,
                         LabelColor=[0.0, 0.0, 0.0],
                         LookupTable=PVLookupTable,
                         TitleFontSize=10,
                         TitleColor=[0.0, 0.0, 0.0],
                         Orientation='Horizontal',
                         Position=[0.04, 0.1],
                         Position2=[0.2, 0.1])


def main(args):
  plot_field_contours(args.field_name,
                      field_range=args.field_range,
                      directory=args.directory,
                      view=args.view,
                      times=args.times,
                      width=args.width,
                      colormap_path=args.colormap_path,
                      display_scalar_bar=args.display_scalar_bar,
                      display_time_text=args.display_time_text,
                      display_mesh=args.display_mesh)


if __name__ == '__main__':
  args = parse_command_line()
  main(args)
