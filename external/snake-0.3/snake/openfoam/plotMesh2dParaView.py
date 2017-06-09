"""
Plots the 2D grid using ParaView.
cli: pvbatch plotMesh2dParaView.py <arguments>
"""

import os
import sys
import argparse

import paraview
from paraview.simple import *

sys.path.append(os.environ['SNAKE'])
from snake import miscellaneous


def parse_command_line():
  """Parses the command-line."""
  print('[info] parsing command-line...'),
  # create the parser
  formatter_class = argparse.ArgumentDefaultsHelpFormatter
  parser = argparse.ArgumentParser(description='Plots the mesh with ParaFOAM',
                                   formatter_class=formatter_class)
  # fill the parser with arguments
  parser.add_argument('--directory',
                      dest='directory',
                      type=str,
                      default=os.getcwd(),
                      help='directory of the OpenFOAM simulation')
  parser.add_argument('--view',
                      dest='view',
                      type=float,
                      nargs=4,
                      default=(-2.0, -2.0, 2.0, 2.0),
                      metavar=('x-bl', 'y-bl', 'x-tr', 'y-tr'),
                      help='bottom-left coordinates followed by top-right '
                           'coordinates of the view to display')
  parser.add_argument('--width',
                      dest='width',
                      type=int,
                      default=800,
                      help='figure width in pixels')
  # parse given options file
  parser.add_argument('--options',
                      type=open,
                      action=miscellaneous.ReadOptionsFromFile,
                      help='path of the file with options to parse')
  print('done')
  return parser.parse_args()


def plot_mesh(directory=os.getcwd(),
              view=(-2.0, -2.0, 2.0, 2.0),
              width=800):
  print('Paraview: \n{}\n'.format(paraview.__path__))
  name = os.path.basename(os.path.normpath(directory))
  openfoam_file_name = '{}.OpenFOAM'.format(name)
  reader = PV4FoamReader(FileName=os.path.join(directory, openfoam_file_name))
  print('[info] plotting the mesh ...')
  reader.MeshParts = ['front - patch']
  # set view
  create_render_view(view=view, width=width)
  view_str = '{:.2f}_{:.2f}_{:.2f}_{:.2f}'.format(*view)
  # create images directory
  images_directory = os.path.join(directory, 'images')
  if not os.path.isdir(images_directory):
    os.makedirs(images_directory)
  print('[info] .png file will be saved in: {}'.format(images_directory))
  # data representation
  data_representation = Show()
  data_representation.Representation = 'Wireframe'
  data_representation.AmbientColor = [0.0, 0.0, 0.0]
  # save
  WriteImage(os.path.join(images_directory, 'mesh_{}.png'.format(view_str)))


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
  render_view.Background = [1.0, 1.0, 1.0]
  Render()
  return render_view


def main(args):
  plot_mesh(directory=args.directory,
            view=args.view,
            width=args.width)


if __name__ == '__main__':
  args = parse_command_line()
  main(args)
