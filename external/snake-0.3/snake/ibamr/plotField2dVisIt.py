"""
Calls VisIt in batch mode to generate .png files of the 2D field contour.
cli: visit -nowin -cli -s plotField2dVisIt.py <arguments>
"""

import os
import sys
import math
import argparse

sys.path.append(os.environ['SNAKE'])
from snake import miscellaneous


def parse_command_line():
  """
  Parses the command-line.
  """
  print('[info] parsing the command-line ...'),
  # create the parser
  formatter_class = argparse.ArgumentDefaultsHelpFormatter
  parser = argparse.ArgumentParser(description='Plots the vorticity field '
                                               'with VisIt',
                                   formatter_class=formatter_class)
  # fill the parser with arguments
  parser.add_argument('--directory',
                      dest='directory',
                      type=str,
                      default=os.getcwd(),
                      help='directory of the IBAMR simulation')
  parser.add_argument('--solution-folder',
                      dest='solution_folder',
                      type=str,
                      default='numericalSolution',
                      help='name of folder containing the solution in time')
  parser.add_argument('--body',
                      dest='body',
                      type=str,
                      default=None,
                      help='name of the body file (without the .vertex '
                           'extension)')
  parser.add_argument('--field',
                      dest='field_name',
                      type=str,
                      choices=['vorticity', 'pressure',
                               'x-velocity', 'x-velocity',
                               'velocity-magnitude'],
                      help='name of the field to plot')
  parser.add_argument('--range',
                      dest='field_range',
                      type=float, nargs=2,
                      default=(-1.0, 1.0),
                      metavar=('min', 'max'),
                      help='Range of the field to plot')
  parser.add_argument('--states',
                      dest='states',
                      type=int, nargs=3,
                      default=(0, 2**10000, 1),
                      metavar=('min', 'max', 'increment'),
                      help='steps to plot')
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


def check_version():
  """
  Check the VisIt version and prints warning if the version has not been
  tested.
  """
  script_version = '2.8.2'
  tested_versions = ['2.8.2', '2.10.2', '2.12.1']
  current_version = Version()
  print('VisIt version: {}\n'.format(Version()))
  if current_version not in tested_versions:
    print('[warning] You are using VisIt-{}'.format(current_version))
    print('[warning] This script was created with '
          'VisIt-{}.'.format(script_version))
    print('[warning] This script was tested with versions: '
          '{}.'.format(tested_versions))
    print('[warning] It may not work as expected')


def plot_field_contours(field_name, field_range,
                        body=None,
                        directory=os.getcwd(),
                        solution_folder='numericalSolution',
                        states=(0, 2**10000, 1),
                        view=(-2.0, -2.0, 2.0, 2.0),
                        width=800):
  """
  Plots the contour of a given field using VisIt.

  Parameters
  ----------
  field_name: string
    Name of field to plot;
    choices: vorticity, pressure, velocity-magnitude, x-velocity, y-velocity.
  field_range: 2-tuple of floats
    Range of the field to plot (min, max).
  body: string, optional
    Name of the immersed body;
    default: None.
  directory: string, optional
    Directory of the IBAMR simulation;
    default: current directory.
  solution_folder: string, optional
    Relative path of the folder containing the numerical solution;
    default: 'numericalSolution'.
  states: 3-tuple of integers, optional
    Limits of index of the states to plot followed by the increment;
    default: (0, 20000, 1).
  view: 4-tuple of floats, optional
    Bottom-left and top-right coordinates of the view to display;
    default: (-2.0, -2.0, 2.0, 2.0).
  width: integer, optional
    Width (in pixels) of the figure;
    default: 800.
  """
  info = {}
  info['vorticity'] = {'variable': 'Omega',
                       'color-table': 'RdBu',
                       'invert-color-table': 1}
  info['pressure'] = {'variable': 'P',
                      'color-table': 'hot',
                      'invert-color-table': 0}
  info['velocity-magnitude'] = {'variable': 'U_magnitude',
                                'color-table': 'RdBu',
                                'invert-color-table': 1}
  info['x-velocity'] = {'variable': 'U_x',
                        'color-table': 'RdBu',
                        'invert-color-table': 1}
  info['y-velocity'] = {'variable': 'U_y',
                        'color-table': 'RdBu',
                        'invert-color-table': 1}
  # define dimensions of domain to plot
  height = int(math.floor(width * (view[3] - view[1]) / (view[2] - view[0])))
  # create images directory
  view_string = '{:.2f}_{:.2f}_{:.2f}_{:.2f}'.format(*view)
  images_directory = os.path.join(directory,
                                  'images',
                                  '_'.join([field_name, view_string]))
  if not os.path.isdir(images_directory):
    print('[info] creating images directory {} ...'.format(images_directory))
    os.makedirs(images_directory)

  ShowAllWindows()

  # display body
  if body:
    OpenDatabase(GetLocalHostName() + ':' + os.path.join(directory,
                                                         solution_folder,
                                                         'lag_data.visit'), 0)
    AddPlot('Mesh', body + '_vertices', 1, 1)
    DrawPlots()
    MeshAtts = MeshAttributes()
    MeshAtts.legendFlag = 0
    MeshAtts.lineStyle = MeshAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    MeshAtts.lineWidth = 0
    MeshAtts.meshColor = (0, 0, 0, 255)
    try:
      MeshAtts.outlineOnlyFlag = 0
      MeshAtts.errorTolerance = 0.01
    except:
      pass
    MeshAtts.meshColorSource = MeshAtts.Foreground  # Foreground, MeshCustom
    MeshAtts.opaqueColorSource = MeshAtts.Background  # Background, OpaqueCustom
    MeshAtts.opaqueMode = MeshAtts.Auto  # Auto, On, Off
    MeshAtts.pointSize = 0.05
    MeshAtts.opaqueColor = (255, 255, 255, 255)
    # MeshAtts.smoothingLevel = MeshAtts.None  # None, Fast, High
    MeshAtts.pointSizeVarEnabled = 0
    MeshAtts.pointSizeVar = 'default'
    MeshAtts.pointType = MeshAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
    MeshAtts.showInternal = 0
    MeshAtts.pointSizePixels = 2
    MeshAtts.opacity = 1
    SetPlotOptions(MeshAtts)

  # display vorticity field
  OpenDatabase(GetLocalHostName() + ':' + os.path.join(directory,
                                                       solution_folder,
                                                       'dumps.visit'), 0)
  HideActivePlots()
  AddPlot('Pseudocolor', info[field_name]['variable'], 1, 1)
  DrawPlots()
  PseudocolorAtts = PseudocolorAttributes()
  PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
  PseudocolorAtts.skewFactor = 1
  PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
  PseudocolorAtts.minFlag = 1
  PseudocolorAtts.min = field_range[0]
  PseudocolorAtts.maxFlag = 1
  PseudocolorAtts.max = field_range[1]
  PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
  PseudocolorAtts.colorTableName = info[field_name]['color-table']
  PseudocolorAtts.invertColorTable = info[field_name]['invert-color-table']
  PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
  PseudocolorAtts.opacityVariable = ''
  PseudocolorAtts.opacity = 1
  PseudocolorAtts.opacityVarMin = 0
  PseudocolorAtts.opacityVarMax = 1
  PseudocolorAtts.opacityVarMinFlag = 0
  PseudocolorAtts.opacityVarMaxFlag = 0
  PseudocolorAtts.pointSize = 0.05
  PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
  PseudocolorAtts.pointSizeVarEnabled = 0
  PseudocolorAtts.pointSizeVar = 'default'
  PseudocolorAtts.pointSizePixels = 2
  PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
  PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
  PseudocolorAtts.lineWidth = 0
  if Version() in ['2.8.2', '2.10.2']:
    PseudocolorAtts.tubeDisplayDensity = 10
  elif Version() in ['2.12.1']:
    PseudocolorAtts.tubeResolution = 10
  else:
    PseudocolorAtts.tubeDisplayDensity = 10
  PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
  PseudocolorAtts.tubeRadiusAbsolute = 0.125
  PseudocolorAtts.tubeRadiusBBox = 0.005
  if Version() in ['2.8.2', '2.10.2']:
    PseudocolorAtts.varyTubeRadius = 0
    PseudocolorAtts.varyTubeRadiusVariable = ''
    PseudocolorAtts.varyTubeRadiusFactor = 10
    # PseudocolorAtts.endPointType = PseudocolorAtts.None
    PseudocolorAtts.endPointStyle = PseudocolorAtts.Spheres
  elif Version() in ['2.12.1']:
    PseudocolorAtts.tubeRadiusVarEnabled = 0
    PseudocolorAtts.tubeRadiusVar = ''
    PseudocolorAtts.tubeRadiusVarRatio = 10
    # PseudocolorAtts.tailStyle = PseudocolorAtts.None
    # PseudocolorAtts.headStyle = PseudocolorAtts.None
  else:
    PseudocolorAtts.varyTubeRadius = 0
    PseudocolorAtts.varyTubeRadiusVariable = ''
    PseudocolorAtts.varyTubeRadiusFactor = 10
    # PseudocolorAtts.endPointType = PseudocolorAtts.None
    PseudocolorAtts.endPointStyle = PseudocolorAtts.Spheres
  PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox
  PseudocolorAtts.endPointRadiusAbsolute = 1
  PseudocolorAtts.endPointRadiusBBox = 0.005
  PseudocolorAtts.endPointRatio = 2
  PseudocolorAtts.renderSurfaces = 1
  PseudocolorAtts.renderWireframe = 0
  PseudocolorAtts.renderPoints = 0
  PseudocolorAtts.smoothingLevel = 0
  PseudocolorAtts.legendFlag = 1
  PseudocolorAtts.lightingFlag = 1
  SetPlotOptions(PseudocolorAtts)
  # colorbar of pseudocolor plot
  legend = GetAnnotationObject(GetPlotList().GetPlots(2).plotName)
  legend.xScale = 1.5
  legend.yScale = 0.5
  legend.numberFormat = '%# -9.2g'
  legend.orientation = legend.HorizontalBottom
  legend.managePosition = 0
  legend.position = (0.10, 0.10)
  legend.fontFamily = legend.Courier
  legend.fontBold = 1
  legend.fontHeight = 0.1
  legend.drawMinMax = 0
  legend.drawTitle = 0
  print('[info] legend settings:')
  print(legend)

  # set up view
  View2DAtts = View2DAttributes()
  View2DAtts.windowCoords = (view[0], view[2], view[1], view[3])
  View2DAtts.viewportCoords = (0, 1, 0, 1)
  View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
  View2DAtts.fullFrameAutoThreshold = 100
  View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
  View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
  View2DAtts.windowValid = 1
  print('[info] view settings:')
  print(View2DAtts)
  SetView2D(View2DAtts)

  # Logging for SetAnnotationObjectOptions is not implemented yet.
  AnnotationAtts = AnnotationAttributes()
  AnnotationAtts.axes2D.visible = 1
  AnnotationAtts.axes2D.autoSetTicks = 1
  AnnotationAtts.axes2D.autoSetScaling = 1
  AnnotationAtts.axes2D.lineWidth = 0
  AnnotationAtts.axes2D.tickLocation = AnnotationAtts.axes2D.Inside  # Inside, Outside, Both
  AnnotationAtts.axes2D.tickAxes = AnnotationAtts.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
  # x-axis
  AnnotationAtts.axes2D.xAxis.title.visible = 0 # hide x-axis title
  AnnotationAtts.axes2D.xAxis.label.visible = 0 # hide x-axis label
  AnnotationAtts.axes2D.xAxis.tickMarks.visible = 0 # hide x-axis tick marks
  AnnotationAtts.axes2D.xAxis.grid = 0 # no grid
  # y-axis
  AnnotationAtts.axes2D.yAxis.title.visible = 0 # hide y-axis title
  AnnotationAtts.axes2D.yAxis.label.visible = 0 # hide y-axis label
  AnnotationAtts.axes2D.yAxis.tickMarks.visible = 0 # hide y-axis tick marks
  AnnotationAtts.axes2D.yAxis.grid = 0 # no grid
  AnnotationAtts.userInfoFlag = 0 # hide text with user's name 
  # settings for legend
  AnnotationAtts.databaseInfoFlag = 0
  AnnotationAtts.timeInfoFlag = 0
  AnnotationAtts.legendInfoFlag = 1
  AnnotationAtts.backgroundColor = (255, 255, 255, 255)
  AnnotationAtts.foregroundColor = (0, 0, 0, 255)
  AnnotationAtts.gradientBackgroundStyle = AnnotationAtts.Radial  # TopToBottom, BottomToTop, LeftToRight, RightToLeft, Radial
  AnnotationAtts.gradientColor1 = (0, 0, 255, 255)
  AnnotationAtts.gradientColor2 = (0, 0, 0, 255)
  AnnotationAtts.backgroundMode = AnnotationAtts.Solid  # Solid, Gradient, Image, ImageSphere
  AnnotationAtts.backgroundImage = ''
  AnnotationAtts.imageRepeatX = 1
  AnnotationAtts.imageRepeatY = 1
  AnnotationAtts.axesArray.visible = 0
  SetAnnotationAttributes(AnnotationAtts)
  print('[info] annotation settings:')
  print(AnnotationAtts)
  SetActiveWindow(1)

  # create time-annotation
  time_annotation = CreateAnnotationObject('Text2D')
  time_annotation.position = (0.05, 0.90)
  time_annotation.fontFamily = 1
  time_annotation.fontBold = 0
  time_annotation.height = 0.05
  print('[info] time-annotation:')
  print(time_annotation)

  # check number of states available
  if states[1] > TimeSliderGetNStates():
    print('[warning] maximum number of states available is '
          '{}'.format(TimeSliderGetNStates()))
    print('[warning] setting new final state ...')
    states[1] = TimeSliderGetNStates()

  # loop over saved time-steps
  for state in xrange(args.states[0], args.states[1], args.states[2]):
    SetTimeSliderState(state)
    time = float(Query('Time')[:-1].split()[-1])
    print('\n[state {}] time: {} - creating and saving the field ...'
          ''.format(state, time))
    time_annotation.text = 'Time: {0:.3f}'.format(time)

    RenderingAtts = RenderingAttributes()
    RenderingAtts.antialiasing = 0
    RenderingAtts.multiresolutionMode = 0
    RenderingAtts.multiresolutionCellSize = 0.002
    RenderingAtts.geometryRepresentation = RenderingAtts.Surfaces  # Surfaces, Wireframe, Points
    RenderingAtts.displayListMode = RenderingAtts.Auto  # Never, Always, Auto
    RenderingAtts.stereoRendering = 0
    RenderingAtts.stereoType = RenderingAtts.CrystalEyes  # RedBlue, Interlaced, CrystalEyes, RedGreen
    RenderingAtts.notifyForEachRender = 0
    RenderingAtts.scalableActivationMode = RenderingAtts.Auto  # Never, Always, Auto
    RenderingAtts.scalableAutoThreshold = 2000000
    RenderingAtts.specularFlag = 0
    RenderingAtts.specularCoeff = 0.6
    RenderingAtts.specularPower = 10
    RenderingAtts.specularColor = (255, 255, 255, 255)
    RenderingAtts.doShadowing = 0
    RenderingAtts.shadowStrength = 0.5
    RenderingAtts.doDepthCueing = 0
    RenderingAtts.depthCueingAutomatic = 1
    RenderingAtts.startCuePoint = (-10, 0, 0)
    RenderingAtts.endCuePoint = (10, 0, 0)
    RenderingAtts.compressionActivationMode = RenderingAtts.Never  # Never, Always, Auto
    RenderingAtts.colorTexturingFlag = 1
    RenderingAtts.compactDomainsActivationMode = RenderingAtts.Never  # Never, Always, Auto
    RenderingAtts.compactDomainsAutoThreshold = 256
    SetRenderingAttributes(RenderingAtts)

    SaveWindowAtts = SaveWindowAttributes()
    SaveWindowAtts.outputToCurrentDirectory = 0
    SaveWindowAtts.outputDirectory = images_directory
    SaveWindowAtts.fileName = '{}{:0>7}'.format(field_name, state)
    SaveWindowAtts.family = 0
    SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
    SaveWindowAtts.width = width
    SaveWindowAtts.height = height
    SaveWindowAtts.screenCapture = 0
    SaveWindowAtts.saveTiled = 0
    SaveWindowAtts.quality = 100
    SaveWindowAtts.progressive = 0
    SaveWindowAtts.binary = 0
    SaveWindowAtts.stereo = 0
    SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
    SaveWindowAtts.forceMerge = 0
    SaveWindowAtts.resConstraint = SaveWindowAtts.NoConstraint  # NoConstraint, EqualWidthHeight, ScreenProportions
    SaveWindowAtts.advancedMultiWindowSave = 0
    SetSaveWindowAttributes(SaveWindowAtts)

    SaveWindow()


def main(args):
  check_version()
  plot_field_contours(args.field_name, args.field_range,
                      directory=args.directory,
                      body=args.body,
                      solution_folder=args.solution_folder,
                      states=args.states,
                      view=args.view,
                      width=args.width)
  os.remove('visitlog.py')


if __name__ == '__main__':
  args = parse_command_line()
  main(args)
  sys.exit(0)
