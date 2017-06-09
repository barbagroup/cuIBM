"""
Contains diverse functions and classes.
"""

import os
import sys
import re
import argparse

try:
  import numpy
except:
  pass

try:
  from IPython.display import HTML, Image, display
  import ipywidgets
except ImportError:
  print('[warning] could not import ipywidgets')


class RedirectOutput(object):
  """
  Context manager to redirect `stdout` to a given file.
  """

  def __init__(self, stdout=''):
    self.stdout = stdout

  def __enter__(self):
    self.sys_stdout = sys.stdout

    if self.stdout:
      sys.stdout = open(self.stdout, 'w')

  def __exit__(self, exc_type, exc_value, traceback):
    sys.stdout = self.sys_stdout


class ReadOptionsFromFile(argparse.Action):
  """
  Container to read parameters from file.
  """

  def __call__(self, parser, namespace, values, option_string=None):
    """
    Fills the name-space with parameters read in file.
    """
    with values as infile:
      lines = [element for line in infile.readlines()
               for element in line.strip().split()
               if not line.startswith('#')]
      lines = [os.path.expandvars(line) if '$' in line else line
               for line in lines[:]]
    parser.parse_args(lines, namespace)


def display_image(figure):
  """
  Display figure into the Jupyter-Notebook.
  """
  if not os.path.isfile(figure):
    print('Image not available')
  display(Image(figure))


def get_images(directory, prefix=None, steps=None):
  """
  Returns the list of image paths of interest.

  Parameters
  ----------
  directory: string
    Directory with images.
  prefix: string, optional
    Prefix shared by all images;default: None.
  steps: list of integers, optional
    List of time-steps to consider;
    default: None.

  Returns
  -------
  images: list of strings
    List containing the absolute path of each image of interest.
  """
  if not os.path.isdir(directory):
    print('[error] {} is not a directory'.format(directory))
    return
  try:
    check = not any(steps)
  except:
    check = not steps
  if check:
    images = sorted([os.path.join(directory, image)
                     for image in os.listdir(directory)])
  else:
    if not prefix:
      prefix = re.match(r"([a-z]+)([0-9]+)",
                        os.listdir(directory)[0],
                        re.I).groups()[0]
    if all(isinstance(step, int) for step in steps):
      images = [os.path.join(directory, '{}{:0>7}.png'.format(prefix, step))
                for step in steps]
    else:
      images = [os.path.join(directory, '{}{:06.2f}.png'.format(prefix, step))
                for step in steps]
  return images


def create_slider(values, description='value'):
  """
  Returns the widget slider.

  Parameters
  ----------
  values: list of floats
    Values to store in the widget.
  description: string, optional
    A description of the widget;
    default: 'value'.
  """
  return ipywidgets.FloatSlider(description=description, value=values[0],
                                min=values[0], max=round(values[-1], 2),
                                step=values[1] - values[0])


def displayer(directories=[os.getcwd()], time=(), openfoam=False):
  """
  Interactively displays multiple images in a Jupyter Notebook with ipywidgets.

  Parameters
  ----------
  directories: list of strings, optional
    Path of folders containing images to display;
    default: current directory.
  time: 3-tuple, optional
    Temporal limits and time-increment to choose which images to display;
    default: empty tuple (all images in the folder).
  openfoam: boolean, optional
    Set 'True' if OpenFOAM simulation;
    default: False (not an OpenFOAM simulation).

  Returns
  -------
  The interactive display.
  """
  # check parameter is a list, if not convert into a list
  try:
    assert isinstance(directories, (list, tuple))
    assert not isinstance(directories, basestring)
  except:
    directories = [directories]
  all_images = []
  if not time:
    for directory in directories:
      all_images.append(get_images(directory))
    slider = create_slider(description='index',
                           values=numpy.arange(len(all_images[0])))
  else:
    times = numpy.arange(time[0], time[1] + time[2] / 2.0, time[2])
    if openfoam:
      steps = times
    else:
      steps = numpy.rint(times / time[3]).astype(int)
    for directory in directories:
      all_images.append(get_images(directory, steps=steps))
    slider = create_slider(description='time', values=times)

  def create_view(tic):
    if not time:
      index = int(round(tic))
    else:
      index = numpy.where(numpy.abs(times - tic) <= 1.0E-06)[0][0]
    for images in all_images:
      display(Image(filename=images[index]))

  ipywidgets.interact(create_view, tic=slider)


def load_style(file_path):
  """
  Loads the style for the Jupyter-Notebook.

  Parameters
  ----------
  file_path: string
    Path of the .css file.
  """
  if not file_path:
    file_path = os.path.join(os.environ['SNAKE'],
                             'snake',
                             'styles',
                             'mesnardo.css')
  style = HTML(open(file_path, 'r').read())
  display(style)
