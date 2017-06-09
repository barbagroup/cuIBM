"""
Collection of classes and function to parse a given PETSc log file.
"""

import re
import os
import collections
import numpy
from matplotlib import pyplot


class Run(object):
  """
  A PetIBM run.
  """

  def __init__(self, directory=os.getcwd(), label='', logpath=None):
    """
    Sets the directory, the label, and the path of the log file.

    Parameters
    ----------
    directory: string, optional
      Directory of the PetIBM run;
      default: current working directory.
    label: string, optional
      Label of the run (to use in the legend of the figures);
      default: ''.
    logpath: string, optional
      Path of the PETSc log file;
      default: None.
    """
    self.directory = directory
    self.label = label
    self.log = Log(filepath=logpath)

  def get_walltime(self):
    """
    Gets the wall-time (in seconds) of the run.

    Returns
    -------
    walltime: float
      The wall-time.
    """
    self.walltime = self.log.get_walltime()
    return self.walltime

  def get_events(self):
    """
    Gets the information of the different PETSc events.

    Returns
    -------
    events: dictionary
      The PETSc events.
    """
    self.events = self.log.get_events()
    return self.events

  def get_resident_set_size(self, unit='GB'):
    """
    Gets the resident set size (memory used) by the application.

    Parameters
    ----------
    unit: string, optional
      Unit to use;
      choices: `KB`, `MB`, `GB`.
      default: `GB`.

    Returns
    -------
    res: float
      The resident set size.
    """
    self.res = self.log.get_resident_set_size(unit=unit)
    return self.res


class Log(object):
  """
  Contains information parsed from a PETSc log file.
  """

  def __init__(self, filepath=None):
    """
    Sets the path of the PETSc log file.

    Parameters
    ----------
    filepath: string, optional
      Path of the PETSc log file;
      default: None.
    """
    self.path = filepath

  def get_walltime(self):
    """
    Parses the wall-time in seconds.

    Returns
    -------
    walltime: float
      The wall-time.
    """
    look_for = 'Time '
    with open(self.path, 'r') as infile:
      for line in infile.readlines():
        if re.search(look_for, line):
          return float(line.split()[2])
    return None

  def get_resident_set_size(self, unit='GB'):
    """
    Parses the resident set size.

    Parameters
    ----------
    unit: string, optional
      Unit to use;
      choices: 'KB', 'MB', 'GB';
      default: 'GB'.

    Returns
    -------
    res: float
      The resident set size.
    """
    units = {'B': 0, 'KB': 1, 'MB': 2, 'GB': 3}
    look_for = 'process memory'
    with open(self.path, 'r') as infile:
      for line in infile.readlines():
        if re.search(look_for, line):
          return float(line.split()[7]) / 1024**units[unit]
    return None

  def get_events(self):
    """
    Parses information about all PETSc events.

    Returns
    -------
    events: dictionary
      The information about all PETSc events.
    """
    events = collections.OrderedDict()
    with open(self.path, 'r') as infile:
      lines = infile.readlines()
    look_for = 'Summary of Stages'
    for index, line in enumerate(lines):
      if re.search(look_for, line):
        index += 2
        while not re.match(r'^\s*$', lines[index]):
          info = re.split(':', re.sub(' +', ' ', lines[index]))
          label = info[1].strip()
          events[label] = {'index': int(info[0]),
                           'label': info[1].strip(),
                           'walltime': float(info[2].split()[0]),
                           'percent': float(info[2].split()[1][:-1]),
                           'flops': float(info[2].split()[2])}
          index += 1
        break
    return events


def plot_breakdown_percents(runs, event_labels=[],
                            title=None, colors=None):
    """
    Plots a bar chart with the percent of the total wall-time of all events for
    multiple runs.

    Parameters
    ----------
    runs: Run object or list of Run objects
      The list of runs to display on the figure.
    event_labels: string or list of strings, optional
      Names of the events to display on the figure;
      default: [].
    title: string, optional
      Title of the figure;
      default: None.
    colors: iterator, optional
      Colors to use;
      default: None.

    Returns
    -------
    fig: Matplotlib Figure object
      The figure.
    ax: Matplotlib Axes object
      Single or array of axes.
    """
    if not isinstance(runs, (list, tuple)):
      runs = [runs]
    if not isinstance(event_labels, (list, tuple)):
      event_labels = [event_labels]
    fig, ax = pyplot.subplots(figsize=(8.0, 6.0))
    ax.yaxis.grid(zorder=0)
    ax.set_ylabel('% of wall-time', fontsize=16)
    indices = numpy.arange(len(runs))
    bar_width = 0.5
    bar_offsets = numpy.zeros(len(runs))
    for label in event_labels:
      if colors:
        color = next(colors)
      else:
        color = next(ax._get_lines.prop_cycler)['color']
      percents = []
      for run in runs:
        if label in run.events.keys():
          percents.append(run.events[label]['percent'])
        else:
          percents.append(0.0)
      ax.bar(indices, percents, bar_width,
             label=label,
             bottom=bar_offsets,
             color=color,
             linewidth=0,
             zorder=0)
      bar_offsets += percents
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)
    ax.set_xticks(indices + 0.25 * bar_width)
    ax.set_xticklabels([run.label for run in runs], rotation=0, fontsize=16)
    ax.set_yticks([0.0, 25.0, 50.0, 75.0, 100.0],
                  ('0', '25', '50', '75', '100'))
    ax.set_xlim(indices[0] - 0.5, indices[-1] + 1.0)
    ax.set_ylim(0.0, 100.0)
    if title:
      fig.set_title(title)
    return fig, ax


def plot_breakdown_walltimes(runs, event_labels=[],
                             title=None, colors=None):
    """
    Plots a bar chart with the wall-time of all events for
    multiple runs.

    Parameters
    ----------
    runs: Run object or list of Run objects
      The list of runs to display on the figure.
    event_labels: string or list of strings, optional
      Names of the events to display on the figure;
      default: [].
    title: string, optional
      Title of the figure;
      default: None.
    colors: iterator, optional
      Colors to use;
      default: None.

    Returns
    -------
    fig: Matplotlib Figure object
      The figure.
    ax: Matplotlib Axes object
      Single or array of axes.
    """
    if not isinstance(runs, (list, tuple)):
      runs = [runs]
    if not isinstance(event_labels, (list, tuple)):
      event_labels = [event_labels]
    fig, ax = pyplot.subplots(figsize=(8.0, 6.0))
    ax.yaxis.grid(zorder=0)
    ax.set_ylabel('wall-time (s)', fontsize=16)
    indices = numpy.arange(len(runs))
    bar_width = 0.5
    bar_offsets = numpy.zeros(len(runs))
    for label in event_labels:
      if colors:
        color = next(colors)
      else:
        color = next(ax._get_lines.prop_cycler)['color']
      walltimes = []
      for run in runs:
        if label in run.events.keys():
          walltimes.append(run.events[label]['walltime'])
        else:
          walltimes.append(0.0)
      ax.bar(indices, walltimes, bar_width,
             label=label,
             bottom=bar_offsets,
             color=color,
             linewidth=0,
             zorder=0)
      bar_offsets += walltimes
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)
    ax.set_xticks(indices + 0.25 * bar_width)
    ax.set_xticklabels([run.label for run in runs], rotation=0, fontsize=16)
    ax.set_xlim(indices[0] - 0.5, indices[-1] + 1.0)
    if title:
      fig.set_title(title)
    return fig, ax
