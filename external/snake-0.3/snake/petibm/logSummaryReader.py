"""
Collection of classes and functions to parse a given PETSc log summary.
"""

import os
import re
import collections

import numpy
from matplotlib import pyplot

# load matplotlib style sheet
style_path = os.path.join(os.environ['SNAKE'],
                          'snake',
                          'styles',
                          'mesnardo.mplstyle')
pyplot.style.use(style_path)


class Series(object):
  def __init__(self, directory, nprocs,
               description='no description'):
    """
    Registers the series.

    Parameters
    ----------
    directory: string
      Directory of the series.
    nprocs: list of integers
      Number of processes used for each simulation in the series.
    description: string, optional
      Description of the series;
      default: 'no description'.
    """
    self.directory = directory
    self.nprocs = numpy.array(nprocs)
    self.description = description
    self.runs = self.get_runs()

  def get_runs(self):
    directories = self.get_runs_directory()
    runs = []
    for directory in directories:
      runs.append(Run(directory))
    return runs

  def get_runs_directory(self):
    if self.nprocs.all():
      return ['{}/n{}'.format(self.directory, n) for n in self.nprocs]
    else:
      return ['{}/{}'.format(self.directory, d)
              for d in os.listdir(self.directory) if os.path.isdir(d)]

  def get_wall_times(self, event=None):
    """Gets the wall-time of a given event for each simulation.

    Parameters
    ----------
    event: string, optional
      Name of the event;
      default: None.
    """
    if event:
      return numpy.array([run.log_summary.events[event].time
                          for run in self.runs])
    else:
      return numpy.array([run.log_summary.wall_time for run in self.runs])

  def print_average_solvers_iterations(self, start=0, end=-1):
    """
    Prints the average number of iterations of each solver
    between two time-steps.

    Parameters
    ----------
    start, end: integers, optional
      Starting and ending time-steps;
      default: 0, -1 (last).
    """
    print('Average number of iterations for each solver:')
    for i, run in enumerate(self.runs):
      n_ite_1, n_ite_2 = run.get_average_solvers_iterations(start=start,
                                                            end=end)
      print('{}: {}\t{}'.format(self.nprocs[i], n_ite_1, n_ite_2))

  def plot_wall_time_vs_process_count(self, other_series=[],
                                      event=None, title=None,
                                      save=None, show=False):
    """
    Plots the wall-time versus the process count.

    Parameters
    ----------
    others_series: list of Series objects, optional
      List of other series to plot;
      default: [].
    event: string, optional
      Name of event to plot;
      default: None (wall-time of the simulation).
    title: string, optional
      Title of the plot;
      default: None (no title).
    save: string, optional
      Name of the .png file to be saved;
      default: None (does not save).
    show: boolean, optional
      Displays the figure is set to `True`;
      default: False.
    """
    # set up figure
    fig, ax = pyplot.subplots(figsize=(8, 6))
    ax.yaxis.grid(zorder=0)
    pyplot.xlabel('process count', fontsize=16)
    pyplot.ylabel('wall-time (s)', fontsize=16)
    # plot master series
    times = self.get_wall_times(event=event)
    pyplot.plot(self.nprocs, times,
                label=self.description, marker='o', zorder=10)
    # plot other series
    for series in other_series:
      pyplot.plot(series.nprocs, series.get_wall_times(event=event),
                  label=series.description, marker='o', zorder=10)
    # plot ideal scaling
    pyplot.plot(self.nprocs, 0.5 * times[0] / self.nprocs,
                label='ideal scaling', color='k', linestyle='--', zorder=10)
    # set plot parameters
    pyplot.legend(fontsize=14)
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.xlim(0, self.nprocs[-1])
    if title:
      pyplot.title(title)
    if save:
      pyplot.savefig(save)
    if show:
      pyplot.show()

  def plot_breakdown(self, events_name=['solvePoisson', 'solveVelocity',
                                        'projectionStep', 'RHSPoisson',
                                        'RHSVelocity'],
                     title=None, save=None, show=False):
    """
    Plots the percentage of wall-time spent for each event in each run of a
    series.

    Parameters
    ----------
    events_name: list of strings, optional
      List of name of events to display in the breakdown chart;
      default: ['solvePoisson', 'solveVelocity',
                'projectionStep', 'RHSPoisson', 'RHSVelocity'].
    title: string, optional
      Title of the plot;
      default: None (no title).
    save: string, optional
      Path of the .png to save;
      default: None (i.e. not saving).
    show: boolean, optional
      Displays the figure if set to 'True';
      default: False.
    """
    fig, ax = pyplot.subplots(figsize=(8, 6))
    ax.yaxis.grid(zorder=0)
    pyplot.xlabel('process count', fontsize=16)
    pyplot.ylabel('% of wall time', fontsize=16)
    color_cycle = ax._get_lines.prop_cycler
    index = numpy.arange(len(self.runs))
    bar_width = 0.5
    bar_offset = numpy.zeros(len(self.runs))
    for name in events_name:
      percents = [run.log_summary.events[name].percent for run in self.runs]
      pyplot.bar(index, percents, bar_width,
                 label=name, bottom=bar_offset,
                 color=next(color_cycle)['color'], linewidth=0, zorder=0)
      bar_offset += percents
    pyplot.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)
    pyplot.xticks(index + 0.5 * bar_width, self.nprocs.astype('str'))
    pyplot.yticks([0.0, 25.0, 50.0, 75.0, 100.0],
                  ('0', '25', '50', '75', '100'))
    pyplot.xlim(index[0] - 0.5, index[-1] + 1.0)
    pyplot.ylim(0.0, 100.0)
    if title:
      pyplot.title(title)
    if save:
      pyplot.savefig(save, bbox_inches='tight')
    if show:
      pyplot.show()


class GroupSeries(Series):
  def __init__(self, directories, nprocs,
               descriptions=None, description='no description'):
    """
    Registers of group of Series.

    Parameters
    ----------
    directories: list of strings
      Directory of each series to consider.
    nprocs: list of integers
      Number of processes used for each simulation of a series.
    descriptions: list of strings, optional
      Description of each series;
      default: None.
    description: string, optional
      Description of the group of series;
      default: 'no description'.
    """
    self.nprocs = numpy.array(nprocs)
    self.series_list = self.get_series(directories, descriptions)
    self.description = description
    # self.average_wall_times = self.get_wall_times()

  def get_series(self, directories, descriptions=None):
    """
    Creates each series.

    Parameters
    ----------
    directories: list of strings
      Directory of each series to create.
    descriptions: list of strings, optional
      Description of each series to create;
      default: None.

    Returns
    -------
    series_list: list of Series objects
      The group of series to consider.
    """
    series_list = []
    for i, directory in enumerate(directories):
      series_list.append(Series(directory, self.nprocs,
                                description=(None if not descriptions
                                             else descriptions[i])))
    return series_list

  def get_wall_times(self, event=None):
    """
    Gets the wall-time of an event for each simulation in a series
    averaged among a group of series.

    Parameters
    ----------
    event: string, optional
      Name of the event to consider;
      default: None.

    Returns
    -------
    wall_times: 1d array of floats
      The averaged wall-time for each simulation.
    """
    wall_times = numpy.zeros(len(self.nprocs))
    for series in self.series_list:
      wall_times += series.get_wall_times(event=event)
    return wall_times / len(self.series_list)


class Run(object):
  def __init__(self, directory, description=None):
    """
    Registers a simulation (a run).

    Parameters
    ----------
    directory: string
      Directory of the run.
    description: string, optional
      Description of the run;
      default: None.
    """
    self.directory = directory
    self.description = description
    self.log_summary = LogSummary(self.get_log_summary_path())

  def get_log_summary_path(self):
    """
    Returns the path of the the file with a PETSc summary of the run.
    """
    return ['{}/{}'.format(self.directory, name)
            for name in os.listdir(self.directory) if '.out' in name][0]

  def get_wall_time(self):
    """
    Returns the wall-time of the run.
    """
    return self.log_summary.get_wall_time()

  def get_nprocs(self):
    """
    Returns the number of processes used for the run.
    """
    return self.log_summary.get_nprocs()

  def get_average_solvers_iterations(self, start=0, end=-1):
    """
    Gets the number of iterations performed for each solver
    between two time-steps.

    Parameters
    ----------
    start, end: integers, optional
      Starting and ending time-steps to consider;
      default: 0, -1 (last).

    Returns
    -------
    data: 2-list of floats
      Averaged number of iterations for each solver.
    """
    with open('{}/iterationCounts.txt'.format(self.directory), 'r') as infile:
      data = numpy.loadtxt(infile, dtype=int)
    return data[start:end, 1].mean(), data[start:end, 2].mean()


class LogSummary(object):
  """
  Contains info about what is happening during a PetIBM run.
  """

  def __init__(self, file_path):
    """
    Parses the log file to get the execution time,
    the number of processes involved,
    and info about the various event that occur during the run.

    Parameters
    ----------
    file_path: string
      Path of the the logging file.
    """
    self.file_path = file_path
    self.nprocs = self.get_nprocs()
    self.wall_time = self.get_wall_time()
    self.events = self.get_events()

  def get_nprocs(self):
    """
    Finds the number of processes involved.
    """
    pattern = 'processor'
    with open(self.file_path, 'r') as infile:
      for line in infile.readlines():
        match = re.search(pattern, line)
        if match:
          data = line.split()
          return int([n for i, n in enumerate(data[:-1])
                      if pattern in data[i + 1]][0])

  def get_wall_time(self):
    """
    Finds the execution-time of the simulation.
    """
    pattern = 'Time '
    with open(self.file_path, 'r') as infile:
      for line in infile.readlines():
        match = re.search(pattern, line)
        if match:
          data = line.split()
          return float(data[2])

  def get_events(self):
    """
    Finds the different events of the run and parses info.
    """
    events = collections.OrderedDict()
    pattern = 'Summary of Stages'
    with open(self.file_path, 'r') as infile:
      lines = infile.readlines()
    for index, line in enumerate(lines):
      match = re.search(pattern, line)
      if match:
        index += 2
        while not re.match(r'^\s*$', lines[index]):
          info = re.split(':', re.sub(' +', ' ', lines[index]))
          name = info[1].strip()
          events[name] = Event(info)
          index += 1
        break
    return events

  def print_events_name(self):
    """
    Prints the list of events that occur during the run.
    """
    print('List of events for {}:'.format(self.file_path))
    for name in self.events.keys():
      print('\t{}'.format(name))


class Phase(object):
  """
  Contains info related to a phase of an event.
  """

  def __init__(self, info):
    """
    Parses the info related to a phase.

    Parameters
    ----------
    info: list of strings
        Data related to the phase.
    """
    self.count = int(info[0])
    self.max_time = float(info[2])
    self.ratio_time = float(info[3])
    self.max_flops = float(info[4])
    self.ratio_flops = float(info[5])


class Event(object):
  """
  Contains info related to an event.
  """

  def __init__(self, info):
    """
    Parses general info about the event.

    Parameters
    ----------
    info: list of strings
        Data about the general info of the event.
    """
    self.index = int(info[0])
    self.name = info[1].strip()
    info = info[2].split()
    self.time = float(info[0])
    self.percent = float(info[1][:-1])
    self.flops = float(info[2])

  def parse_phases(self, info):
    """
    Gets info about all phases that occur during the event.

    Parameters
    ----------
    info: list of strings
        Data about all phases that occur during the event.
    """
    self.phases = collections.OrderedDict()
    for data in info:
      name = data[0]
      self.phases[name] = Phase(data[1:])

  def print_phase_names(self):
    """
    Prints the list of phases that occur during the event.
    """
    print('List of phases for event {}:'.format(self.name))
    for name in self.phases.keys():
      print('\t{}'.format(name))


def plot_phases_event(logs,
                      event_name, variable_name,
                      phase_names=[],
                      xlabel='process count', ylabel='time (s)',
                      save=None):
  """
  Plots info about the phases of a given event.

  Parameters
  ----------
  logs: list of LogSummary objects
      List of parsed logging file.
  event_name: string
      Name of the event to plot.
  variable_name: string
      Name of the variable of the phase to consider.
  phase_names: list of strings, optional
      List of phase names to consider;
      default [] (all phases).
  xlabel: string, optional
      x-label;
      default: 'process count'.
  ylabel: string, optional
      y-label;
      default: 'time (s)'.
  save: string, optional
      name of the .png file to be save;
      default: None (does not save).
  """
  fig, ax = pyplot.subplots(figsize=(8, 6))
  ax.yaxis.grid(zorder=0)
  pyplot.xlabel(xlabel, fontsize=16)
  pyplot.ylabel(ylabel, fontsize=16)
  pyplot.title(event_name)
  nprocs = [log.nprocs for log in logs]
  data = collections.OrderedDict()
  index = 0
  if not phase_names:
    phase_names = logs[0].events[event_name].phases.keys()
  for name in phase_names:
    index += 1
    marker = ('o' if index <= 7 else 's')
    data[name] = []
    for log in logs:
      data[name].append(getattr(log.events[event_name].phases[name],
                                variable_name))
    pyplot.plot(nprocs, data[name],
                label=name,
                marker=marker, zorder=10)
  pyplot.xscale('log')
  if 'max' in variable_name:
    pyplot.yscale('log')
  pyplot.xlim(min(nprocs), max(nprocs))
  pyplot.legend(loc='best', fontsize=14,
                bbox_to_anchor=(1.0, 1.0), frameon=False)
  if save:
    pyplot.savefig('./images/{}.png'.format(save))


def plot_phases_event_average(series_list,
                              event_name, variable_name,
                              phase_names=[],
                              xlabel='process count',
                              ylabel='time (s)',
                              save=None):
  """
  Plots info about the phases of a given event.

  Parameters
  ----------
  series_list: list of Series objects
      List of of series,
      each element of a series being a parsed log file.
  event_name: string
      Name of the event to plot.
  variable_name: string
      Name of the variable of the phase to consider.
  phase_names: list of strings, optional
      List of phase names to consider;
      default [] (all phases).
  xlabel: string, optional
      x-label;
      default: 'process count'.
  ylabel: string, optional
      y-label;
      default: 'time (s)'.
  save: string, optional
      name of the .png file to be save;
      default: None (does not save).
  """
  # grab info common to all series
  series_model = series_list[0]
  nprocs = [log.nprocs for log in series_model.logs]
  if not phase_names:
    phase_names = series_model.logs[0].events[event_name].phases.keys()
  # set up fure parameters
  fig, ax = pyplot.subplots(figsize=(8, 6))
  ax.yaxis.grid(zorder=0)
  pyplot.xlabel(xlabel, fontsize=16)
  pyplot.ylabel(ylabel, fontsize=16)
  pyplot.title(event_name)
  # get data to plot
  data = collections.OrderedDict()
  index = 0
  for name in phase_names:
    index += 1
    marker = ('o' if index <= 7 else 's')
    data[name] = [0] * len(series_model.logs)
    # average data over all series in list
    for series in series_list:
      for i, log in enumerate(series.logs):
        data[name][i] += getattr(log.events[event_name].phases[name],
                                 variable_name)
    data[name] = numpy.array(data[name]) / len(series_list)
    # plot data
    pyplot.plot(nprocs, data[name],
                label=name,
                marker=marker, zorder=10)
  pyplot.xscale('log')
  pyplot.xlim(min(nprocs), max(nprocs))
  if 'max' in variable_name:
    pyplot.yscale('log')
  pyplot.legend(loc='best', fontsize=14,
                bbox_to_anchor=(1.0, 1.0), frameon=False)
  if save:
    pyplot.savefig('./images/{}.png'.format(save))
