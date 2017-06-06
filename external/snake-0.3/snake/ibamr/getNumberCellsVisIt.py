"""
Calls VisIt in batch mode to grab the number of cells at each requested state
and computes mean.
cli: visit -nowin -cli -s getNumberCellsVisIt.py <arguments>
"""

import os
import sys
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
  parser = argparse.ArgumentParser(description='Counts the number of cells in '
                                               'IBAMR mesh',
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
  parser.add_argument('--states',
                      dest='states',
                      type=int,
                      nargs=3,
                      default=(0, 20000, 1),
                      metavar=('start', 'end', 'increment'),
                      help='states to consider')
  parser.add_argument('--time-limits',
                      dest='time_limits',
                      type=float, nargs=2,
                      default=(0.0, float('inf')),
                      metavar=('time-start', 'time-end'),
                      help='temporal limits within which the mean number of '
                           'cells is calculated')
  # parse given options file
  parser.add_argument('--options',
                      type=open,
                      action=miscellaneous.ReadOptionsFromFile,
                      help='path of the file with options to parse')
  print('done')
  return parser.parse_args()


def check_version():
  """
  Checks the VisIt version and prints warning if the version has not been
  tested.
  """
  script_version = '2.8.2'
  tested_versions = ['2.8.2', '2.10.2']
  current_version = Version()
  if current_version not in tested_versions:
    print('[warning] You are using VisIt-{}'.format(current_version))
    print('[warning] This script was created with '
          'VisIt-{}.'.format(script_version))
    print('[warning] This script was tested with versions: '
          '{}.'.format(tested_versions))
    print('[warning] It may not work as expected')


def get_number_cells(directory=os.getcwd(),
                     solution_folder='numericalSolution',
                     states=(0, 20000, 1)):
  """
  Gets the number of a cells in the AMR mesh at requested states.

  Parameters
  ----------
  directory: string, optional
    Directory of the IBAMR simulation;
    default: <current working directory>.
  solution_folder: string, optional
    Relative path of the folder containing the numerical solution;
    default: 'numericalSolution'.
  states: 3-tuple of integers, optional
    Starting and ending states followed by the increment;
    default: (0, 20000, 1).

  Returns
  -------
  times: list of floats
    Time of each state requested.
  n_cells: list of floats
    Number of cells at each state requested.
  """
  OpenDatabase(GetLocalHostName() + ':' + os.path.join(directory,
                                                       solution_folder), 0)
  AddPlot('Mesh', 'amr_mesh', 1, 1)
  DrawPlots()
  SetQueryFloatFormat('%g')
  times, n_cells = [], []
  # check number of states available
  if states[1] > TimeSliderGetNStates():
    print('[warning] maximum number of states available is '
          '{}'.format(TimeSliderGetNStates()))
    print('[warning] setting new final state ...')
    states[1] = TimeSliderGetNStates()
  for state in xrange(states[0], states[1] + 1, states[2]):
    SetTimeSliderState(state)
    times.append(float(Query('Time')[:-1].split()[-1]))
    n_cells.append(int(Query('NumZones')[:-1].split()[-1]))
    print('[step {}] time: {} ; number of cells: {}'.format(state,
                                                            times[-1],
                                                            n_cells[-1]))
  return times, n_cells


def get_mean(values,
             times=None,
             time_limits=(0.0, float('inf'))):
  """
  Calculates the mean value of a given list.
  If the list is associated with a list of discrete time values, it is possible
  to calculate the mean within a certain time-interval.

  Parameters
  ----------
  values: list of floats
    The list for which the mean is computed.
  times: list of floats, optional
    Time associated with each value;
    default: None.
  time_limits: 2-tuple of floats, optional
    Time-interval to consider for the mean calculation;
    default: (0.0, inf).

  Returns
  -------
  mean: float
    The mean value
  """
  sum_values, count = 0, 0
  if any(times):
    assert len(times) == len(values)
  for index, value in enumerate(values):
    if any(times):
      if time_limits[0] <= times[index] <= time_limits[1]:
        sum_values += value
        count += 1
    else:
      sum_values += value
      count += 1
  return sum_values / count


def main(args):
  check_version()
  n_cells, times = get_number_cells(directory=args.directory,
                                    solution_folder=args.solution_folder,
                                    states=args.states)
  mean = get_mean(n_cells, times,
                  time_limits=args.time_limits)
  print('[info] The AMR grid has on average {} cells '
        'between {} and {} time-units'.format(mean, *args.time_limits))
  os.remove('visitlog.py')


if __name__ == '__main__':
  args = parse_command_line()
  main(args)
  sys.exit(0)
