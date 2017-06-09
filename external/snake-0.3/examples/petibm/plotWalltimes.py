"""
Plots bar charts with the wall-time (in seconds and percent of the total
wall-time) of all PETSc events registered for a PetIBM run.
"""

import itertools
from matplotlib import pyplot
from snake.petibm import logViewReader


run = logViewReader.Run(logpath='run.log')
run.get_walltime()
run.get_events()
run.get_resident_set_size(unit='GB')
print('RES = {:6f} GB'.format(run.res))

pyplot.style.use('seaborn-dark')
colors = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c',
          '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00',
          '#cab2d6', '#6a3d9a']
event_labels = list(run.events.keys())
logViewReader.plot_breakdown_walltimes(run,
                                       event_labels=event_labels,
                                       colors=itertools.cycle(colors))

logViewReader.plot_breakdown_percents(run,
                                      event_labels=event_labels,
                                      colors=itertools.cycle(colors))

pyplot.show()
