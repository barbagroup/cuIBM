"""
Get observed order of convergence using the three first grids and the three
last grids.
"""

import os
import numpy

from snake.cuibm.simulation import CuIBMSimulation
from snake.convergence import get_observed_orders


def read_mask(filepath, nx, ny):
  """
  Reads the mask for x- and y- velocities from file.
  The mask can be applied to set the field values to 0.0 inside the body.

  Parameters
  ----------
  filepath: string
    Path the file containing the mask (output of cuIBM).
  nx: integer
    Number of grid divisions in the x-direction.
  ny: integer
    Number of grid divisions in the y-direction.

  Returns
  -------
  mask_u: 2d array of integers
    Mask for the x-velocity.
  mask_v: 2d array of integers
    Mask for the y-velocity.
  """
  with open(filepath, 'r') as infile:
    mask = numpy.loadtxt(infile, dtype=numpy.int32)
  offset = (nx - 1) * ny
  return mask[:offset].reshape(ny, nx - 1), mask[offset:].reshape(ny - 1, nx)


data = {'20x20': {'directory': os.path.join('circ_T00.25_20_0.00050',
                                            'linear'),
                  'folders': ['20', '60', '180', '540']},
        '25x25': {'directory': os.path.join('circ_T00.25_25_0.00050',
                                            'linear'),
                  'folders': ['25', '75', '225', '675']},
        '30x30': {'directory': os.path.join('circ_T00.25_30_0.00050',
                                            'linear'),
                  'folders': ['30', '90', '270', '810']},
        'uquad': {'directory': os.path.join('uquad_T00.25_20_0.00050',
                                            'linear'),
                  'folders': ['20', '60', '180', '540']},
        'vquad': {'directory': os.path.join('vquad_T00.25_20_0.00050',
                                            'linear'),
                  'folders': ['20', '60', '180', '540']}}

time_step = 500
field_names = ['x-velocity', 'y-velocity']

for key, series in data.items():
  series['cases'] = []
  for folder in series['folders']:
    case = CuIBMSimulation(directory=os.path.join(series['directory'], folder),
                           description=folder)
    case.read_grid()
    case.read_fields(field_names, time_step)
    nx, ny = case.grid[0].size - 1, case.grid[1].size - 1
    masks = read_mask(os.path.join(case.directory, 'mask.txt'), nx, ny)
    for i, name in enumerate(field_names):
      case.fields[name].values = numpy.multiply(case.fields[name].values,
                                                masks[i])

    series['cases'].append(case)
  series['first'] = get_observed_orders(series['cases'][:-1],
                                        field_names,
                                        series['cases'][0],
                                        save_name=None)
  series['last'] = get_observed_orders(series['cases'][1:],
                                       field_names,
                                       series['cases'][0],
                                       save_name=None)

print()
labels = ['20x20', '25x25', '30x30', 'uquad', 'vquad']
for label in labels:
  print('- ' + label)
  for name in field_names:
    print('\t- {}: {:.3f}\t{:.3f}'.format(name,
                                          data[label]['first'][name],
                                          data[label]['last'][name]))
