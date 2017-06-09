"""
Computes the observed order of convergence.
"""

import math
import argparse
import numpy


parser = argparse.ArgumentParser()
parser.add_argument('--filepath',
                    dest='filepath',
                    type=str,
                    required=True,
                    help='file containing relative errors')
args = parser.parse_args()

with open(args.filepath, 'r') as infile:
  data = numpy.loadtxt(infile, unpack=True)

ratio = data[0][1] / data[0][0]
for d in data[2:]:
  for i in range(data[0].size - 2):
    p = math.log((d[i + 1] - d[i]) / (d[i + 2] - d[i + 1])) / math.log(ratio)
    print(p)
