#!/usr/bin/env python
import argparse
from numpy import *

# Generate optionsfile as per command line options
parser = argparse.ArgumentParser(description="Generates a NACA 00xx airfoil using the specified options.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Command line options
parser.add_argument("--t", type=int, dest="t", help="thickness of symmetric NACA airfoil", default=12)
parser.add_argument("--ds", type=float, dest="ds", help="x-distance between segments", default=0.004)

args = parser.parse_args()

# axes
t = args.t/100.
# number of boundary points
n = int(ceil(1/args.ds))
ds = 1./n

xx = zeros(2*n);
y = zeros(2*n);

for i in range(n):
	x = 1.-i*ds
	xx[i] = x
	y[i] = t/0.2*( 0.2969*sqrt(x) - 0.1260*x - 0.3516*x*x + 0.2843*x*x*x - 0.1015*x*x*x*x );

for i in range(n):
	x = i*ds
	xx[i+n] = x
	y[i+n] = -t/0.2*( 0.2969*sqrt(x) - 0.1260*x - 0.3516*x*x + 0.2843*x*x*x - 0.1015*x*x*x*x );

filename = "naca00%02d.bdy" % args.t
print filename
f = open(filename, 'w')
f.write("%d\n" % (2*n))
for i in range(2*n):
	f.write("%f\t%f\n" % ( xx[i], y[i] ) )
f.close()
