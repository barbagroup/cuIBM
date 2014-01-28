#!/usr/bin/env python
import argparse
from numpy import *

# Generate optionsfile as per command line options
parser = argparse.ArgumentParser()
parser.add_argument("--alpha", dest="alpha", help="angle of rotation", default="-4")
parser.add_argument("--xc", type=float, dest="xc", help="x-coordinate of the center of rotation", default=0.0)
parser.add_argument("--yc", type=float, dest="yc", help="y-coordinate of the center of rotation", default=0.0)
parser.add_argument("--infile", dest="infile", help="name of input file", default="naca0008_0.008.bdy")
args = parser.parse_args()

alpha = float(args.alpha)*pi/180.
xC = args.xc
yC = args.yc
outFile = args.infile.split(".bdy")[0] + "_" + args.alpha + ".bdy"

f = open(args.infile, 'r')
g = open(outFile, 'w')

a = f.readline().strip().split()
nb = int(a[0])
g.write("%d\n" % nb)

xb = zeros(nb)
yb = zeros(nb)

for i in range(nb):
	a = f.readline().strip().split()
	x = float(a[0])
	y = float(a[1])
	xb[i] = xC + ( (x-xC)*cos(alpha) - (y-yC)*sin(alpha) )
	yb[i] = yC + ( (x-xC)*sin(alpha) + (y-yC)*cos(alpha) )
	g.write("%f\t%f\n" % (xb[i], yb[i]))

f.close()
g.close()
