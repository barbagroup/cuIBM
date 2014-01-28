#!/usr/bin/env python
import argparse
import numpy as np

# Generate optionsfile as per command line options
parser = argparse.ArgumentParser()

# Command line options
parser.add_argument("--filename", dest="filename", help="name of file containing body data", default="naca0012.bdy")

args = parser.parse_args()

filename = args.filename
x=[]
y=[]
nb = 0
f = open(filename, 'r')
a = f.readline().strip().split()
a = f.readline().strip().split()
x0 = float(a[0])
y0 = float(a[1])
xmin = x0
xmax = x0
a = f.readline().strip().split()
x1 = float(a[0])
y1 = float(a[1])
if x1 < xmin:
	xmin = x1
if x1 > xmax:
	xmax = x1
dsmin = np.sqrt( (x0-x1)**2 + (y0-y1)**2 )
dsmax = dsmin
while True:
	x0 = x1
	y0 = y1
	a = f.readline().strip().split()
	if a==[]:
		break
	x1 = float(a[0])
	y1 = float(a[1])
	if x1 < xmin:
		xmin = x1
	if x1 > xmax:
		xmax = x1
	temp = np.sqrt( (x0-x1)**2 + (y0-y1)**2 )
	if temp < dsmin:
			dsmin = temp;
	if temp > dsmax:
			dsmax = temp;
print dsmin, dsmax
print xmax-xmin
f.close()
	
