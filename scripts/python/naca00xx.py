#!/usr/bin/env python
import argparse
from numpy import *

def yNACA00xx(x, t=0.12, c=1.):
	return t*c/0.2*( 0.2969*sqrt(x/c) - 0.1260*(x/c) - 0.3516*(x/c)**2 + 0.2843*(x/c)**3 - 0.1036*(x/c)**4 )

# Read command line options
parser = argparse.ArgumentParser(description="Generates a NACA 00xx airfoil with uniform segments using the specified options.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--t", type=int, dest="t", help="thickness of symmetric NACA airfoil", default=8)
parser.add_argument("--c", type=float, dest="c", help="length of the chord", default=1.)
parser.add_argument("--ds", dest="ds", help="segment length", default='0.004')
args = parser.parse_args()

t = args.t/100.
c = args.c

# calculate the semi-circumference of the airfoil
nCirc = int(2*c/float(args.ds))+1
h = c/nCirc
C = 0.
x1 = 0.
y1 = -yNACA00xx(x1, t, c)
for i in range(1, nCirc+1):
	x2 = i*h
	y2 = -yNACA00xx(x2, t, c)
	dC = sqrt((x2-x1)**2 + (y2-y1)**2)
	C = C + dC
	x1 = x2
	y1 = y2

print "Type of airfoil: NACA 00%02d" % args.t
print "Chord Length   :", args.c
print "circumference  :", C
print "Segment Length :", args.ds

# number of boundary segments
n = int(C/float(args.ds))+1

# end points
x0 = 0
xN = args.c

width = xN-x0
temp = 0.
mean = args.c

x = linspace(x0, xN, n+1)
dx = zeros(n)
ds = zeros(n)

while fabs(temp-mean)/mean > 1e-6:
	temp = mean
	mean = 0.
	x1 = x[0]
	y1 = -yNACA00xx(x1, t, c)
	for i in range(n):
		x2 = x[i+1]
		y2 = -yNACA00xx(x2, t, c)
		ds[i] = sqrt((x2-x1)**2 + (y2-y1)**2)
		mean += ds[i]
		x1 = x2
		y1 = y2
	mean = mean/n
		
	total = 0.
	for i in range(n):
		dx[i] = mean/ds[i]*(x[i+1]-x[i])
		total += dx[i]
		
	omega = width/total
	
	for i in range(n-1):
		x[i+1] = x[i] + omega*dx[i]


filename = "naca%04d_%s.bdy" % (args.t, args.ds)
f = open(filename, 'w')
f.write("%d\n" % (2*n))
for i in range(n):
	f.write("%f\t%f\n" % ( x[i]-0.5*c, -yNACA00xx(x[i], t, c) ) )
for i in range(n):
	f.write("%f\t%f\n" % ( x[n-i]-0.5*c, yNACA00xx(x[n-i], t, c) ) )
f.close()
