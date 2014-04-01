#!/usr/bin/env python
import argparse
from numpy import *

# Generate optionsfile as per command line options
parser = argparse.ArgumentParser()

# Command line options
parser.add_argument("--a", type=float, dest="a", help="semi-axis along X", default=0.5)
parser.add_argument("--b", type=float, dest="b", help="semi-axis along Y", default=0.06)
parser.add_argument("--x0", type=float, dest="x0", help="X coordinate of ellipse center", default=0.)
parser.add_argument("--y0", type=float, dest="y0", help="Y coordinate of ellipse center", default=0.)
parser.add_argument("--ds", type=float, dest="ds", help="segment length", default=0.025)
parser.add_argument("--theta0", type=float, dest="t0", help="starting parameter", default=0.)
parser.add_argument("--thetaN", type=float, dest="tN", help="ending parameter", default=2*pi)
parser.add_argument("--filename", dest="filename", help="name of file generated", default="ell_0.025.bdy")

args = parser.parse_args()

# axes
a = args.a
b = args.b
x0 = args.x0
y0 = args.y0
# number of boundary points
C = pi*(a+b)*(1 + 3*((a-b)/(a+b))**2/(10 + sqrt(4 - 3*((a-b)/(a+b))**2))) # circumference
n = int(ceil(C/args.ds))
# end points
theta0 = args.t0
thetaN = args.tN

angle = thetaN - theta0
temp = 0.
mean = 1.

theta = linspace(theta0, thetaN, n+1)
dtheta = zeros(n)
ds = zeros(n)

while fabs(temp-mean)/mean > 1e-6:
	
	temp = mean
	mean = 0.
	for i in range(n):
		dx = a*( cos(theta[i+1]) - cos(theta[i]) )
		dy = b*( sin(theta[i+1]) - sin(theta[i]) )
		ds[i] = sqrt(dx*dx+dy*dy)
		mean += ds[i]
	mean = mean/n
		
	total = 0.
	for i in range(n):
		dtheta[i] = mean/ds[i]*(theta[i+1]-theta[i])
		total += dtheta[i]
		
	omega = angle/total
	
	for i in range(n-1):
		theta[i+1] = theta[i] + omega*dtheta[i]

f = open(args.filename, 'w')
f.write("%d\n" % n)
for i in range(n):
	f.write("%f\t%f\n" % ( x0+a*cos(theta[i]), y0+b*sin(theta[i]) ) )
f.close()