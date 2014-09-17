#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/ellipse.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: script to generate an ellipse


import argparse
import math

import numpy
from matplotlib import pyplot


def read_inputs():
	"""Creates the parser ans parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Generates an elliptical body '
									 'with uniform segment legnths using the '
									 'the specified options',
									 formatter_class= \
									 argparse.ArgumentDefaultsHelpFormatter)
	# add arguments to the parser
	parser.add_argument('--a', dest='a', type=float, default=0.5,
						help='x-direction semi-axis length')
	parser.add_argument('--b', dest='b', type=float, default=0.06,
						help='y-direction semis-axis length')
	parser.add_argument('--x0', dest='x0', type=float, default=0.0,
						help='x-coordinate of the center of the ellipse')
	parser.add_argument('--y0', dest='y0', type=float, default=0.0,
						help='y-coordinate of the center of the ellipse')
	parser.add_argument('--ds', dest='ds', type=float, default=0.025,
						help='segment length')
	parser.add_argument('--theta0', dest='theta0', type=float, default=0.0,
						help='starting parameter')
	parser.add_argument('--thetaN', dest='thetaN', type=float, 
						default= 2.*math.pi,
						help='ending parameter')
	parser.add_argument('-f', '--filename', dest='filename', type=str,
						default='ell_0.025.bdy',
						help='name of the file generated')
	parser.add_argument('--show', dest='show', action='store_true',
						help='plots the ellipse using Matplotlib')
	return parser.parse_args()

def plot_ellipse(x, y):
	"""Plots the ellipse.
	
	Arguments
	---------
	x, y -- coordinates of the ellipse.
	"""
	pyplot.figure()
	pyplot.grid(True)
	pyplot.xlabel(r'$x$', fontsize=18)
	pyplot.ylabel(r'$y$', fontsize=18)
	pyplot.plot(x, y, color='b', ls='-', lw=2.0, marker='o', markersize=6)
	pyplot.xlim(x.min()-1, x.max()+1)
	pyplot.ylim(y.min()-1, y.max()+1)
	pyplot.show()

def main():
	"""Generates an ellipse."""
	# parse the command-line
	args = read_inputs()

	# parameters of the ellipse
	a, b = args.a, args.b
	x0, y0 = args.x0, args.y0

	# perimeter of the ellipse (Ramanujan approximation)
	h = ((a-b)/(a+b))**2
	p = math.pi*(a+b)*(1. + 3.*h/(10.+math.sqrt(4.-3.*h)))

	# number of boundary points
	n = int(math.ceil(p/args.ds))

	# end parameter points
	theta0, thetaN = args.theta0, args.thetaN

	angle = thetaN - theta0

	theta = numpy.linspace(theta0, thetaN, n+1)

	previous_mean, mean = 0.0, 1.0
	tol = 1.0E-06
	while math.fabs(previous_mean - mean)/mean > tol:
		previous_mean = mean

		ds = numpy.sqrt((a*(numpy.cos(theta[1:])-numpy.cos(theta[:-1])))**2
					   +(b*(numpy.sin(theta[1:])-numpy.sin(theta[:-1])))**2)
		mean = ds.sum()/n

		dtheta = mean/ds[:]*(theta[1:]-theta[:-1])
		total = dtheta.sum()

		omega = angle/total

		for i in xrange(n-1):
			theta[i+1] = theta[i] + omega*dtheta[i]

	# write the coordinates into a file
	with open(args.filename, 'w') as outfile:
		outfile.write('%d\n' % n)
		numpy.savetxt(outfile, numpy.c_[x0+a*numpy.cos(theta)[:-1], 
										y0+b*numpy.sin(theta)[:-1]],
					  fmt='%.6f', delimiter='\t')

	# plot the ellipse
	if args.show:
		plot_ellipse(x0+a*numpy.cos(theta), y0+b*numpy.sin(theta))

if __name__ == '__main__':
	main()
