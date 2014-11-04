#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/naca00xx.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: generates the coordinates of a NACA 00XX airfoil


import argparse
import math

import numpy as np
from matplotlib import pyplot as plt


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Generates a NACA00xx airfoil '
						'with uniform segments given the specified arguments',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('--t', dest='t', type=int, default=12,
						help='maximum thickness of the NACA airfoil (percent)')
	parser.add_argument('--c', dest='c', type=float, default=1.0,
						help='chord-length of the symmetric NACA airfoil')
	parser.add_argument('--ds', dest='ds', type=float, default=0.004,
						help='segment-length')
	parser.add_argument('--show', dest='show', action='store_true',
						help='plots the NACA00xx profile')
	return parser.parse_args()


def get_y_naca00xx(x, t, c):
	"""Returns the y-coordinate of the upper surface 
	of a symmetric NACA airfoil.
	
	Arguments
	---------
	x -- x-coordinate
	t -- maximum thickness as a fraction of the chord.
	c -- chord-length of the airfoil.
	"""
	t /= 100.
	return t/0.2*c*(  0.2969*np.sqrt(x/c)
					- 0.1260*(x/c)
					- 0.3516*(x/c)**2
					+ 0.2843*(x/c)**3
					- 0.1036*(x/c)**4 )


def plot_airfoil(x, y, t):
	"""Plots the profile in a figure.
	
	Arguments
	---------
	x, y -- coordinates of the profile.
	t -- maximum thickness of the arfoil as a percent of the chord-length.
	"""
	plt.figure()
	plt.grid(True)
	plt.xlabel(r'$x$', fontsize=18)
	plt.ylabel(r'$y$', fontsize=18)
	plt.plot(x, y, 
			 color='b', ls='-', lw=2, marker='o', markersize=6)
	plt.plot(x[::-1], -y[::-1], 
			 color='b', ls='-', lw=2, marker='o', markersize=6)
	plt.axis('equal')
	plt.title('NACA 00%02d' % t)
	plt.show()


def main():
	"""Generates the NACA 00xx coordinates."""
	# parse the command-line
	args = read_inputs()
	
	# NACA parameters
	t = args.t
	c = args.c
	ds = args.ds

	# compute the coordinates of the airfoil
	factor = 2						# refinement factor
	n = factor*int(c/ds)			# number of segments
	x = np.linspace(0.0, c, n+1)	# x-coordinates
	y = get_y_naca00xx(x, t, c)		# y-coordinates

	# length of the upper semi-airfoil
	length = np.sum(np.sqrt((x[1:]-x[:-1])**2 + (y[1:]-y[:-1])**2))

	# print parameters
	print 'Type of airfoil: NACA 00%02d' % t
	print 'Chord-length: %g' % c
	print 'Segment-length: %g' % ds
	print 'Length of semi-airfoil: %g' % length

	# number of boundary segments
	n = int(length/ds)+1

	# x-coordinate of end-points
	x0, xn = 0.0, c

	width = xn - x0

	# compute the coordinates of upper side of arfoil
	x = np.linspace(x0, xn, n+1)
	y = get_y_naca00xx(x, t, c)

	previous_mean, mean = 0.0, c
	tol = 1.0E-06
	while math.fabs(previous_mean - mean)/mean > tol:
		previous_mean = mean

		# compute segment-lengths
		ds = np.sqrt( (x[1:]-x[:-1])**2 + (y[1:]-y[:-1])**2 )

		mean = ds.sum()/n

		dx = mean/ds*(x[1:]-x[:-1])
		total = dx.sum()

		omega = width/total
		
		# update x-coordinates
		for i in xrange(n-1):
			x[i+1] = x[i] + omega*dx[i]
		y = get_y_naca00xx(x, t, c)

	# write coordinates into a file
	file_name = 'naca00%02d_%g.bdy' % (t, args.ds)
	with open(file_name, 'w') as outfile:
		outfile.write('%d\n' % (2*n))
		np.savetxt(outfile, np.c_[x-0.5*c, y], 
				   fmt='%.6f', delimiter='\t')
		np.savetxt(outfile, np.c_[x[-2:0:-1]-0.5*c, -y[-2:0:-1]], 
				   fmt='%.6f', delimiter='\t')

	# plot airfoil
	if args.show:
		plot_airfoil(x, y, t)


if __name__ == '__main__':
	main()
