#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/rotateBody.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: script to rotate a body


import argparse
import math

import numpy as np
from matplotlib import pyplot as plt


def read_inputs():
	# create the parser
	parser = argparse.ArgumentParser(description='Rotates a given body '
						'about a given center and generates a new body file',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('--alpha', dest='alpha', type=float, default=-4.0,
						help='angle of rotation')
	parser.add_argument('--xc', dest='xc', type=float, default=0.0,
						help='x-center of the center of rotation')
	parser.add_argument('--yc', dest='yc', type=float, default=0.0,
						help='y-center of the center of rotation')
	parser.add_argument('--infile', dest='infile', type=str, 
						default='naca0012_0.004.bdy',
						help='input file name')
	parser.add_argument('--show', dest='show', action='store_true',
						help='plots the geometry in both configurations')
	
	return parser.parse_args()


def plot_geometry(x, y, x_new, y_new):
	"""Plots the geometry in both the initial and rotated configuration.
	
	Arguments
	---------
	x, y -- coordinates of the initial geometry.
	x_new, y_new -- coordinates of the rotated geometry.
	"""
	plt.figure()
	plt.grid(True)
	plt.xlabel(r'$x$', fontsize=18)
	plt.ylabel(r'$y$', fontsize=18)
	plt.plot(x, y, label='initial', 
			 color='b', ls='-', lw=2, marker='o', markersize=6)
	plt.plot(x_new, y_new, label='rotated', 
			 color='r', ls='-', lw=2, marker='o', markersize=6)
	plt.legend()
	plt.show()


def main():
	"""Rotates a given geometry."""
	# parse the command-line
	args = read_inputs()

	# read the coordinates from a file
	with open(args.infile, 'r') as infile:
		n = int(infile.readline())
		x, y = np.loadtxt(infile, dtype=float, delimiter='\t', unpack=True)

	# parameters of the rotation
	x_c, y_c = args.xc, args.yc
	alpha = args.alpha*math.pi/180.

	# rotation
	x_new = x_c + (x-x_c)*math.cos(alpha) - (y-y_c)*math.sin(alpha)
	y_new = y_c + (x-x_c)*math.sin(alpha) + (y-y_c)*math.cos(alpha)

	# write the coordinates into a file
	file_name = '%s_%g.bdy' % (args.infile.split('.bdy')[0], args.alpha)
	with open(file_name, 'w') as outfile:
		outfile.write('%d\n' % n)
		np.savetxt(outfile, np.c_[x_new, y_new], fmt='%.6f', delimiter='\t')

	# plot the geometry in both configuration
	if args.show:
		plot_geometry(x, y, x_new, y_new)
		

if __name__ == '__main__':
	main()
