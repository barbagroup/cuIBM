#!/usr/bin/env python

# file: $CUIBM_DIR/script/python/elements.py
# author: Anush Krishnan (anush@bu.edu), Oliver Mesnard (mesnardo@gwu.edu)
# description: display some properties of a given body


import os
import argparse
import math

import numpy as np
from matplotlib import pyplot as plt


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Display the properties '
						'of elements of a specified body file',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('--path', dest='path', type=str,
						help='path of the body file')
	parser.add_argument('--show', dest='show', action='store_true',
						help='plots the body')
	return parser.parse_args()


def plot_body(x, y, file_name):
	"""Plots the body.
	
	Arguments
	---------
	x, y -- coordinates of the body.
	file_name -- name of the body file.
	"""
	plt.figure()
	plt.grid(True)
	plt.xlabel(r'$x$', fontsize=18)
	plt.ylabel(r'$y$', fontsize=18)
	plt.plot(x, y, color='b', ls='-', lw=2, marker='o', markersize=4)
	plt.axis('equal')
	plt.title(file_name, fontsize=16)
	plt.show()


def main():
	"""Displays the properties of elements of a body."""
	# parse the command-line
	args = read_inputs()

	# read the body file
	with open(args.path, 'r') as infile:
		n = int(infile.readline())
		x, y = np.loadtxt(infile, dtype=float, delimiter='\t', unpack=True)

	# compute segment-lengths
	ds = np.append( np.sqrt((x[1:]-x[:-1])**2 + (y[1:]-y[:-1])**2),
					math.sqrt((x[-1]-x[0])**2 + (y[-1]-y[0])**2) )
	
	file_name = os.path.basename(args.path)	# basename of the path of the file

	# print properties of the body
	print '\tFile name: %s' % file_name
	print '\tFile path: %s' % args.path
	print '\tNumber of elements: %d' % n
	print '\tSize of largest element: %g' % ds.max()
	print '\tSize of smallest element: %g' % ds.min()
	print '\tx-width of the body: %g' % (x.max() - x.min())

	# plot the body
	if args.show:
		plot_body(x, y, file_name)


if __name__ == '__main__':
	main()
