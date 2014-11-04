#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/plotPressure.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: plot the contours of pressure at every saved timestep

import os
import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from readData import readSimulationParameters, readGridData, readPressureData

def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='plots the pressure field '
						'at every save point for a given simulation',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('-folder', dest='folder', type=str, default='.',
						help='folder of the simulation')
	parser.add_argument("-xmin", type=float, dest="xmin", help="lower x-limit of the plotting region", default=-2.0)
	parser.add_argument("-xmax", type=float, dest="xmax", help="upper x-limit of the plotting region", default=4.0)
	parser.add_argument("-ymin", type=float, dest="ymin", help="lower y-limit of the plotting region", default=-3.0)
	parser.add_argument("-ymax", type=float, dest="ymax", help="upper y-limit of the plotting region", default=3.0)
	parser.add_argument('-plim', dest='plim', type=float, default=1.0,
						help='pressure cutoff on the plot')
	parser.add_argument('-numlevels', dest='numlevels', type=int, default=21,
						help='number of pressure contour line levels (choose an odd number)')
	return parser.parse_args()


def main():
	"""Plots the contour of vorticity at every time saved."""
	# parse the command-line
	args = read_inputs()

	folder = args.folder	# name of the folder

	# read the parameters of the simulation
	nt, start_step, nsave, _ = readSimulationParameters(folder)

	# calculate the mesh characteristics
	nx, ny, dx, dy, _, yu, xv, _ = readGridData(folder)

	# calculate appropriate array boundaries
	i_start = np.where(xv >= args.xmin)[0][0]
	i_end = np.where(xv <= args.xmax)[0][-1]
	j_start = np.where(yu >= args.ymin)[0][0]
	j_end = np.where(yu <= args.ymax)[0][-1]

	x = np.zeros(i_end+1-i_start)
	y = np.zeros(j_end+1-j_start)
	x[:] = xv[i_start:i_end+1]
	y[:] = yu[j_start:j_end+1]
	X, Y = np.meshgrid(x, y)

	P = np.zeros((j_end+1-j_start, i_end+1-i_start))

	# time-loop
	for ite in xrange(start_step+nsave, nt+1, nsave):
		# read the velocity data at the given time-step
		p = readPressureData(folder, ite, nx, ny)
		p = p.reshape((ny, nx))

		# calculate and write the vorticity
		for j in xrange(j_start, j_end+1):
			P[j-j_start, :] = p[j, i_start:i_end+1]

		CS = plt.contour(X, Y, P, levels=np.linspace(-args.plim, args.plim, args.numlevels))
		plt.title("Pressure")
		plt.colorbar(CS)
		plt.axis([xv[i_start], xv[i_end], yu[j_start], yu[j_end]])
		plt.gca().set_aspect('equal', adjustable='box')
		plt.savefig('{}/p{:0>7}.png'.format(folder,ite))
		plt.clf()
		print "Saved figure {}/p{:0>7}.png".format(folder,ite)
	
	print 'DONE!'

if __name__ == '__main__':
	main()
