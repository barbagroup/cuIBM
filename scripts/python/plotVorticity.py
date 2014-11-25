#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/plotVorticity.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: plot the contours of vorticity at every saved timestep

import os
import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from readData import readSimulationParameters, readGridData, readVelocityData

def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='plots the vorticity field '
						'at every save point for a given simulation',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('-folder', dest='folder', type=str, default='.',
						help='folder of the simulation')
	parser.add_argument("-xmin", type=float, dest="xmin", help="lower x-limit of the plotting region", default=-2.0)
	parser.add_argument("-xmax", type=float, dest="xmax", help="upper x-limit of the plotting region", default=4.0)
	parser.add_argument("-ymin", type=float, dest="ymin", help="lower y-limit of the plotting region", default=-3.0)
	parser.add_argument("-ymax", type=float, dest="ymax", help="upper y-limit of the plotting region", default=3.0)
	parser.add_argument('-vortlim', dest='vortlim', type=float, default=3.0,
						help='vorticity cutoff on the plot')
	parser.add_argument('-numlevels', dest='numlevels', type=int, default=16,
						help='number of vortex contour line levels (choose an even number)')
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

	y = np.zeros(j_end-j_start)
	x = np.zeros(i_end-i_start)
	x[:] = 0.5*(xv[i_start:i_end] + xv[i_start+1:i_end+1])
	y[:] = 0.5*(yu[j_start:j_end] + yu[j_start+1:j_end+1])
	X, Y = np.meshgrid(x, y)

	Omg = np.zeros((j_end-j_start, i_end-i_start))

	# time-loop
	for ite in xrange(start_step+nsave, nt+1, nsave):
		# read the velocity data at the given time-step
		u, v = readVelocityData(folder, ite, nx, ny, dx, dy)
		if u is None or v is None:
			break

		# calculate the vorticity
		for j in xrange(j_start, j_end):
			Dy = 0.5 * (dy[j] + dy[j+1])
			for i in xrange(i_start, i_end):
				Dx = 0.5 * (dx[i] + dx[i+1])
				Omg[j-j_start, i-i_start] = (v[j*nx+i+1] - v[j*nx+i]) / Dx \
					  - (u[(j+1)*(nx-1)+i] - u[j*(nx-1)+i]) / Dy

		CS = plt.contour(X, Y, Omg, levels=np.linspace(-args.vortlim, args.vortlim, args.numlevels))
		plt.title("Vorticity")
		plt.colorbar(CS)
		plt.axis([xv[i_start], xv[i_end], yu[j_start], yu[j_end]])
		plt.gca().set_aspect('equal', adjustable='box')
		plt.savefig('{}/o{:0>7}.png'.format(folder,ite))
		plt.clf()
		print "Saved figure {}/O{:0>7}.png".format(folder,ite)
	
	print 'DONE!'

if __name__ == '__main__':
	main()
