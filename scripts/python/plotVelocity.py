#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/plotVelocity.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: plot the contour of velocity


import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from readData import readSimulationParameters, readGridData, readVelocityData


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Plots the velocity field at '
						'all save points for a given simulation',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('--folder', dest='folder', type=str, default='.',
						help='folder of the simulation')
	parser.add_argument('--blx', dest='bl_x', type=float, default=-2.0,
						help='x-coordinate of the bottom-left corner'
							 'of the plot region')
	parser.add_argument('--bly', dest='bl_y', type=float, default=-3.0,
						help='y-coordinate of the bottom-left corner'
							 'of the plot region')
	parser.add_argument('--trx', dest='tr_x', type=float, default=4.0,
						help='x-coordinate of the top-right corner'
							 'of the plot region')
	parser.add_argument('--try', dest='tr_y', type=float, default=3.0,
						help='y-coordinate of the top-right corner'
							 'of the plot region')
	parser.add_argument('--ulim', dest='u_lim', type=float, default=1.0,
						help='x-velocity cutoff on the plot')
	parser.add_argument('--vlim', dest='v_lim', type=float, default=1.0,
						help='y-velocity cutoff on the plot')
	
	return parser.parse_args()


def main():
	"""Plots the contour of velocity (u and v) at every time saved."""
	# parse the command-line
	args = read_inputs()

	folder = args.folder	# name of the folder

	# read the parameters of the simulation
	nt, start_step, nsave, _ = readSimulationParameters(folder)

	# calculate the mesh characteristics
	nx, ny, dx, dy, xu, yu, xv, yv = readGridData(folder)

	# generate a mesh grid for u- and v- velocities
	Xu, Yu = np.meshgrid(xu, yu)
	Xv, Yv = np.meshgrid(xv, yv)

	for ite in xrange(start_step+nsave, nt+1, nsave):
		print 'iteration %d' % ite
		# read the velocity data at the given time-step
		u, v = readVelocityData(folder, ite, nx, ny, dx, dy)
		if u == None or v == None:
			break

		# plot u-velocity contour
		plt.contour(Xu, Yu, u.reshape((ny, nx-1)),
					levels=np.linspace(-args.u_lim, args.u_lim, 21))
		plt.axis('equal')
		plt.axis([args.bl_x, args.tr_x, args.bl_y, args.tr_y])
		plt.colorbar()
		plt.savefig('%s/u%07d.png' % (folder, ite))
		plt.clf()
	
		# plot v-velocity contour
		plt.contour(Xv, Yv, v.reshape((ny-1, nx)), 
					levels=np.linspace(-args.v_lim, args.v_lim, 21))
		plt.axis('equal')
		plt.axis([args.bl_x, args.tr_x, args.bl_y, args.tr_y])
		plt.colorbar()
		plt.savefig('%s/v%07d.png' % (folder, ite))
		plt.clf()

if __name__ == '__main__':
	main()
