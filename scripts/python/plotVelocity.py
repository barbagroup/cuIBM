#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/plotVelocity.py
# author: Krishnan, A. (anush@bu.edu)
# description: plot the contour of velocity


import argparse

import numpy as np
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
	parser.add_argument('--vortlim', dest='vort_lim', type=float, default=15.0,
						help='vorticity cutoff on the plot')
	parser.add_argument('--plot-only', dest='plot_only', action='store_true',
						help='only generate the vorticity plots '
							 'and do not recalculate the data')
	
	return parser.parse_args()


def main():
	"""Plots the contour of velocity (u and v) at every time saved."""
	# parse the command-line
	args = read_inputs()

	folder = args.folder	# name of the folder

	# read the parameters of the simulation
	nt, start_step, nsave, _ = readSimulationParameters(folder)

	if not args.plot_only:
		# calculate the mesh characteristics
		nx, ny, dx, dy, xu, yu, xv, yv = readGridData(folder)

		# generate a mesh grid for u- and v- velocities
		Xu, Yu = np.meshgrid(xu, yu)
		Xv, Yv = np.meshgrid(xv, yv)

		# WARNING: useless in this code
		i_start = np.where(xv >= args.bl_x)[0][0]
		i_end = np.where(xv <= args.tr_x)[0][0]
		j_start = np.where(yu >= args.bl_y)[0][0]
		j_end = np.where(yu <= args.tr_y)[0][0]

		ite = start_step + nsave
		while ite < nt+1:
			print 'iteration %d' % ite
			# read the velocity data at the given time-step
			u, v = readVelocityData(folder, ite, nx, ny, dx, dy)
			if u == None or v == None:
				break

			# plot u-velocity contour
			plt.contour(Xu, Yu, u.reshape((ny, nx-1)),
						levels=np.arange(-1.0, 1.0, 0.1))
			plt.axis('equal')
			plt.axis([args.bl_x, args.tr_x, args.bl_y, args.tr_y])
			plt.colorbar()
			plt.savefig('%s/u%07d.png' % (folder, ite))
			plt.clf()
		
			# plot v-velocity contour
			plt.contour(Xv, Yv, v.reshape((ny-1, nx)), 
						levels=np.arange(-1.0, 1.0, 0.1))
			plt.axis('equal')
			plt.axis([args.bl_x, args.tr_x, args.bl_y, args.tr_y])
			plt.colorbar()
			plt.savefig('%s/v%07d.png' % (folder, ite))
			plt.clf()

			ite += nsave


if __name__ == '__main__':
	main()
