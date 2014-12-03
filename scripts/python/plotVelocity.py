#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/plotVelocity.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: plot the contours of velocity


import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.font_manager as font_manager

from readData import readSimulationParameters, readGridData, readVelocityData


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Plots the velocity field at '
						'all save points for a given simulation',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('-folder', dest='folder', type=str, default='.',
						help='folder of the simulation')
	parser.add_argument("-xmin", type=float, dest="xmin", help="lower x-limit of the plotting region", default=-2.0)
	parser.add_argument("-xmax", type=float, dest="xmax", help="upper x-limit of the plotting region", default=4.0)
	parser.add_argument("-ymin", type=float, dest="ymin", help="lower y-limit of the plotting region", default=-3.0)
	parser.add_argument("-ymax", type=float, dest="ymax", help="upper y-limit of the plotting region", default=3.0)
	parser.add_argument('-ulim', dest='u_lim', type=float, default=1.0,
						help='x-velocity cutoff on the plot')
	parser.add_argument('-vlim', dest='v_lim', type=float, default=1.0,
						help='y-velocity cutoff on the plot')
	
	return parser.parse_args()


def main():

	#view available fonts
	#for font in font_manager.findSystemFonts():
	#	print font

	try:
		fontpath = '/usr/share/fonts/truetype/msttcorefonts/Georgia.ttf'

		prop = font_manager.FontProperties(fname=fontpath)
		matplotlib.rcParams['font.family'] = prop.get_name()
	except:
		pass

	"""Plots the contour of velocity (u and v) at every time saved."""
	# parse the command-line
	args = read_inputs()

	folder = args.folder	# name of the folder

	# read the parameters of the simulation
	nt, start_step, nsave, _ = readSimulationParameters(folder)

	# calculate the mesh characteristics
	nx, ny, dx, dy, xu, yu, xv, yv = readGridData(folder)

	# calculate appropriate array boundaries
	i_start = np.where(xu >= args.xmin)[0][0]
	i_end = np.where(xu <= args.xmax)[0][-1]
	j_start = np.where(yu >= args.ymin)[0][0]
	j_end = np.where(yu <= args.ymax)[0][-1]

	x_start = xu[i_start] - dx[i_start]
	x_end   = xu[i_end] + dx[i_end+1]
	y_start = yu[j_start] - dy[j_start]/2.
	y_end   = yu[j_end] + dy[j_end]/2.

	# generate a mesh grid for u- and v- velocities
	Xu, Yu = np.meshgrid(xu, yu)
	Xv, Yv = np.meshgrid(xv, yv)

	for ite in xrange(start_step+nsave, nt+1, nsave):
		print 'iteration %d' % ite
		# read the velocity data at the given time-step
		u, v = readVelocityData(folder, ite, nx, ny, dx, dy)
		if u is None or v is None:
			break

		# plot u-velocity contour
		fig = plt.figure()
		ax = fig.add_subplot(111)
		CS = ax.contour(Xu, Yu, u.reshape((ny, nx-1)),
					levels=np.linspace(-args.u_lim, args.u_lim, 21))
		fig.gca().set_aspect('equal', adjustable='box')
		ax.set_xlim([x_start,x_end])
		ax.set_ylim([y_start,y_end])
		ax.tick_params(labelsize=16)
		cbar = fig.colorbar(CS)
		cbar.ax.tick_params(labelsize=16) 
		fig.savefig('%s/u%07d.pdf' % (folder, ite))
		fig.clf()
	
		# plot v-velocity contour
		fig = plt.figure()
		ax = fig.add_subplot(111)
		CS = ax.contour(Xv, Yv, v.reshape((ny-1, nx)), 
					levels=np.linspace(-args.v_lim, args.v_lim, 21))
		fig.gca().set_aspect('equal', adjustable='box')
		ax.set_xlim([x_start,x_end])
		ax.set_ylim([y_start,y_end])
		ax.tick_params(labelsize=16)
		cbar = fig.colorbar(CS)
		cbar.ax.tick_params(labelsize=16) 
		fig.savefig('%s/v%07d.pdf' % (folder, ite))
		fig.clf()

if __name__ == '__main__':
	main()
