#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/plotVorticity.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Plots the contour of vorticity at saved time-steps.


import os
import argparse

import numpy
from matplotlib import pyplot

from readData import readSimulationParameters, readGridData, readVelocityData


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Plots the contour of '
						'vorticity at saved time-steps',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('--folder', dest='folder_path', type=str, 
						default=os.getcwd(),
						help='directory of the simulation')
	parser.add_argument('--iterations', '-i', dest='iterations', type=float, 
						nargs='+', default=[None, None, None],
						help='time-steps to plot (min, max, interval)')
	parser.add_argument('--bottom-left', '-bl', dest='bottom_left', type=float,
						default=[-2.0, -3.0],
						help='bottom-left coordinates of the rectangular view')
	parser.add_argument('--top-right', '-tr', dest='top_right', type=float,
						default=[4.0, 3.0],
						help='top-right coordinates of the rectangular view')
	parser.add_argument('--vorticity-limit', '-vl', dest='vorticity_limit', 
						type=float, default=3.0,
						help='upper limit of zero-symmetric vorticity range')
	parser.add_argument('--levels', '-l', dest='levels', type=int, default=16,
						help='number of contour line levels '
							 '(choose an even number)')
	return parser.parse_args

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
	"""Plots the contour of vorticity at saved time-steps."""
	# parse the command-line
	args = read_inputs()

	folder_path = args.folder_path	# name of the folder

	# get the time-steps to plot
	if not any(args.time_steps):
		time_steps = [args.time_steps[0]:args.time_steps[1]:args.time_steps[2]]
	else:
		time_steps = sorted(int(folder) for folder in os.listdir(folder_path)
										if folder[0]=='0')

	# calculate the mesh characteristics
	nx, ny, dx, dy, _, yu, xv, _ = readGridData(args.folder_path)

	# calculate appropriate array boundaries
	i_start = numpy.where(xv >= args.bottom_left[0])[0][0]
	i_end = numpy.where(xv <= args.top_right[0])[0][-1]
	j_start = numpy.where(yu >= args.bottom_left[1])[0][0]
	j_end = numpy.where(yu <= args.top_right[1])[0][-1]

	# create a mesh-grid
	x = 0.5*(xv[i_start:i_end] + xv[i_start+1:i_end+1])
	y = 0.5*(yu[j_start:j_end] + yu[j_start+1:j_end+1])
	X, Y = numpy.meshgrid(x, y)

	vorticity = numpy.zeros((y.size, x.size))

	for time_step in time_steps:
		# read the velocity data at the given time-step
		u, v = readVelocityData(folder_path, time_step, nx, ny, dx, dy)

	Omg = np.zeros((j_end-j_start, i_end-i_start))

	# time-loop
	for ite in xrange(start_step+nsave, start_step+nt+1, nsave):
		# read the velocity data at the given time-step
		u, v = readVelocityData(folder, ite, nx, ny, dx, dy)
		if u == None or v == None:
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
