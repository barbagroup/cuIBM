#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/plotPressure.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Plots the contour of pressure at saved time-steps


import os
import argparse

import numpy
from matplotlib import pyplot

from readData import read_grid, read_pressure


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Plots the contour of '
						'pressure at saved time-steps',
						formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# fill the parser with arguments
	parser.add_argument('--folder', dest='folder_path', type=str, 
						default=os.getcwd(),
						help='directory of the simulation')
	parser.add_argument('--time-steps', '-t', dest='time_steps', type=int,
						nargs='+', default=[None, None, None],
						help='time-steps to plot (min, max, interval)')
	parser.add_argument('--bottom-left', '-bl', dest='bottom_left', type=float,
						nargs='+', default=[-2.0, -3.0],
						help='bottom-left coordinates of the rectangular view')
	parser.add_argument('--top-right', '-tr', dest='top_right', type=float,
						nargs='+', default=[4.0, 3.0],
						help='top-right coordinates of the rectangular view')
	parser.add_argument('--pressure-limits', '-pl', dest='pressure_limits',
						type=float, nargs='+', default=[-1.0, 1.0],
						help='limits of the pressure range to plot')
	parser.add_argument('--levels', '-l', dest='levels', type=int, default=21,
						help='number of contour line levels')
	return parser.parse_args()


def main():
	"""Plots the contour of vorticity at every time saved."""
	# parse the command-line
	args = read_inputs()

	# get the time-steps to plot
	if any(args.time_steps):
		# if provided by command-line arguments
		time_steps = range(args.time_steps[0],
						   args.time_steps[1]+1,
						   args.time_steps[2])
	else:
		# if not, list solution folders
		time_steps = sorted(int(folder) for folder
										in os.listdir(args.folder_path)
										if folder[0]=='0')
	
	# create images folder if does not exit to store PNG files
	images_path = '%s/images' % args.folder_path
	if not os.path.isdir(images_path):
		os.makedirs(images_path)

	# calculate the mesh characteristics
	nx, ny, dx, dy, _, yu, xv, _ = read_grid(args.folder_path)

	# calculate appropriate array boundaries
	i_start = numpy.where(xv >= args.bottom_left[0])[0][0]
	i_end = numpy.where(xv <= args.top_right[0])[0][-1]
	j_start = numpy.where(yu >= args.bottom_left[1])[0][0]
	j_end = numpy.where(yu <= args.top_right[1])[0][-1]

	# create a mesh-grid
	x = xv[i_start:i_end]
	y = yu[j_start:j_end]
	X, Y = numpy.meshgrid(x, y)

	for time_step in time_steps:
		# read the velocity data at the given time-step
		pressure = read_pressure(args.folder_path, time_step, nx, ny)
		pressure = pressure.reshape((ny, nx))[j_start:j_end, i_start:i_end]

		# plot the contour of pressure
		print 'Generating PNG file at time-step %d...' % time_step
		pyplot.figure()
		pyplot.xlabel(r'$x$', fontsize=18)
		pyplot.ylabel(r'$y$', fontsize=18)
		pyplot.xlim(args.bottom_left[0], args.top_right[0])
		pyplot.ylim(args.bottom_left[1], args.top_right[1])
		pyplot.axis('equal')
		levels = numpy.linspace(args.pressure_limits[0], 
								args.pressure_limits[1],
								args.levels)
		cont = pyplot.contour(X, Y, pressure, levels)
		cont_bar = pyplot.colorbar(cont)
		cont_bar.set_label('pressure')
		pyplot.savefig('{}/p{:0>7}.png'.format(images_path, time_step))
		pyplot.clf()
		pyplot.close()

	print '\nPressure contours: DONE!\n'


if __name__ == '__main__':
	main()
