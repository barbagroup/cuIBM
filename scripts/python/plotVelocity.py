#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/plotVelocity.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Plots the contour of velocity at saved time-steps.


import os
import argparse

import numpy
from matplotlib import pyplot

from readData import read_grid, read_velocity


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='Plots the contour of '
						'velocity at saved time-steps',
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
	parser.add_argument('--velocity-limits', '-vl', dest='velocity_limits',
						type=float, nargs='+', default=[-1.0, 1.0, -1.0, 1.0],
						help='limits of the u-velocity and v-velocity ranges')
	parser.add_argument('--levels', '-l', dest='levels', type=int, nargs='+',
						default=[21, 21],
						help='number of contour line levels '
							 'for u- and v- velocities')
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

	# create image folder if does not exist to store PNG files
	images_path = '%s/images' % args.folder_path
	if not os.path.isdir(images_path):
		os.makedirs(images_path)

	# calculate the mesh characteristics
	nx, ny, dx, dy, xu, yu, xv, yv = read_grid(args.folder_path)

	# generate a mesh grid for u- and v- velocities
	Xu, Yu = numpy.meshgrid(xu, yu)
	Xv, Yv = numpy.meshgrid(xv, yv)

	for time_step in time_steps:
		# read the velocity data at the given time-step
		u, v = read_velocity(args.folder_path, time_step, nx, ny, dx, dy)

		print 'Generating PNG file at time-step %d...' % time_step
		# plot u-velocity contour
		print '\tu-velocity...'
		pyplot.figure()
		pyplot.xlabel(r'$x$', fontsize=18)
		pyplot.ylabel(r'$y$', fontsize=18)
		levels = numpy.linspace(args.velocity_limits[0],
								args.velocity_limits[1],
								args.levels[0])
		cont = pyplot.contour(Xu, Yu, u.reshape((ny, nx-1)), levels)
		cont_bar = pyplot.colorbar(cont)
		cont_bar.set_label('u-velocity')
		pyplot.xlim(args.bottom_left[0], args.top_right[0])
		pyplot.ylim(args.bottom_left[1], args.top_right[1])
		pyplot.savefig('{}/u{:0>7}.png'.format(images_path, time_step))
		pyplot.clf()
		pyplot.close()
		# plot v-velocity contour
		print '\tv-velocity...'
		pyplot.figure()
		pyplot.xlabel(r'$x$', fontsize=18)
		pyplot.ylabel(r'$y$', fontsize=18)
		levels = numpy.linspace(args.velocity_limits[2],
								args.velocity_limits[3],
								args.levels[1])
		cont = pyplot.contour(Xv, Yv, v.reshape((ny-1, nx)), levels)
		cont_bar = pyplot.colorbar(cont)
		cont_bar.set_label('v-velocity')
		pyplot.xlim(args.bottom_left[0], args.top_right[0])
		pyplot.ylim(args.bottom_left[1], args.top_right[1])
		pyplot.savefig('{}/v{:0>7}.png'.format(images_path, time_step))
		pyplot.clf()
		pyplot.close()

	print '\nVelocity contours: DONE!\n'


if __name__ == '__main__':
	main()
