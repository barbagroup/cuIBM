#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/analyzeVorticity.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Plots calculated and theoretical circulations at saved time-steps.


import os
import argparse
import math

import numpy
from matplotlib import pyplot

from readData import read_simulation_parameters, read_grid


def read_inputs():
	"""Parse the command-line"""
	# create the parser
	parser = argparse.ArgumentParser()
	# fill the parser with arguments
	parser.add_argument('--folder', dest='folder_path', type=str, 
						default=os.getcwd(),
						help='directory of the simulation')
	parser.add_argument('--time-increment', '-dt', dest='dt', type=float,
						help='time-increment of the simulation')
	parser.add_argument('--time-steps', '-t', dest='time_steps', type=int,
						nargs='+', default=[None, None, None],
						help='time-steps to plot (min, max, interval)')
	parser.add_argument('--limits', dest='limits', type=float, nargs='+', 
						default=[0.0, 10.0, -1.5, 1.5],
						help='x- and y- limits of the plot')
	parser.add_argument('--no-save', dest='save', action='store_false',
						help='does not save the figure')
	parser.add_argument('--show', dest='show', action='store_true',
						help='displays the figure')
	parser.set_defaults(save=True)
	return parser.parse_args()


def main():
	"""Plots the calculated and theoretical circulations at saved time-steps."""
	# parse the command-line
	args = read_inputs()

	# get the time-steps
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

	# read and generate the computational mesh
	nx, ny, dx, dy, _, yu, xv, _ = read_grid(args.folder_path)

	# grid-spacing
	h = dx[0]	# grid-spacing (grid should be uniformly spaced)
	# total number of iterations
	nt = time_steps[-1] - time_steps[0]

	# initialization
	times = args.dt * numpy.array(time_steps)
	total_vorticity = numpy.empty_like(times, dtype=float)
	circulation = numpy.empty_like(times, dtype=float)
	
	# loop over the saved time-steps
	for index, time_step in enumerate(time_steps):
		# read the vorticity file
		vorticity_file = '%s/%07d/vorticity' % (args.folder_path, time_step)
		with open(vorticity_file, 'r') as infile:
			vorticity = numpy.loadtxt(infile, dtype=float, usecols=(2,))
		# calculate the total vorticity
		total_vorticity[index] = h**2 * vorticity.sum()
		# calculate the theoretical circulation
		circulation[index] = -math.sin(math.pi*time_step/nt)

	# plot the total vorticity and theoretical circulation
	pyplot.figure()
	pyplot.grid(True)
	pyplot.xlabel('Time', fontsize=16)
	pyplot.ylabel('Total circulation', fontsize=16)
	pyplot.plot(times, total_vorticity, label='Total vorticity',
			 	color='r', ls='.', marker='o', markersize=10)
	pyplot.plot(times, circulation, label='Theoretical circulation',
			 	color='b', ls='-', lw=2)
	pyplot.xlim(args.limits[0], args.limits[1])
	pyplot.ylim(args.limits[2], args.limits[3])
	pyplot.title('Comparison of Total Vorticity and Theoretical Circulation '
			  	 'at Re=10,000')
	pyplot.legend(loc='best', prop={'size': 16})
	# save the figure
	if args.save:
		images_path = '%s/images' % args.folder_path
		if not os.path.isdir(images_path):
			os.makedirs(images_path)
		pyplot.savefig('%s/total_vorticity.png' % images_path)
	# display the figure
	if args.show:
		pyplot.show()


if __name__ == '__main__':
	main()
