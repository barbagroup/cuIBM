#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/plotVorticity.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Plots the contour of vorticity at saved time-steps.


import os
import argparse

import numpy
from matplotlib import pyplot

from readData import read_grid, read_velocity


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
	parser.add_argument('--time-steps', '-t', dest='time_steps', type=int, 
						nargs='+', default=[None, None, None],
						help='time-steps to plot (min, max, interval)')
	parser.add_argument('--bottom-left', '-bl', dest='bottom_left', type=float,
						nargs='+', default=[-2.0, -3.0],
						help='bottom-left coordinates of the rectangular view')
	parser.add_argument('--top-right', '-tr', dest='top_right', type=float,
						nargs='+', default=[4.0, 3.0],
						help='top-right coordinates of the rectangular view')
	parser.add_argument('--vorticity-limit', '-vl', dest='vorticity_limit', 
						type=float, default=3.0,
						help='upper limit of zero-symmetric vorticity range')
	parser.add_argument('--levels', '-l', dest='levels', type=int, default=16,
						help='number of contour line levels')
	parser.add_argument('--no-save', dest='save', action='store_false',
						help='saves the vorticity data in the time-step folder')
	parser.add_argument('--gnuplot', dest='using_gnuplot', action='store_true',
						help='plots contours using Gnuplot instead of Pyplot')
	parser.set_defaults(save=True)
	return parser.parse_args()


def main():
	"""Plots the contour of vorticity at saved time-steps."""
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

	# create images folder if does not exist to store PNG files
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
	x = 0.5*(xv[i_start:i_end] + xv[i_start+1:i_end+1])
	y = 0.5*(yu[j_start:j_end] + yu[j_start+1:j_end+1])
	X, Y = numpy.meshgrid(x, y)

	# initialize vorticity
	vorticity = numpy.empty((y.size, x.size))

	for time_step in time_steps:
		# read the velocity data at the given time-step
		u, v = read_velocity(args.folder_path, time_step, nx, ny, dx, dy)

		# calculate the vorticity
		for j in xrange(j_start, j_end):
			Dy = 0.5 * (dy[j] + dy[j+1])
			for i in xrange(i_start, i_end):
				Dx = 0.5 * (dx[i] + dx[i+1])
				vorticity[j-j_start, i-i_start] = (v[j*nx+i+1]-v[j*nx+i])/Dx \
										- (u[(j+1)*(nx-1)+i]-u[j*(nx-1)+i])/Dy
		
		if args.save:
			# save the vorticity data in the time-step folder
			print 'Saving vorticity data at time-step %d...' % time_step
			vorticity_file = '{}/{:0>7}/vorticity'.format(args.folder_path, 
			                                              time_step)
			with open(vorticity_file, 'w') as outfile:
				for j in xrange(y.size):
					for i in xrange(x.size):
						outfile.write('%f\t%f\t%f\n' 
									  % (x[i], y[j], vorticity[j, i]))
					outfile.write('\n')
		if not args.using_gnuplot:
			# plot the contour of vorticity using pyplot
			print ('Generating PNG file with Pyplot at time-step %d...' 
				   % time_step)
			pyplot.figure()
			pyplot.xlabel(r'$x$', fontsize=18)
			pyplot.ylabel(r'$y$', fontsize=18)
			pyplot.xlim(args.bottom_left[0], args.top_right[0])
			pyplot.ylim(args.bottom_left[1], args.top_right[1])
			pyplot.axis('equal')
			levels = numpy.linspace(-args.vorticity_limit, args.vorticity_limit,
									args.levels)
			cont = pyplot.contour(X, Y, vorticity, levels)
			cont_bar = pyplot.colorbar(cont)
			cont_bar.set_label('vorticity')
			pyplot.savefig('{}/o{:0>7}.png'.format(images_path, time_step))
			pyplot.clf()
			pyplot.close()

	if args.using_gnuplot:
		# create the gnuplot script file
		print 'Creating Gnuplot vorticity script...'
		gnuplot_file = '{}/vorticity_script.plt'.format(args.folder_path)
		with open(gnuplot_file, 'w') as outfile:
			outfile.write('reset;\n')
			outfile.write('set terminal pngcairo enhanced '
						  'font "Times, 15" size 900,600;\n')
			for time_step in time_steps:
				outfile.write('\nset output "%s/o%07d.png"\n' % (images_path,
																	time_step))
				outfile.write('set multiplot;\n')
				outfile.write('reset;\n')
				outfile.write('set view map; set size ratio -1; unset key;\n')
				outfile.write('set pm3d map;\n')
				outfile.write('set palette defined '
							  '(-2 "dark-blue", -1 "light-blue", 0 "white", '
							  '1 "light-red", 2 "dark-red");\n')
				outfile.write('set cbrange [%f:%f];\n' 
							  % (-args.vorticity_limit, args.vorticity_limit))
				outfile.write('splot [%f:%f] [%f:%f] "%s/%07d/vorticity";\n'
							  % (args.bottom_left[0], args.top_right[0],
							     args.bottom_left[1], args.top_right[1],
								 args.folder_path, time_step))
				outfile.write('unset multiplot;\n')
		print 'Generating PNG files with Gnuplot...'
		os.system('gnuplot %s' % gnuplot_file)

	print '\nVorticity contours: DONE!\n'


if __name__ == '__main__':
	main()
