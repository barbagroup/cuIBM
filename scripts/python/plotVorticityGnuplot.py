#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/plotVorticity.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: plot the contour of vorticity at every time saved


import os
import argparse

import numpy as np

from readData import readSimulationParameters, readGridData, readVelocityData


def read_inputs():
	"""Parses the command-line."""
	# create the parser
	parser = argparse.ArgumentParser(description='plots the vorticity field '
						'at every save points for a given simulation',
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
	"""Plots the contour of vorticity at every time saved."""
	# parse the command-line
	args = read_inputs()

	folder = args.folder	# name of the folder

	# read the parameters of the simulation
	nt, start_step, nsave, _ = readSimulationParameters(folder)

	if not args.plot_only:
		# calculate the mesh characteristics
		nx, ny, dx, dy, _, yu, xv, _ = readGridData(folder)

		# calculate appropriate array boundaries
		i_start = np.where(xv >= args.bl_x)[0][0]
		i_end = np.where(xv <= args.tr_x)[0][-1]
		j_start = np.where(yu >= args.bl_y)[0][0]
		j_end = np.where(yu <= args.tr_y)[0][-1]

		# time-loop
		for ite in xrange(start_step+nsave, nt+1, nsave):
			# read the velocity data at the given time-step
			u, v = readVelocityData(folder, ite, nx, ny, dx, dy)
			if u is None or v is None:
				break

			# calculate and write the vorticity
			vort_file = '%s/%07d/vorticity' % (folder, ite)
			print vort_file
			with open(vort_file, 'w') as outfile:
				for j in xrange(j_start, j_end):
					y = 0.5 * (yu[j] + yu[j+1])
					Dy = 0.5 * (dy[j] + dy[j+1])
					for i in xrange(i_start, i_end):
						x = 0.5 * (xv[i] + xv[i+1])
						Dx = 0.5 * (dx[i] + dx[i+1])
						vort = (v[j*nx+i+1] - v[j*nx+i]) / Dx \
							 - (u[(j+1)*(nx-1)+i] - u[j*(nx-1)+i]) / Dy
						outfile.write('%f\t%f\t%f\n' % (x, y, vort))
					outfile.write('\n')

	# create the gnuplot file
	print 'Creating gnuplot file...'
	gnuplot_file = '%s/vorticity.plt' % folder
	with open(gnuplot_file, 'w') as outfile:
		outfile.write('reset;\n')
		outfile.write('set terminal pngcairo enhanced '
					  'font "Times, 15" size 900,600;\n')

		for ite in xrange(start_step+nsave, nt+1, nsave):
			# image file name
			outfile.write('\nset output "%s/vort%07d.png"\n' % (folder, ite))
			outfile.write('set multiplot;\n')
			# vorticity
			outfile.write('reset;\n')
			outfile.write('set view map; set size ratio -1; unset key;\n')
			outfile.write('set pm3d map;\n')
			outfile.write('set palette defined (-2 "dark-blue", '
											   '-1 "light-blue", '
											   '0 "white", '
											   '1 "light-red", '
											   '2 "dark-red");\n')
			outfile.write('set cbrange [%f:%f];\n' 
						  % (-args.vort_lim, args.vort_lim))
			outfile.write('splot [%f:%f] [%f:%f] "%s/%07d/vorticity";\n'
						  % (xv[i_start], xv[i_end], yu[j_start], yu[j_end], folder, ite))
			# bodies
			#outfile.write('reset;\n')
			#outfile.write('set view map; set size ratio -1; unset key;\n')
			#outfile.write('splot [%f:%f] [%f:%f] "%s/%07d/bodies" '
			#			  'u 1:2:(0) wl lw 2 lc "black";\n'
			#			  % (args.bl_x, args.tr_x, args.bl_y, args.tr_y, 
			#			  	 folder, ite))
			outfile.write('unset multiplot;\n')
	
	print 'DONE!'

	# run the plotting script
	print 'Generating PNGs...'
	os.system('gnuplot %s' % gnuplot_file)
	
	print 'DONE!'

if __name__ == '__main__':
	main()
