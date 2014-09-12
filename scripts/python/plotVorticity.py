#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/plotVorticity.py
# author: Krishnan, A. (anush@bu.edu)
# description: plot the contour of vorticity at every time saved


import os
import argparse

import numpy as np

from readData import readSimulationParameters, readGridData


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
	args = read_inputs()

	folder = args.folder	# name of the folder

	# read the parameters of the simulation
	nt, start_step, nsave, _ = readSimulationParameters(folder)

	if not args.plot_only:
		# calculate the mesh characteristics
		nx, ny, dx, dy, _, yu, xv, _ = readGridData(folder)

		# calculate appropriate array boundaries
		i_start = np.where(xv >= args.bl_x)[0][0]
		i_end = np.where(xv <= args.tr_x)[0][0]
		j_start = np.where(yu >= args.bl_y)[0][0]
		j_end = np.where(yu <= args.tr_y)[0][0]

		ite = start_step + nsave
		while ite < nt+1:
			# read the velocity data at the given time-step
			u, v = readVelocityData(folder, ite, nx, ny, dx, dy)
			if u == None or v == None:
				break

			

			vort_file = '%s/%07d/vorticity' % (folder, ite)
			print vort_file
			

if __name__ == '__main__':
	main()

	k = startStep + nsave
	while k < (nt+1):
		u, v = rd.readVelocityData(folder, k, nx, ny, dx, dy)
		if u==None or v==None:
			break
		vortFile = args.folder + ("/%07d" % k) + "/vorticity"
		print vortFile
		f = open(vortFile, 'w')
		for j in range(Jstart, Jend):
			y  = (yu[j] + yu[j+1])/2.
			Dy = (dy[j] + dy[j+1])/2.
			for i in range(Istart, Iend):
				# calculate the vorticity
				# w = dv/dx-du/dy
				x  = (xv[i] + xv[i+1])/2.
				Dx = (dx[i] + dx[i+1])/2.
				vort = (v[j*nx+i+1] - v[j*nx+i])/Dx - (u[(j+1)*(nx-1)+i] - u[j*(nx-1)+i])/Dy
				f.write("%f\t%f\t%f\n" % (x, y, vort))
			f.write("\n")
		f.close()
		k = k + nsave

vortFile = "%s/vorticity.plt" % folder
f = open( vortFile, 'w')

print 'Creating gnuplot file... ';

f.write("reset;\n")
f.write("""set terminal pngcairo enhanced font "Times, 15" size 900,600;\n""")

k = startStep+nsave
while k < (nt+1):
	# image file name
	f.write("\nset output '%s/plot%07d.png'\n" % (folder, k))
	f.write("set multiplot;\n")
	# vorticity
	f.write("reset;\n")
	f.write("set view map; set size ratio -1; unset key;\n")
	f.write("set pm3d map;\n")
	f.write("set palette defined (-2 'dark-blue', -1 'light-blue', 0 'white', 1 'light-red', 2 'dark-red');\n")
	f.write("set cbrange [%f:%f];\n" % (-args.vortlim, args.vortlim))
	f.write("splot [%f:%f] [%f:%f] '%s/%07d/vorticity';\n" % (args.bl_x, args.tr_x, args.bl_y, args.tr_y, folder, k) )
	# bodies
	#f.write("reset;\n")
	#f.write("set view map; set size ratio -1; unset key;\n")
	#f.write("splot [%f:%f] [%f:%f] '%s/%07d/bodies' u 1:2:(0) w l lw 2 lc 'black';\n" % (args.bl_x, args.tr_x, args.bl_y, args.tr_y, folder, i) )
	f.write("unset multiplot;\n")
	k = k + nsave
f.close()

print 'DONE!'

# Run the plotting script
print 'Generating PNGs... '
os.system("gnuplot %s" % vortFile)
print 'DONE!'
