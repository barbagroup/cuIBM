#!/usr/bin/env python
import os
import sys
import argparse
import numpy as np
import readData as rd

# Generate optionsfile as per command line options
parser = argparse.ArgumentParser(description="Plots the vorticity field at all save points for a specified simulation", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--folder", dest="folder", help="Case folder", default=".")
parser.add_argument("--blx", type=float, dest="bl_x", help="x-coordinate of the bottom left of the plot region", default=-2)
parser.add_argument("--bly", type=float, dest="bl_y", help="y-coordinate of the bottom left of the plot region", default=-3)
parser.add_argument("--trx", type=float, dest="tr_x", help="x-coordinate of the top right of the plot region", default=4)
parser.add_argument("--try", type=float, dest="tr_y", help="y-coordinate of the top right of the plot region", default=3)
parser.add_argument("--vortlim", type=float, dest="vortlim", help="cutoff vorticity on the plot", default=15)
parser.add_argument("--plot-only", dest="plot_only", help="only generate the vorticity plots and do no recalculate the data", action='store_true')
parser.set_defaults(plot_only=False)
args = parser.parse_args()

folder = args.folder
nt, startStep, nsave, _ = rd.readSimulationParameters(folder)

if not args.plot_only:
	nx, ny, dx, dy, _, yu, xv, _ = rd.readGridData(folder)

	Istart = 0
	while args.bl_x > xv[Istart]:
		Istart = Istart+1
	Iend = nx-1
	while args.tr_x < xv[Iend]:
		Iend = Iend-1

	Jstart = 0
	while args.bl_y > yu[Jstart]:
		Jstart = Jstart+1
	Jend = ny-1
	while args.tr_y < yu[Jend]:
		Jend = Jend-1

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
