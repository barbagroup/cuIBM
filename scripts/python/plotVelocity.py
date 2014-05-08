#!/usr/bin/env python
import os
import sys
import argparse
import numpy as np
import readData as rd
import matplotlib.pyplot as plt

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
	nx, ny, dx, dy, xu, yu, xv, yv = rd.readGridData(folder)

	Xu, Yu = np.meshgrid(xu, yu)
	Xv, Yv = np.meshgrid(xv, yv)
	

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

		u = np.reshape(u, (ny, nx-1))
		plt.contour(Xu, Yu, u, levels=np.arange(-1., 1., 0.1))
		plt.axis('equal')
		plt.axis([args.bl_x, args.tr_x, args.bl_y, args.tr_y])
		plt.colorbar()
		plt.savefig("%s/u%07d.png" % (folder,k))
		plt.clf()

		v = np.reshape(v, (ny-1, nx))
		plt.contour(Xv, Yv, v, levels=np.arange(-1., 1., 0.1))
		plt.axis('equal')
		plt.axis([args.bl_x, args.tr_x, args.bl_y, args.tr_y])
		plt.colorbar()
		plt.savefig("%s/v%07d.png" % (folder,k))
		plt.clf()

		k = k + nsave
