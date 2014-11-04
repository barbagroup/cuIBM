#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/analyzeVorticity.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: script to analyse the vorticity


import argparse
import math

import numpy as np
from matplotlib import pyplot as plt

from readData import readSimulationParameters, readGridData


def read_inputs():
	"""Parse the command-line"""
	# create the parser
	parser = argparse.ArgumentParser()
	# fill the parser
	parser.add_argument('--folder', dest='folder', type=str, default='.',
						help='path of the simualtion case')
	parser.add_argument('--show', dest='show', action='store_true',
						help='boolean to plot the results in a figure')
	return parser.parse_args()


def main():
	"""Analyses the vorticity given a case."""
	# parse the command-line
	args = read_inputs()

	folder = args.folder	# path of the case

	# read parameters of the simulation
	nt, start_step, nsave, dt = readSimulationParameters(folder)

	# read and generate the computational mesh
	nx, ny, dx, dy, _, yu, xv, _ = readGridData(folder)
	
	# number of saved points
	save_points = int( (nt-start_step)/nsave )

	T = 10.0	# total time of the simulation
	h = dx[0]	# uniform cell-width

	# initialization
	total_vorticity = np.empty(save_points, dtype=float)
	circulation = np.empty(save_points, dtype=float)
	time = np.empty(save_points, dtype=float)

	# time-loop
	ite = start_step + nsave
	idx = 0
	while ite < nt+1:
		time[idx] = ite*dt
		# read the vorticity
		vort_file = '%s/%07d/vorticity' % (folder, ite)
		with open(vort_file, 'r') as infile:
			vorticity = np.loadtxt(infile, dtype=float, usecols=(2,))
		# calculate the total vorticity
		total_vorticity[idx] = h**2*vorticity.sum()
		# calculate the theoretical circulation
		circulation[idx] = -math.sin(math.pi*time[idx]/T)
		idx += 1
		ite += nsave

	# plot the total vorticity and theretical circulation
	plt.figure()
	plt.grid(True)
	plt.xlabel('Time', fontsize=16)
	plt.ylabel('Total circulation', fontsize=16)
	plt.plot(time, total_vorticity, label='Total vorticity',
			 color='r', ls='.', marker='o', markersize=10)
	plt.plot(time, circulation, label='Theoretical circulation',
			 color='b', ls='-', lw=2)
	plt.xlim(0, 80)
	plt.ylim(-1.5, 1.5)
	plt.title('Comparison of Total Vorticity and Theoretical Circulation '
			  'at Re=10,000')
	plt.legend(loc='best', prop={'size': 16})
	plt.savefig('%s/total_vorticity.png' % folder)
	if args.show:
		plt.show()


if __name__ == '__main__':
	main()
