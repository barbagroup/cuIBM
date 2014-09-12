#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/readData.py
# author: Krishnan, A. (anush@bu.edu)
# description: definitions of functions to read numerical data


import argparse

import numpy as np


def readSimulationParameters(folder):
	"""Reads the parameters of the simulation.
	
	Arguments
	---------
	folder -- path of the simulation case.

	Returns
	-------
	nt -- total number of time-steps.
	start_step -- initial time-step.
	nsave -- saving interval.
	dt -- time-increment.
	"""
	# open and read the info file
	with open(folder+'/run.info', 'r') as infile:
		option_list = infile.read().split()
	
	# create the parser
	parser = argparse.ArgumentParser()
	# fill the parser with arguments
	parser.add_argument('--nx', dest='nx', type=int,
						help='number of cells in the x-direction')
	parser.add_argument('--ny', dest='ny', type=int,
						help='number of cells in the y-direction')
	parser.add_argument('--nt', dest='nt', type=int,
						help='total number of time-steps')
	parser.add_argument('--nsave', dest='nsave', type=int,
						help='saving interval')
	parser.add_argument('--startStep', dest='start_step', type=int,
						help='starting time-step')
	parser.add_argument('--dt', dest='dt', type=float,
						help='time-increment')
	parser.add_argument('--vortlim', dest='vort_lim', type=float,
						help='vorticity limit')
	parser.add_argument('--folder', dest='folder_name', type=str,
						help='name of the folder in which data is saved')
	parser.add_argument('--nu', dest='nu', type=float,
						help='kinematic viscosity')
	parser.add_argument('--input', dest='body_file', type=str,
						help='body input file')
	parser.add_argument('--flowFile', dest='flow_file', type=str,
						help='flow description file')
	parser.add_argument('--domainFile', dest='domain_file', type=str,
						help='mesh description file')
	parser.add_argument('--simulation_file', dest='simulation_file', type=str,
						help='simulation parameters file')
	parser.add_argument('--bodyFile', dest='body_file', type=str,
						help='body file')
	# parse the list
	args = parser.parse_args(option_list)

	return args.nt, args.start_step, args.nsave, args.dt


def readGridData(folder):
	# open the file with the grid data
	gridFile = folder + '/grid'
	f = open(gridFile, 'r')

	# read the number of cells in the x-direction
	a = f.readline().strip().split()
	nx = int(a[0])
	
	# arrays to store grid information
	x = np.zeros(nx+1)
	dx = np.zeros(nx)
	xu = np.zeros(nx-1)
	xv = np.zeros(nx)

	# read the x-coordinates of the grid nodes
	a = f.readline().strip().split()
	x[0] = float(a[0])
	for i in range(nx):
		a = f.readline().strip().split()
		if a==[]:
			break
		x[i+1] = float(a[0])
		dx[i] = x[i+1] - x[i]
		xv[i] = 0.5*(x[i+1] + x[i]);
		if i < nx-1:
			xu[i] = x[i+1]

	a = f.readline() # skip the empty lines

	# read the number of cells in the y-direction
	a = f.readline().strip().split()
	ny = int(a[0])
	
	# arrays to store y-grid information
	y = np.zeros(ny+1)
	dy = np.zeros(ny)
	yu = np.zeros(ny)
	yv = np.zeros(ny-1)

	# read the y-coordinates of the grid nodes
	a = f.readline().strip().split()
	y[0] = float(a[0])
	for j in range(ny):
		a = f.readline().strip().split()
		if a==[]:
			break
		y[j+1] = float(a[0])
		dy[j] = y[j+1] - y[j]
		yu[j] = 0.5*(y[j+1] + y[j])
		if j < ny-1:
			yv[j] = y[j+1]
	
	f.close()
	
	return nx, ny, dx, dy, xu, yu, xv, yv


def readVelocityData(folder, time_step, nx, ny, dx, dy):
	"""Reads the velocity data at a given time-step.
	
	Arguments
	---------
	folder -- path of the simulation case.
	time-step -- current time-step.
	nx, ny -- number of cells on the x- and y- directions.
	dx, dy -- cell widths in the x- and y- directions.
	"""
	flux_file = '%s/%07d/q' % (folder, time_step)
	with open(flux_file, 'r') as infile:
		nq = int(infile.readline())			# length of flux vector q
		q = np.loadtxt(infile, dtype=float)	# store flux vector q
	
	# store u-velocities
	n_u = (nx-1) * ny				# number of u-velocity points
	u = np.empty(n_u, dtype=float)
	for j in xrange(ny):
		for i in xrange(nx-1):
			u[j*(nx-1)+i] = q[j*(nx-1)+i]/dy[j]

	# store v-velocities
	v = np.empty(nx*(ny-1), dtype=float)
	for j in xrange(ny-1):
		for i in xrange(nx):
			v[j*nx+i] = q[n_u+j*nx+i]/dx[i]
	
	return u, v
