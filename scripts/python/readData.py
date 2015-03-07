#!/usr/bin/env python

# file: $CUIBM_DIR/scripts/python/readData.py
# author: Anush Krishnan (anush@bu.edu), Olivier Mesnard (mesnardo@gwu.edu)
# description: Definitions of functions to read numerical data files.


import sys
import struct
import argparse

import numpy


def read_simulation_parameters(folder_path):
	"""Reads the parameters of the simulation.
	
	Arguments
	---------
	folder_path -- directory of the simulation.

	Returns
	-------
	nt -- total number of time-steps.
	start_step -- initial time-step.
	nsave -- saving interval.
	dt -- time-increment.
	"""
	# open and read the run.info file
	run_file = '%s/run.info' % folder_path
	with open(run_file, 'r') as infile:
		parameters = infile.read().split()
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
	args = parser.parse_args(parameters)
	return args.nt, args.start_step, args.nsave, args.dt


def read_grid(folder_path):
	"""Reads the computational grid from the binary grid file.
	
	Arguments
	---------
	folder_path -- directory of the simulation.

	Returns
	-------
	nx, ny -- number of cells in the x- and y- directions.
	dx, dy -- cell-widths in the x- and y- directions.
	xu, yu -- coordinates of the u-velocity points.
	xv, yv -- coordinates of the v-velocity points.
	"""
	grid_file = '%s/grid' % folder_path
	with open(grid_file, 'rb') as infile:
		# x-direction
		nx = struct.unpack('i', infile.read(4))[0]
		x = numpy.array(struct.unpack('d'*(nx+1), infile.read(8*(nx+1))))
		dx = x[1:] - x[:-1]
		xu = x[1:-1]
		xv = 0.5 * (x[1:] + x[:-1])
		# y-direction
		ny = struct.unpack('i', infile.read(4))[0]
		y = numpy.array(struct.unpack('d'*(ny+1), infile.read(8*(ny+1))))
		dy = y[1:] - y[:-1]
		yu = 0.5 * (y[1:] + y[:-1])
		yv = y[1:-1]
	return nx, ny, dx, dy, xu, yu, xv, yv


def read_velocity(folder_path, time_step, nx, ny, dx, dy):
	"""Reads velocity data from binary flux file at a given time-step.
	
	Arguments
	---------
	folder_path -- directory of the simulation.
	time_step -- time-step to read the data.
	nx, ny -- number of cells in the x- and y- directions.
	dx, dy -- cell-widths in the x- and y- directions.

	Returns
	-------
	u, v -- u-velocity values and v-velocity values.
	"""
	flux_file = '%s/%07d/q' % (folder_path, time_step)
	with open(flux_file, 'rb') as infile:
		nq = struct.unpack('i', infile.read(4))[0]
		q = numpy.array(struct.unpack('d'*nq, infile.read(8*nq)))
	# u-velocities
	u = numpy.empty((nx-1)*ny, dtype=float)
	for j in xrange(ny):
		for i in xrange(nx-1):
			u[j*(nx-1)+i] = q[j*(nx-1)+i] / dy[j]
	# v-velocities
	offset = u.size
	v = numpy.empty(nx*(ny-1), dtype=float)
	for j in xrange(ny-1):
		for i in xrange(nx):
			v[j*nx+i] = q[offset+j*nx+i] / dx[i]
	return u, v


def read_pressure(folder_path, time_step, nx, ny):
	"""Reads pressure data from binary lambda file at given time-step.
	
	Arguments
	---------
	folder_path -- directory of the simulation.
	time_step -- time-step to read the data.
	nx, ny -- number of cells in the x- and y- directions.

	Returns
	-------
	p -- pressure values.
	"""
	lambda_file = '%s/%07d/lambda' % (folder_path, time_step)
	with open(lambda_file, 'rb') as infile:
		nlambda = struct.unpack('i', infile.read(4))[0]
		lamBda = numpy.array(struct.unpack('d'*nlambda, infile.read(8*nlambda)))
	# return only the pressure data
	return lamBda[:nx*ny]


def read_mask(folder_path, nx, ny):
	"""Reads the mask file.
	
	Arguments
	---------
	folder_path -- directory of the simulation.
	nx, ny -- number of cells in the x- and y- directions.

	Returns
	-------
	mask -- the mask in the x- and y- directions.
	"""
	mask_file = '{}/mask.txt'.format(folder_path)
	mask = numpy.loadtxt(mask_file)
	offset = (nx-1) * ny
	return mask[:offset], mask[offset:]


def main():
	pass


if __name__ == '__main__':
	main()