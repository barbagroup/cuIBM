#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

def interpolateFromStructuredGrid(X, Y, u, x, nx, y, ny):
	i=0;
	j=0;
	while X > x[i+1]:
		i = i+1
	while Y > y[j+1]:
		j = j+1
	xi  = (X-x[i])/(x[i+1]-x[i]);
	eta = (Y-y[j])/(y[j+1]-y[j]);
	idx = j*nx+i
	uInterp = u[idx]*(1.-xi)*(1.-eta) + u[idx+1]*xi*(1.-eta) + u[idx+nx+1]*eta*xi + u[idx+nx]*(1.-xi)*eta;
	return uInterp

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
	#xv = np.zeros(nx)

	# read the x-coordinates of the grid nodes
	a = f.readline().strip().split()
	x[0] = float(a[0])
	for i in range(nx):
		a = f.readline().strip().split()
		if a==[]:
			break
		x[i+1] = float(a[0])
		dx[i] = x[i+1] - x[i]
		#xv[i] = 0.5*(x[i+1] + x[i]);
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
	#yv = np.zeros(ny-1)

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
		#if j < ny-1:
		#	yv[j] = y[j+1]

	f.close()
	
	print "%d, %d" % (nx, ny)
	return nx, ny, dx, dy, xu, yu

def readSimulationParameters(infoFile):
	# create options list
	f = open(infoFile, 'r')
	optionsList = f.read().split()
	f.close()
	
	# Command line options used in simulation
	parser = argparse.ArgumentParser()
	# create list of options
	parser.add_argument("--nx", type=int,dest="nx", help="number of cells in the x-direction")
	parser.add_argument("--ny", type=int,dest="ny", help="number of cells in the y-direction")
	parser.add_argument("--nt", type=int,dest="nt", help="total number of time steps")
	parser.add_argument("--nsave", type=int, dest="nsave", help="save interval")
	parser.add_argument("--startStep", type=int, dest="startStep", help="start time step")
	parser.add_argument("--dt", type=float, dest="dt", help="time step")
	parser.add_argument("--vortlim", type=float,dest="vortlim", help="vorticity limit")
	parser.add_argument("--folder", dest="foldername", help="name of folder in which data is saved")
	parser.add_argument("--nu", type=float, dest="nu", help="kinematic viscosity")
	parser.add_argument("--input", dest="input", help="Body input file")
	parser.add_argument("--flowFile", dest="flowFile", help="Flow description file")
	parser.add_argument("--domainFile", dest="domainFile", help="Mesh description file")
	parser.add_argument("--simulationFile", dest="simulationFile", help="Simulation parameters file")
	parser.add_argument("--bodyFile", dest="bodyFile", help="body file")
	# parse options list
	args = parser.parse_args(optionsList)
	
	return args.nsave, args.nt, args.startStep

if __name__ == "__main__":
	# list of folders from which velocity data is to be obtained
	folders = ['400', '600', '900', '1350']
	numFolders = len(folders)

	# for each of the above case folders with different mesh sizes
	# simulate flow over a cylinder at Re=40
	for case in folders:
		runCommand = os.path.expandvars("${CUIBM_DIR}/bin/cuIBM -caseFolder ${CUIBM_DIR}/orderOfConvergence/cylinderRe40/%s" % case)
		os.system(runCommand)

	# define the grid on which the values are interpolated
	# it is uniform and of size 41 x 41 and encompasses the region [-1,1]x[-1,1]
	Ni = 41

	# create arrays to store the interpolated values
	Ui = np.zeros((numFolders, Ni*Ni))
	Xi = np.linspace(-1, 1, Ni)
	Yi = np.linspace(-1, 1, Ni)
	errNorm  = np.zeros(numFolders-1)
	meshSize = np.zeros(numFolders-1)

	for fIdx in range(numFolders):
		# path to folder
		folder = os.path.expandvars("${CUIBM_DIR}/orderOfConvergence/cylinderRe40/%s"%(folders[fIdx]));
		# path to info file
		infoFile = folder + '/run.info'
		# read simulation information
		nsave, nt, startStep = readSimulationParameters(infoFile)

		# read the grid data
		# nx and ny are the number of cells
		# dx and dy are the cell widths
		# xu and yu are the coordinates of the locations where U is calculated
		nx, ny, dx, dy, xu, yu = readGridData(folder)
		
		if fIdx>0:
			meshSize[fIdx-1] = nx

		numU = (nx-1)*ny    # number of elements in U 
		u = np.zeros(numU)  # vector to store U
	
		timeStep = nt
		fluxFile = folder + '/' + ("%07d" % timeStep) + "/q"
		f = open(fluxFile, 'r')

		# length of vector q
		a = f.readline().strip().split()
		nq = int(a[0])

		# read the u velocities
		for j in range(ny):
			for i in range(nx-1):
				a = f.readline().strip().split()
				u[j*(nx-1)+i] = float(a[0])/dy[j]
	
		idx = 0
		for Y in Yi:
			for X in Xi:
				Ui[fIdx, idx] = interpolateFromStructuredGrid(X, Y, u, xu, nx-1, yu, ny)
				idx = idx+1

		#print Ui[fIdx]
		print 'Completed folder %s' % folder

	for idx in range(numFolders-1):
		errNorm[idx] = la.norm(Ui[idx+1]-Ui[idx])

	print meshSize
	print errNorm

	outFile = os.path.expandvars("${CUIBM_DIR}/orderOfConvergence/cylinderRe40/errNorms.txt");
	f = open(outFile, 'w')
	for c1, c2 in zip(meshSize, errNorm):
		f.write("%d\t%f\n" % (c1, c2))
	f.close()

	runCommand = os.path.expandvars("gnuplot ${CUIBM_DIR}/scripts/convergence/cylinderRe40.plt")
	os.system(runCommand)