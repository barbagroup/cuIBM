#!/usr/bin/env python

import os
import argparse
import numpy as np
import numpy.linalg as la
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.expandvars("${CUIBM_DIR}/scripts/python"))
from readData import readSimulationParameters, readGridData, readVelocityData
import subprocess

def main():
	# Command line options
	parser = argparse.ArgumentParser(description="Calculates the order of convergence.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-folder", dest="caseDir", help="folder in which the cases for different mesh sizes are present", default=os.path.expandvars("${CUIBM_DIR}/cases/convergence/cavityRe100/NavierStokes/20x20"))
	parser.add_argument("-tolerance", dest="tolerance", help="folder in which the cases for different mesh sizes are present", default=1.e-8)
	parser.add_argument("-run_simulations", dest="runSimulations", help="run the cases if this flag is used", action='store_true', default=False)
	args = parser.parse_args()

	# list of folders from which velocity data is to be obtained
	folders = sorted(os.walk(args.caseDir).next()[1], key=int)
	numFolders = len(folders)

	# run the cases in each of the folders
	if args.runSimulations:
		for folder in folders:
			runCommand = [os.path.expandvars("${CUIBM_DIR}/bin/cuIBM"),
							'-caseFolder', "{}/{}".format(args.caseDir, folder),
							'-velocityTol', "{}".format(args.tolerance),
							'-poissonTol', "{}".format(args.tolerance)]
			print " ".join(runCommand)
			subprocess.call(runCommand)

	# create arrays to store the required values
	U = []
	V = []
	errNormU  = np.zeros(numFolders-1)
	errNormV  = np.zeros(numFolders-1)
	meshSize = np.zeros(numFolders-1, dtype=int)

	stride = 1
	for fIdx, folder in enumerate(folders):
		# path to folder
		folderPath = os.path.expandvars("{}/{}".format(args.caseDir, folder));
		# read simulation information
		nt, _, _, _ = readSimulationParameters(folderPath)

		# read the grid data
		# nx and ny are the number of cells
		# dx and dy are the cell widths
		# xu and yu are the coordinates of the locations where U is calculated
		nx, ny, dx, dy, xu, yu, xv, yv = readGridData(folderPath)
		
		if fIdx>0:
			meshSize[fIdx-1] = nx

		# read velocity data
		u, v = readVelocityData(folderPath, nt, nx, ny, dx, dy)

		U.append(np.reshape(u, (ny, nx-1))[stride/2::stride,stride-1::stride])
		V.append(np.reshape(v, (ny-1, nx))[stride-1::stride,stride/2::stride])
		print U[fIdx].shape, V[fIdx].shape
		
		print 'Completed folder {}'.format(folder)
		stride = stride*3

	for idx in range(numFolders-1):
		errNormU[idx] = la.norm(U[idx+1]-U[idx])
		errNormV[idx] = la.norm(V[idx+1]-V[idx])
	
	orderOfConvergenceU = -np.polyfit(np.log10(meshSize), np.log10(errNormU), 1)[0]
	orderOfConvergenceV = -np.polyfit(np.log10(meshSize), np.log10(errNormV), 1)[0]
	
	print "\nMesh sizes: {}".format(meshSize)
	
	print "\nU:"
	print "errNorms: {}".format(errNormU)
	print "Convergence rates: {:.3f}, {:.3f}".format(np.log(errNormU[0]/errNormU[1])/np.log(3), np.log(errNormU[1]/errNormU[2])/np.log(3))
	print "Linear fit convergence rate: {:.3f}".format(orderOfConvergenceU)

	print "\nV:"
	print "errNorms: {}".format(errNormV)
	print "Convergence rates: {:.3f}, {:.3f}".format(np.log(errNormV[0]/errNormV[1])/np.log(3), np.log(errNormV[1]/errNormV[2])/np.log(3))
	print "Linear fit convergence rate: {:.3f}".format(orderOfConvergenceV)
	
	
	plt.loglog(meshSize, errNormU, 'o-b', label="L-2 norm of difference in $u$\nOrder of convergence={:.3f}".format(orderOfConvergenceU))
	plt.loglog(meshSize, errNormV, 'o-r', label="L-2 norm of difference in $v$\nOrder of convergence={:.3f}".format(orderOfConvergenceV))
	plt.axis([1, 1e4, 1e-4, 100])
	x  = np.linspace(1, 1e4, 2)
	x1 = 1/x
	x2 = 1/x**2
	plt.loglog(x, x1, '--k', label="First-order convergence")
	plt.loglog(x, x2, ':k', label="Second-order convergence")
	plt.legend()
	ax = plt.axes()
	ax.set_xlabel("Mesh size")
	ax.set_ylabel("L-2 Norm of difference between solutions on consecutive grids")
	plt.savefig("{}/convergence.png".format(args.caseDir))

if __name__ == "__main__":
	main()