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
	parser.add_argument("--folder", dest="caseDir", help="folder in which the cases for different mesh sizes are present", default=os.path.expandvars("${CUIBM_DIR}/cases/convergence/cavityRe100/NavierStokes/20x20"))
	parser.add_argument("--run-cases", dest="runCases", help="run the cases if this flag is used", action='store_true', default=False)
	args = parser.parse_args()

	# list of folders from which velocity data is to be obtained
	folders = sorted(os.walk(args.caseDir).next()[1], key=int)
	numFolders = len(folders)

	# run the cases in each of the folders
	if args.runCases:
		for folder in folders:
			runCommand = [os.path.expandvars("${CUIBM_DIR}/bin/cuIBM"), '-caseFolder', "%s/%s" % (args.caseDir, folder)]
			print " ".join(runCommand)
			subprocess.call(runCommand)

	# create arrays to store the required values
	U  = []
	errNorm  = np.zeros(numFolders-1)
	meshSize = np.zeros(numFolders-1)

	stride = 1
	for fIdx, folder in enumerate(folders):
		# path to folder
		folderPath = os.path.expandvars("%s/%s" % (args.caseDir, folder));
		# read simulation information
		nt, _, _, _ = readSimulationParameters(folderPath)

		# read the grid data
		# nx and ny are the number of cells
		# dx and dy are the cell widths
		# xu and yu are the coordinates of the locations where U is calculated
		nx, ny, dx, dy, xu, yu, _, _ = readGridData(folderPath)
		
		if fIdx>0:
			meshSize[fIdx-1] = nx

		# read velocity data
		u, _ = readVelocityData(folderPath, nt, nx, ny, dx, dy)

		U.append(np.reshape(u, (ny, nx-1))[stride/2::stride,stride-1::stride])
		print U[fIdx].shape
		
		print 'Completed folder %s' % folder
		stride = stride*3

	for idx in range(numFolders-1):
		errNorm[idx] = la.norm(U[idx+1]-U[idx])
	
	orderOfConvergence = -np.polyfit(np.log10(meshSize), np.log10(errNorm), 1)[0]
	print meshSize
	print errNorm
	print np.log(errNorm[0]/errNorm[1])/np.log(3), np.log(errNorm[1]/errNorm[2])/np.log(3)
	print orderOfConvergence
	
	plt.loglog(meshSize, errNorm, 'o-b', label="L-2 norm of difference\nOrder of convergence=%1.3f" % orderOfConvergence)
	plt.axis([1, 1e4, 1e-4, 10])
	x  = np.linspace(1, 1e4, 2)
	x1 = 1/x
	x2 = 1/x**2
	plt.loglog(x, x1, '--k', label="First-order convergence")
	plt.loglog(x, x2, ':k', label="Second-order convergence")
	plt.legend()
	ax = plt.axes()
	ax.set_xlabel("Mesh size")
	ax.set_ylabel("L-2 Norm of difference between solutions on consecutive grids")
	plt.savefig("%s/convergence.png" % args.caseDir)

if __name__ == "__main__":
	main()