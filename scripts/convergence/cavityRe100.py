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
from readData import readSimulationParameters, readGridData, readVelocityData, readMask
import subprocess
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def line():
	print '-'*80

def main():
	# Command line options
	parser = argparse.ArgumentParser(description="Calculates the order of convergence.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-folder", dest="caseDir", help="folder in which the cases for different mesh sizes are present", default=os.path.expandvars("${CUIBM_DIR}/cases/convergence/cavityRe100/NavierStokes/20x20"))
	parser.add_argument("-tolerance", dest="tolerance", help="folder in which the cases for different mesh sizes are present", default="")
	parser.add_argument("-interpolation_type", dest="interpolation_type", help="the type of interpolation used in the direct forcing method", default="")
	parser.add_argument("-run_simulations", dest="runSimulations", help="run the cases if this flag is used", action='store_true', default=False)
	parser.add_argument("-remove_mask", dest="removeMask", help="use values from the entire domain", action='store_true', default=False)
	args = parser.parse_args()

	commandOptions = []
	if args.tolerance:
		commandOptions.extend(["-velocityTol", "{}".format(args.tolerance), "-poissonTol", "{}".format(args.tolerance)])
	if args.interpolation_type:
		commandOptions.extend(["-interpolationType", "{}".format(args.interpolation_type)])

	# list of folders from which velocity data is to be obtained
	folders = sorted(os.walk(args.caseDir).next()[1], key=int)
	numFolders = len(folders)

	# run the cases in each of the folders
	if args.runSimulations:
		for folder in folders:
			runCommand = [os.path.expandvars("${CUIBM_DIR}/bin/cuIBM"),
							'-caseFolder', "{}/{}".format(args.caseDir, folder)]
			runCommand.extend(commandOptions)
			print " ".join(runCommand)
			subprocess.call(runCommand)

	# create arrays to store the required values
	U = []
	V = []
	errL2NormU  = np.zeros(numFolders-1)
	errL2NormV  = np.zeros(numFolders-1)
	errLInfNormU  = np.zeros(numFolders-1)
	errLInfNormV  = np.zeros(numFolders-1)
	meshSize = np.zeros(numFolders-1, dtype=int)

	line()

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
		
		if fIdx==0:
			initialMeshSpacing = dx[0]
		else:
			meshSize[fIdx-1] = nx

		# read velocity data
		u, v = readVelocityData(folderPath, nt, nx, ny, dx, dy)

		if not args.removeMask:
			# read mask
			mask_u, mask_v = readMask(folderPath, nx, ny)
			u[:] = u[:]*mask_u[:]
			v[:] = v[:]*mask_v[:]

		U.append(np.reshape(u, (ny, nx-1))[stride/2::stride,stride-1::stride])
		V.append(np.reshape(v, (ny-1, nx))[stride-1::stride,stride/2::stride])
		
		print 'Completed folder {}. u:{}, v:{}'.format(folder, U[fIdx].shape, V[fIdx].shape)
		stride = stride*3

	for idx in range(numFolders-1):
		errL2NormU[idx] = la.norm(U[idx+1]-U[idx])
		errL2NormV[idx] = la.norm(V[idx+1]-V[idx])

		errLInfNormU[idx] = la.norm(U[idx+1]-U[idx], np.inf)
		errLInfNormV[idx] = la.norm(V[idx+1]-V[idx], np.inf)
		
		if idx==0:
			h = initialMeshSpacing
			x = np.arange(h/2., 1., h)
			y = np.arange(h, 1., h)
			X, Y = np.meshgrid(x,y)
			plt.ioff()
			fig = plt.figure()
			ax = fig.add_subplot(111)
			diffV = np.abs(V[idx+1]-V[idx] )
			CS = ax.pcolor(X, Y, diffV, norm=LogNorm(vmin=1e-10, vmax=1))
			fig.gca().set_aspect('equal', adjustable='box')
			fig.colorbar(CS)
			if args.removeMask:
				fig.savefig("{}/diff_nomask.png".format(args.caseDir))
			else:
				fig.savefig("{}/diff.png".format(args.caseDir))
	
	L2OrderOfConvergenceU = -np.polyfit(np.log10(meshSize), np.log10(errL2NormU), 1)[0]
	L2OrderOfConvergenceV = -np.polyfit(np.log10(meshSize), np.log10(errL2NormV), 1)[0]

	LInfOrderOfConvergenceU = -np.polyfit(np.log10(meshSize), np.log10(errLInfNormU), 1)[0]
	LInfOrderOfConvergenceV = -np.polyfit(np.log10(meshSize), np.log10(errLInfNormV), 1)[0]
	
	line()
	print "Mesh sizes: {}".format(meshSize)
	line()
	print "U:"
	print "errL2Norms: {}".format(errL2NormU)
	print "Convergence rates: {:.3f}, {:.3f}".format(np.log(errL2NormU[0]/errL2NormU[1])/np.log(3), np.log(errL2NormU[1]/errL2NormU[2])/np.log(3))
	print "Linear fit convergence rate: {:.3f}".format(L2OrderOfConvergenceU)

	print "\nV:"
	print "errL2Norms: {}".format(errL2NormV)
	print "Convergence rates: {:.3f}, {:.3f}".format(np.log(errL2NormV[0]/errL2NormV[1])/np.log(3), np.log(errL2NormV[1]/errL2NormV[2])/np.log(3))
	print "Linear fit convergence rate: {:.3f}".format(L2OrderOfConvergenceV)
	
	plt.clf()
	plt.loglog(meshSize, errL2NormU, 'o-b', label="L-2 norm of difference in $u$\nOrder of convergence={:.3f}".format(L2OrderOfConvergenceU))
	plt.loglog(meshSize, errL2NormV, 'o-r', label="L-2 norm of difference in $v$\nOrder of convergence={:.3f}".format(L2OrderOfConvergenceV))
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
	plt.savefig("{}/L2Convergence.png".format(args.caseDir))

	line()
	print "U:"
	print "errLInfNorms: {}".format(errLInfNormU)
	print "Convergence rates: {:.3f}, {:.3f}".format(np.log(errLInfNormU[0]/errLInfNormU[1])/np.log(3), np.log(errLInfNormU[1]/errLInfNormU[2])/np.log(3))
	print "Linear fit convergence rate: {:.3f}".format(LInfOrderOfConvergenceU)

	print "\nV:"
	print "errLInfNorms: {}".format(errLInfNormV)
	print "Convergence rates: {:.3f}, {:.3f}".format(np.log(errLInfNormV[0]/errLInfNormV[1])/np.log(3), np.log(errLInfNormV[1]/errLInfNormV[2])/np.log(3))
	print "Linear fit convergence rate: {:.3f}".format(LInfOrderOfConvergenceV)
	line()

	plt.clf()
	plt.loglog(meshSize, errLInfNormU, 'o-b', label="L-Inf norm of difference in $u$\nOrder of convergence={:.3f}".format(LInfOrderOfConvergenceU))
	plt.loglog(meshSize, errLInfNormV, 'o-r', label="L-Inf norm of difference in $v$\nOrder of convergence={:.3f}".format(LInfOrderOfConvergenceV))
	plt.axis([1, 1e4, 1e-4, 100])
	x  = np.linspace(1, 1e4, 2)
	x1 = 1/x
	x2 = 1/x**2
	plt.loglog(x, x1, '--k', label="First-order convergence")
	plt.loglog(x, x2, ':k', label="Second-order convergence")
	plt.legend()
	ax = plt.axes()
	ax.set_xlabel("Mesh size")
	ax.set_ylabel("L-Inf Norm of difference between solutions on consecutive grids")
	plt.savefig("{}/LInfConvergence.png".format(args.caseDir))

if __name__ == "__main__":
	main()