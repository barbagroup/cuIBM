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
from matplotlib.font_manager import FontProperties
import matplotlib.font_manager as font_manager
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

def line():
	print '-'*80

def main():
	
	#view available fonts
	#for font in font_manager.findSystemFonts():
	#	print font

	try:
		fontpath = '/usr/share/fonts/truetype/msttcorefonts/Georgia.ttf'

		prop = font_manager.FontProperties(fname=fontpath)
		matplotlib.rcParams['font.family'] = prop.get_name()
	except:
		pass
	# Command line options
	parser = argparse.ArgumentParser(description="Calculates the order of convergence.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-folder", dest="caseDir", help="folder in which the cases for different mesh sizes are present", default=os.path.expandvars("${CUIBM_DIR}/cases/convergence/cavityRe100/NavierStokes/20x20"))
	parser.add_argument("-tolerance", dest="tolerance", help="folder in which the cases for different mesh sizes are present", default="")
	parser.add_argument("-interpolation_type", dest="interpolation_type", help="the type of interpolation used in the direct forcing method", default="")
	parser.add_argument("-device", dest="device", help="CUDA device number", default="")
	parser.add_argument("-run", dest="runSimulations", help="run the cases if this flag is used", action='store_true', default=False)
	parser.add_argument("-remove_mask", dest="removeMask", help="use values from the entire domain", action='store_true', default=False)
	args = parser.parse_args()

	commandOptions = []
	if args.tolerance:
		commandOptions.extend(["-velocityTol", "{}".format(args.tolerance), "-poissonTol", "{}".format(args.tolerance)])
	if args.interpolation_type:
		commandOptions.extend(["-interpolationType", "{}".format(args.interpolation_type)])
	if args.device:
		commandOptions.extend(["-deviceNumber", "{}".format(args.device)])

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
	Q = []
	errL2NormU  = np.zeros(numFolders-1)
	errL2NormV  = np.zeros(numFolders-1)
	errL2NormQ  = np.zeros(numFolders-1)
	errLInfNormU  = np.zeros(numFolders-1)
	errLInfNormV  = np.zeros(numFolders-1)
	meshSize = np.zeros(numFolders-1, dtype=int)

	line()
	print "Case directory: {}".format(args.caseDir)
	print "Folders: ",

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
			try:
				mask_u, mask_v = readMask(folderPath, nx, ny)
				u[:] = u[:]*mask_u[:]
				v[:] = v[:]*mask_v[:]
			except:
				pass

		U.append(np.reshape(u, (ny, nx-1))[stride/2::stride,stride-1::stride])
		V.append(np.reshape(v, (ny-1, nx))[stride-1::stride,stride/2::stride])
		Q.append(np.concatenate([U[fIdx].flatten(),V[fIdx].flatten()]))
		
		print '{},'.format(folder),
		stride = stride*3

	for idx in range(numFolders-1):
		errL2NormU[idx] = la.norm(U[idx+1]-U[idx])
		errL2NormV[idx] = la.norm(V[idx+1]-V[idx])
		errL2NormQ[idx] = la.norm(Q[idx+1]-Q[idx])

		errLInfNormU[idx] = la.norm(U[idx+1]-U[idx], np.inf)
		errLInfNormV[idx] = la.norm(V[idx+1]-V[idx], np.inf)
		
		if idx==0:
			plt.ioff()
			h = initialMeshSpacing
			# u diffs
			x = np.arange(0., 1.+h/2., h)
			y = np.arange(h/2., 1., h)
			X, Y = np.meshgrid(x,y)
			fig = plt.figure()
			ax = fig.add_subplot(111)
			ax.set_xlim([0,1])
			ax.set_ylim([0,1])
			diffV = np.abs(V[idx+1]-V[idx])
			ax.tick_params(labelsize=16)
			CS = ax.pcolormesh(X, Y, diffV, norm=LogNorm(vmin=1e-10, vmax=1))
			fig.gca().set_aspect('equal', adjustable='box')
			cbar = fig.colorbar(CS)
			cbar.ax.tick_params(labelsize=16) 
			if args.removeMask:
				fig.savefig("{}/diffV_nomask.pdf".format(args.caseDir))
			else:
				fig.savefig("{}/diffV.pdf".format(args.caseDir))
			# v diffs
			x = np.arange(h/2., 1., h)
			y = np.arange(0., 1.+h/2., h)
			X, Y = np.meshgrid(x,y)
			fig = plt.figure()
			ax = fig.add_subplot(111)
			ax.set_xlim([0,1])
			ax.set_ylim([0,1])
			diffU = np.abs(U[idx+1]-U[idx])
			ax.tick_params(labelsize=16)
			CS = ax.pcolormesh(X, Y, diffU, norm=LogNorm(vmin=1e-10, vmax=1))
			fig.gca().set_aspect('equal', adjustable='box')
			cbar = fig.colorbar(CS)
			cbar.ax.tick_params(labelsize=16) 
			if args.removeMask:
				fig.savefig("{}/diffU_nomask.pdf".format(args.caseDir))
			else:
				fig.savefig("{}/diffU.pdf".format(args.caseDir))

	L2OrderOfConvergenceU = -np.polyfit(np.log10(meshSize), np.log10(errL2NormU), 1)[0]
	L2OrderOfConvergenceV = -np.polyfit(np.log10(meshSize), np.log10(errL2NormV), 1)[0]
	L2OrderOfConvergenceQ = -np.polyfit(np.log10(meshSize), np.log10(errL2NormQ), 1)[0]

	LInfOrderOfConvergenceU = -np.polyfit(np.log10(meshSize), np.log10(errLInfNormU), 1)[0]
	LInfOrderOfConvergenceV = -np.polyfit(np.log10(meshSize), np.log10(errLInfNormV), 1)[0]
	
	print " "
	line()
	
	print "U:",
	#print "errL2Norms: {}".format(errL2NormU)
	print "Convergence rates:",
	for i in range(len(errL2NormU)-1):
		print "{:.3f}, ".format(np.log(errL2NormU[i]/errL2NormU[i+1])/np.log(3)),
	print "Linear fit convergence rate: {:.3f}".format(L2OrderOfConvergenceU)

	print "V:",
	#print "errL2Norms: {}".format(errL2NormV)
	print "Convergence rates:",
	for i in range(len(errL2NormV)-1):
		print "{:.3f}, ".format(np.log(errL2NormV[i]/errL2NormV[i+1])/np.log(3)),
	print "Linear fit convergence rate: {:.3f}".format(L2OrderOfConvergenceV)

	'''
	print "Q:",
	#print "errL2Norms: {}".format(errL2NormQ)
	print "Convergence rates:",
	for i in range(len(errL2NormQ)-1):
		print "{:.3f}, ".format(np.log(errL2NormQ[i]/errL2NormQ[i+1])/np.log(3)),
	print "Linear fit convergence rate: {:.3f}".format(L2OrderOfConvergenceQ)
	'''
	
	line()

	plt.clf()
	plt.loglog(meshSize, errL2NormU, 's-k', label="$u$".format(L2OrderOfConvergenceU))
	plt.loglog(meshSize, errL2NormV, 'o-k', label="$v$".format(L2OrderOfConvergenceV))
	plt.axis([1, 1e4, 1e-4, 10])
	x  = np.linspace(1, 1e4, 2)
	x1 = 1/x
	x2 = 1/x**2
	plt.loglog(x, x1, '--k', label="First-order convergence")
	plt.loglog(x, x2, ':k', label="Second-order convergence")
	plt.legend(prop={'size':15})
	ax = plt.axes()
	ax.set_xlabel("Mesh size", fontsize=18)
	ax.set_ylabel("$L^2$-norm of differences", fontsize=18)
	ax.tick_params(labelsize=18)
	plt.savefig("{}/L2Convergence.pdf".format(args.caseDir))

	'''
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
	plt.savefig("{}/LInfConvergence.pdf".format(args.caseDir))
	'''

if __name__ == "__main__":
	main()
