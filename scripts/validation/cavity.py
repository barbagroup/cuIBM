#!/usr/bin/env python

import argparse
import numpy as np

import os
import os.path
cuibmFolder = os.path.expandvars("${CUIBM_DIR}")

import sys
sys.path.insert(0, cuibmFolder+'/scripts/python')

import readData as rd

# Parse command line options
parser = argparse.ArgumentParser(description="Runs the validation case for flow in a lid-driven cavity for a specified Reynolds number", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--Re", dest="Re", help="Reynolds number", default='100')
args = parser.parse_args()
Re = args.Re

if Re=='100':
	uCol           = '2'
	vCol           = '7'
elif Re=='1000':
	uCol           = '3'
	vCol           = '8'
elif Re=='10000':
	uCol           = '5'
	vCol           = '10'
else:
	print "Unavailable option for Reynolds number. Choose 100, 1000 or 10000."
	sys.exit()

validationData = '/cavity-GGS82.txt'
caseFolder     = cuibmFolder + '/validation/lidDrivenCavity/Re' + Re
validationData = cuibmFolder + '/validation-data' + validationData
execPath       = cuibmFolder + '/bin/cuIBM'

print "\n"+"-"*120
print "Running the validation case for flow in a lid-driven cavity at Reynolds number %s" % Re
print "-"*120+"\n"

runCommand = "%s -caseFolder %s" % (execPath, caseFolder)
print runCommand+"\n"
os.system(runCommand)

nt, _, _, _ = rd.readSimulationParameters(caseFolder)
nx, ny, dx, dy, _, yu, xv, _ = rd.readGridData(caseFolder)
u, v = rd.readVelocityData(caseFolder, nt, nx, ny, dx, dy)

f = open(caseFolder+"/u.txt", 'w')
for j in range(len(yu)):
	f.write("%f\t%f\n" % (yu[j], u[j*(nx-1)+nx/2-1]) )
f.close()

f = open(caseFolder+"/v.txt", 'w')
for i in range(len(xv)):
	f.write("%f\t%f\n" % (xv[i], v[(ny/2-1)*nx + i]) )
f.close()

print "-"*120
print "Plotting the centerline velocities for flow in a lid-driven cavity at Reynolds number %s" % Re
print "-"*120

gnuplotFile    = caseFolder + '/cavityRe' + Re + '.plt'
uOutFile       = caseFolder + '/uRe' + Re + '.pdf'
vOutFile       = caseFolder + '/vRe' + Re + '.pdf'

f = open(gnuplotFile, 'w')
f.write("reset;\nset terminal pdf enhanced color font 'Palatino, 11' size 15cm, 15cm;\n\n");
f.write("set title 'Velocity along the vertical centerline (Re=%s)'\n" % Re)
f.write("set xlabel 'y-coordinate'\n")
f.write("set ylabel 'Centerline u-velocity'\n")
f.write("set output '%s'\n" % uOutFile)
f.write("plot [0:1] [-0.7:1.3] \\\n")
f.write("'%s' u 1:%s w p pt 13 ps 1 lc rgb '#4B5ED7' title 'Ghia et al, 1982', \\\n" % (validationData, uCol) )
f.write("'%s/u.txt' u 1:2 w l lw 5 lc rgb '#228A4C' title 'cuIBM (Taira and Colonius, 2005)'\n" % caseFolder)

f.write("\n")
f.write("reset;\nset terminal pdf enhanced color font 'Palatino, 11' size 15cm, 15cm;\n\n");
f.write("set title 'Velocity along the horizontal centerline (Re=%s)'\n" % Re)
f.write("set xlabel 'x-coordinate'\n")
f.write("set ylabel 'Centerline v-velocity'\n")
f.write("set output '%s'\n" % vOutFile)
f.write("plot [0:1] [-0.7:1.3] \\\n")
f.write("'%s' u 6:%s w p pt 13 ps 1 lc rgb '#4B5ED7' title 'Ghia et al, 1982', \\\n" % (validationData, vCol) )
f.write("'%s/v.txt' u 1:2 w l lw 5 lc rgb '#228A4C' title 'cuIBM (Taira and Colonius, 2005)'\n" % caseFolder)

f.close()

print "\nCreated gnuplot script "+gnuplotFile
runCommand = "gnuplot "+gnuplotFile
os.system(runCommand)
print "\nDone plotting! Plots saved in folder " + caseFolder + "\n"
