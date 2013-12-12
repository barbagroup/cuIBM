#!/usr/bin/env python

import argparse
import os
import os.path
import sys
import numpy as np
import readData as rd

cuibmFolder = os.path.expandvars("${CUIBM_DIR}")

# Parse command line options
parser = argparse.ArgumentParser()
parser.add_argument("--Re", dest="Re", help="Reynolds number", default='100')
args = parser.parse_args()
Re = args.Re

if Re=='100':
	domainFile     = '/cavity0100.yaml'
	simulationFile = '/cavity0100.yaml'
	nu             = '0.01'
	uCol           = '2'
	vCol           = '6'
elif Re=='1000':
	domainFile     = '/cavity1000.yaml'
	simulationFile = '/cavity1000.yaml'
	nu             = '0.001'	
	uCol           = '3'
	vCol           = '7'
else:
	print "Unavailable option for Reynolds number. Choose 100 or 1000."
	sys.exit()

validationData = '/cavity-GGS82.txt'
folderName     = cuibmFolder + '/validation/lidDrivenCavity/Re' + Re
flowFile       = cuibmFolder + '/flows/cavity.yaml'
bodyFile       = cuibmFolder + '/bodies/empty.yaml'
domainFile     = cuibmFolder + '/domains' + domainFile
simulationFile = cuibmFolder + '/simParams' + simulationFile
validationData = cuibmFolder + '/validation-data' + validationData
execPath       = cuibmFolder + '/bin/cuIBM'

print "\n"+"-"*120
print "Running the validation case for flow in a lid-driven cavity at Reynolds number %s" % Re
print "-"*120+"\n"

runCommand = "%s -flowFile %s -domainFile %s -bodyFile %s -simulationFile %s -folderName %s -nu %s" % (execPath, flowFile, domainFile, bodyFile, simulationFile, folderName, nu)
print runCommand+"\n"
os.system(runCommand)

nt, startStep, nsave = rd.readSimulationParameters(folderName)
nx, ny, dx, dy, xu, yu, xv, yv = rd.readGridData(folderName)
u, v = rd.readVelocityData(folderName, nt, nx, ny, dx, dy)

f = open(folderName+"/u.txt", 'w')
for j in range(len(yu)):
	f.write("%f\t%f\n" % (yu[j], u[j*(nx-1)+nx/2-1]) )
f.close()

f = open(folderName+"/v.txt", 'w')
for i in range(len(xv)):
	f.write("%f\t%f\n" % (xv[i], v[(ny/2-1)*nx + i]) )
f.close()

print "-"*120
print "Plotting the centerline velocities for flow in a lid-driven cavity at Reynolds number %s" % Re
print "-"*120

gnuplotFile    = cuibmFolder + '/scripts/validation/cavityRe' + Re + '.plt'
uOutFile       = folderName + '/uRe' + Re + '.pdf'
vOutFile       = folderName + '/vRe' + Re + '.pdf'

f = open(gnuplotFile, 'w')
f.write("reset;\nset terminal pdf enhanced color font 'Palatino, 11' size 15cm, 15cm;\n\n");
f.write("set title 'Velocity along the vertical centerline %s'\n" % Re)
f.write("set xlabel 'y-coordinate'\n")
f.write("set ylabel 'Centerline u-velocity'\n")
f.write("set output '%s'\n" % uOutFile)
f.write("plot [0:1] [-0.7:1.3] \\\n")
f.write("'%s' u 1:%s w p pt 13 ps 1 lc rgb '#4B5ED7' title 'Ghia et al, 1982', \\\n" % (validationData, uCol) )
f.write("'%s/u.txt' u 1:2 w l lw 5 lc rgb '#228A4C' title 'cuIBM (Taira and Colonius, 2005)'\n" % folderName)

f.write("\n")
f.write("reset;\nset terminal pdf enhanced color font 'Palatino, 11' size 15cm, 15cm;\n\n");
f.write("set title 'Velocity along the horizontal centerline %s'\n" % Re)
f.write("set xlabel 'x-coordinate'\n")
f.write("set ylabel 'Centerline v-velocity'\n")
f.write("set output '%s'\n" % vOutFile)
f.write("plot [0:1] [-0.7:1.3] \\\n")
f.write("'%s' u 5:%s w p pt 13 ps 1 lc rgb '#4B5ED7' title 'Ghia et al, 1982', \\\n" % (validationData, vCol) )
f.write("'%s/v.txt' u 1:2 w l lw 5 lc rgb '#228A4C' title 'cuIBM (Taira and Colonius, 2005)'\n" % folderName)

f.close()

print "\nCreated gnuplot script "+gnuplotFile
runCommand = "gnuplot "+gnuplotFile
os.system(runCommand)
print "\nDone plotting! Plots saved in folder "+folderName+"\n"
