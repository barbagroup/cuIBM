#!/usr/bin/env python

import argparse
import os
import os.path
import sys

cuibmFolder = os.path.expandvars("${CUIBM_DIR}")

# Parse command line options
parser = argparse.ArgumentParser()
parser.add_argument("--h", dest="h", help="Cell width near the body", default='0.004')
args = parser.parse_args()

if args.h == '0.004':
	domainFile     = '/0.004.yaml'
	bodyFile       = '/04_0.004.yaml'
	dt             = '0.0004'
	nt             = '30000'
elif args.h == '0.002':
	domainFile     = '/0.002.yaml'
	bodyFile       = '/04_0.002.yaml'
	dt             = '0.0002'
	nt             = '60000'
else:
	print "Unavailable option for h. Choose 0.004 or 0.002."
	sys.exit()

folderName     = cuibmFolder + '/validation/naca0008/04_' + args.h
flowFile       = cuibmFolder + '/flows/openFlow.yaml'
simulationFile = cuibmFolder + '/simParams/naca0008.yaml'
domainFile     = cuibmFolder + '/domains/naca0008' + domainFile
bodyFile       = cuibmFolder + '/bodies/naca0008' + bodyFile
execPath       = cuibmFolder + '/bin/cuIBM'

print "\n"+"-"*120
print "Running the validation case for Re=2000 flow over a NACA 0008 airfoil at angle of attack 4 degrees (h = %s)" % args.h
print "-"*120+"\n"

runCommand = "%s -flowFile %s -domainFile %s -bodyFile %s -simulationFile %s -folderName %s -nu 0.0005 -dt %s -nt %s" % (execPath, flowFile, domainFile, bodyFile, simulationFile, folderName, dt, nt)
print runCommand+"\n"
os.system(runCommand)

print "-"*120
print "Plotting the drag coefficient for Re=2000 flow over NACA 0008 airfoil at angle of attack 4 degrees (h = %s) " % args.h
print "-"*120

gnuplotFile    = cuibmFolder + '/scripts/validation/naca0008_' + args.h + '.plt'
outFile        = folderName + '/NACA0008_' + args.h + '.pdf'

f = open(gnuplotFile, 'w')
f.write("reset;\nset terminal pdf enhanced color font 'Palatino, 11' size 20cm, 15cm;\n\n");
f.write("set title 'Flow over a NACA 0008 airfoil at angle of attach 4^o'\n")
f.write("set xlabel 'Non-dimensional time'\n")
f.write("set ylabel 'Force Coefficient'\n")
f.write("set output '%s'\n" % outFile)
f.write("plot [0:12] [0:0.5] \\\n")
f.write("'%s/forces' u 1:(2*$3) w l lw 5 lc rgb '#4B5ED7' title 'Lift coefficient', \\\n" % folderName)
f.write("'%s/forces' u 1:(2*$2) w l lw 5 lc rgb '#228A4C' title 'Drag coefficient'\n" % folderName)
f.close()
print "\nCreated gnuplot script "+gnuplotFile

runCommand = "gnuplot "+gnuplotFile
os.system(runCommand)
print "\nDone plotting! Plot saved to file "+outFile+"\n"
