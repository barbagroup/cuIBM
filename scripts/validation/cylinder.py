#!/usr/bin/env python

import argparse
import os
import os.path
import sys

cuibmFolder = os.path.expandvars("${CUIBM_DIR}")

# Parse command line options
parser = argparse.ArgumentParser(description="Runs the validation case for impulsively started flow over a circular cylinder for a specified Reynolds number", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--Re", dest="Re", help="Reynolds number", default='40')
args = parser.parse_args()
Re = args.Re

if Re=='40':
	validationData = '/cylinderRe40-KL95.txt'
	yMax           = '6' 
elif Re=='550':
	validationData = '/cylinderRe550-KL95.txt'
	yMax           = '2' 
elif Re=='3000':
	validationData = '/cylinderRe3000-KL95.txt'
	yMax           = '2' 
else:
	print "Unavailable option for Reynolds number. Choose 40, 550 or 3000."
	sys.exit()

execPath       = cuibmFolder + '/bin/cuIBM'
caseFolder     = cuibmFolder + '/validation/cylinder/Re' + Re
validationData = cuibmFolder + '/validation-data' + validationData

print "\n"+"-"*120
print "Running the validation case for flow over an impulsively started circular cylinder at Reynolds number %s" % Re
print "-"*120+"\n"

runCommand = "%s -caseFolder %s" % (execPath, caseFolder)
print runCommand+"\n"
os.system(runCommand)

print "-"*120
print "Plotting the drag coefficient for flow over an impulsively started circular cylinder at Reynolds number %s" % Re
print "-"*120

gnuplotFile    = caseFolder + '/cylinderRe' + Re + '.plt'
outFile        = caseFolder + '/cylRe' + Re + 'Drag.pdf'

f = open(gnuplotFile, 'w')
f.write("reset;\nset terminal pdf enhanced color font 'Palatino, 11' size 20cm, 15cm;\n\n");
f.write("set title 'Flow over an impulsively started cylinder at Reynolds number %s'\n" % Re)
f.write("set xlabel 'Non-dimensional time'\n")
f.write("set ylabel 'Drag Coefficient'\n")
f.write("set output '%s'\n" % outFile)
f.write("plot [0:3] [0:%s] \\\n" % yMax)
f.write("'%s' u (0.5*$1):2 w p pt 13 ps 1 lc rgb '#4B5ED7' title 'Koumoutsakos and Leonard, 1995', \\\n" % validationData)
f.write("'%s/forces' u 1:(2*$2) w l lw 5 lc rgb '#228A4C' title 'cuIBM (Taira and Colonius, 2005)'\n" % caseFolder)
f.close()

print "\nCreated gnuplot script "+gnuplotFile
runCommand = "gnuplot "+gnuplotFile
os.system(runCommand)
print "\nDone plotting! Plot saved to file "+outFile+"\n"
