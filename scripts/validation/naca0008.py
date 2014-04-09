#!/usr/bin/env python

import argparse
import os
import os.path
import sys

cuibmFolder = os.path.expandvars("${CUIBM_DIR}")

# Parse command line options
parser = argparse.ArgumentParser(description="Runs the validation case for flow over a NACA 0008 airfoil at angle of attack 4 degrees at Reynolds number 2000", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--h", dest="h", help="Cell width near the body", default='0.004')
args = parser.parse_args()

if args.h!='0.004' and args.h!='0.002':
	print "Unavailable option for h. Choose 0.004 or 0.002."
	sys.exit()

caseFolder     = cuibmFolder + '/validation/NACA0008/AoA04_' + args.h
execPath       = cuibmFolder + '/bin/cuIBM'

print "\n"+"-"*120
print "Running the validation case for Re=2000 flow over a NACA 0008 airfoil at angle of attack 4 degrees (h = %s)" % args.h
print "-"*120+"\n"

runCommand = "%s -caseFolder %s" % (execPath, caseFolder)
print runCommand+"\n"
os.system(runCommand)

print "-"*120
print "Plotting the drag coefficient for Re=2000 flow over NACA 0008 airfoil at angle of attack 4 degrees (h = %s) " % args.h
print "-"*120

gnuplotFile    = caseFolder + '/AoA04_' + args.h + '.plt'
outFile        = caseFolder + '/AoA04_' + args.h + '.pdf'

f = open(gnuplotFile, 'w')
f.write("reset;\nset terminal pdf enhanced color font 'Palatino, 11' size 20cm, 15cm;\n\n");
f.write("set title 'Flow over a NACA 0008 airfoil at angle of attach 4^o'\n")
f.write("set xlabel 'Non-dimensional time'\n")
f.write("set ylabel 'Force Coefficient'\n")
f.write("set output '%s'\n" % outFile)
f.write("plot [0:12] [0:0.5] \\\n")
f.write("'%s/forces' u 1:(2*$3) w l lw 5 lc rgb '#4B5ED7' title 'Lift coefficient', \\\n" % caseFolder)
f.write("'%s/forces' u 1:(2*$2) w l lw 5 lc rgb '#228A4C' title 'Drag coefficient'\n" % caseFolder)
f.close()
print "\nCreated gnuplot script "+gnuplotFile

runCommand = "gnuplot "+gnuplotFile
os.system(runCommand)
print "\nDone plotting! Plot saved to file "+outFile+"\n"
