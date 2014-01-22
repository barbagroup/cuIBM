#!/usr/bin/env python
import argparse
import numpy as np
import sys

# Generate optionsfile as per command line options
parser = argparse.ArgumentParser()

# Command line options
parser.add_argument("--outfile", dest="outfile", help="name of file generated", default="domain.yaml")
parser.add_argument("--d-blx", type=float, dest="d_blx", help="bottom-left x-coordinate of the domain", default=-15)
parser.add_argument("--d-bly", type=float, dest="d_bly", help="bottom-left y-coordinate of the domain", default=-15)
parser.add_argument("--d-trx", type=float, dest="d_trx", help="top-right x-coordinate of the domain", default=15)
parser.add_argument("--d-try", type=float, dest="d_try", help="top-right y-coordinate of the domain", default=15)
parser.add_argument("--u-blx", type=float, dest="u_blx", help="bottom-left x-coordinate of the uniform region", default=-0.52)
parser.add_argument("--u-bly", type=float, dest="u_bly", help="bottom-left y-coordinate of the uniform region", default=-2)
parser.add_argument("--u-trx", type=float, dest="u_trx", help="top-right x-coordinate of the uniform region", default=3.48)
parser.add_argument("--u-try", type=float, dest="u_try", help="top-right y-coordinate of the uniform region", default=2)
parser.add_argument("--h", type=float, dest="h", help="cell width in uniform region", default=0.004)
parser.add_argument("--sr", type=float, dest="sr", help="stretching ratio", default=1.01)

args = parser.parse_args()

print '-'*120
print "Domain:"
print "\t(%s, %s) to (%s, %s)" % (str(args.d_blx), str(args.d_bly), str(args.d_trx), str(args.d_try))
print 'Uniform region:'
print "\t(%s, %s) to (%s, %s)\n" % (str(args.u_blx), str(args.u_bly), str(args.u_trx), str(args.u_try))

unx = int((args.u_trx-args.u_blx)/args.h + 0.5)
chx = abs((args.u_trx-args.u_blx)/args.h - unx)

uny = int((args.u_try-args.u_bly)/args.h + 0.5)
chy = abs((args.u_try-args.u_bly)/args.h - uny)

if chx > 1e-6 or chy > 1e-6:
	print "Non-integer number of cells in the uniform region! Choose a different h or change the extent of the uniform region.\n"
	sys.exit(0)

f = open(args.outfile, 'w')

# domain description in the x-direction

nx = unx

f.write("- direction: x\n")
f.write("  start: %s\n" % str(args.d_blx))
f.write("  subDomains:\n")

n = 2
L = (args.u_blx - args.d_blx)
h1 = L * (args.sr - 1)/(args.sr**n - 1)
while h1>args.h:
	n = n+1
	h1 = L * (args.sr - 1)/(args.sr**n - 1)
n = n-1
nx = nx + n
h1 = L * (args.sr - 1)/(args.sr**n - 1)
h2 = h1 * args.sr**(n-1)
print "-X :  h1: %.6f    h2: %.6f    AR: %.1f" % (h1, h2, h2/h1)

f.write("    - end: %s\n" % str(args.u_blx))
f.write("      cells: %d\n" % n)
f.write("      stretchRatio: %s\n" % str(1.0/args.sr))

f.write("    - end: %s\n" % str(args.u_trx))
f.write("      cells: %d\n" % unx)
f.write("      stretchRatio: 1.0\n")

n = 2
L = (args.d_trx - args.u_trx)
h1 = L * (args.sr - 1)/(args.sr**n - 1)
while h1>args.h:
	n = n+1
	h1 = L * (args.sr - 1)/(args.sr**n - 1)
n = n-1
nx = nx + n
h1 = L * (args.sr - 1)/(args.sr**n - 1)
h2 = h1 * args.sr**(n-1)
print "+X :  h1: %.6f    h2: %.6f    AR: %.1f" % (h1, h2, h2/h1)

f.write("    - end: %s\n" % str(args.d_trx))
f.write("      cells: %d\n" % n)
f.write("      stretchRatio: %s\n\n" % str(args.sr))

# domain description in the y-direction

ny = uny

f.write("- direction: y\n")
f.write("  start: %s\n" % str(args.d_bly))
f.write("  subDomains:\n")

n = 2
L = (args.u_bly - args.d_bly)
h1 = L * (args.sr - 1)/(args.sr**n - 1)
while h1>args.h:
	n = n+1
	h1 = L * (args.sr - 1)/(args.sr**n - 1)
n = n-1
ny = ny + n
h1 = L * (args.sr - 1)/(args.sr**n - 1)
h2 = h1 * args.sr**(n-1)
print "-Y :  h1: %.6f    h2: %.6f    AR: %.1f" % (h1, h2, h2/h1)

f.write("    - end: %s\n" % str(args.u_bly))
f.write("      cells: %d\n" % n)
f.write("      stretchRatio: %s\n" % str(1.0/args.sr))

f.write("    - end: %s\n" % str(args.u_try))
f.write("      cells: %d\n" % uny)
f.write("      stretchRatio: 1.0\n")

n = 2
L = (args.u_bly - args.d_bly)
h1 = L * (args.sr - 1)/(args.sr**n - 1)
while h1>args.h:
	n = n+1
	h1 = L * (args.sr - 1)/(args.sr**n - 1)
n = n-1
ny = ny + n
h1 = L * (args.sr - 1)/(args.sr**n - 1)
h2 = h1 * args.sr**(n-1)
print "+Y :  h1: %.6f    h2: %.6f    AR: %.1f" % (h1, h2, h2/h1)
print ''

f.write("    - end: %s\n" % str(args.d_try))
f.write("      cells: %d\n" % n)
f.write("      stretchRatio: %s\n" % str(args.sr))

f.close()

print "h  : " + str(args.h)
print "sr : " + str(args.sr) + '\n'
print "Mesh size : " + str(nx) + " x " + str(ny) + '\n'

print "Domain information written to file " + args.outfile

print '-'*120

#time bin/cuIBM -flowFile flows/openFlow.yaml -domainFile domains/gI0203.yaml -bodyFile bodies/snake_0.003.yaml -simulationFile simParams/openFlow.yaml -folderName gI0203 -nu 0.0005 -ibmScheme TairaColonius -dt 0.0003 -nt 266665 -nsave 53333 -deviceNumber 3 -uPerturb 0.0
#time bin/cuIBM -flowFile flows/openFlow.yaml -domainFile domains/gI0103.yaml -bodyFile bodies/snake_0.003.yaml -simulationFile simParams/openFlow.yaml -folderName gI0103 -nu 0.0005 -ibmScheme TairaColonius -dt 0.0003 -nt 266665 -nsave 53333 -deviceNumber 3 -uPerturb 0.0

