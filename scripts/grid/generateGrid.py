#!/usr/bin/env python
import argparse
import numpy as np
import sys
import os
import re

# parse command line options
parser = argparse.ArgumentParser(description="Generates a cartesian mesh with a uniform region surrounded by a stretched grid", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--output", dest="output", help="name of file generated", default="domain.yaml")
parser.add_argument("--input", dest="input", help="name of input file", default="gridOptions")
args = parser.parse_args()

inFile = os.path.expandvars("${CUIBM_DIR}/scripts/grid/"+args.input);

f = open(inFile, 'r')
for line in f:
	b = filter(None, re.split('\[|\]|\n|:|,| ', line))
	if b != []:
		if b[0] == 'DomainBottomLeft':
			d_blx = float(b[1])
			d_bly = float(b[2])
		elif b[0] == 'DomainTopRight':
			d_trx = float(b[1])
			d_try = float(b[2])
		elif b[0] == 'UniformRegionBottomLeft':
			u_blx = float(b[1])
			u_bly = float(b[2])
		elif b[0] == 'UniformRegionTopRight':
			u_trx = float(b[1])
			u_try = float(b[2])
		elif b[0] == 'FinestMeshSpacing':
			h = float(b[1])
		elif b[0] == 'StretchingRatio':
			sr = float(b[1])

print '-'*120
print "Domain:"
print "\t(%s, %s) to (%s, %s)\n" % (str(d_blx), str(d_bly), str(d_trx), str(d_try))
print 'Uniform region:'
print "\t(%s, %s) to (%s, %s)\n" % (str(u_blx), str(u_bly), str(u_trx), str(u_try))

unx = int((u_trx-u_blx)/h + 0.5)
chx = abs((u_trx-u_blx)/h - unx)

uny = int((u_try-u_bly)/h + 0.5)
chy = abs((u_try-u_bly)/h - uny)

if chx > 1e-6 or chy > 1e-6:
	print "Non-integer number of cells in the uniform region! Choose a different h or change the extent of the uniform region.\n"
	sys.exit(0)

f = open(args.output, 'w')

# domain description in the x-direction

nx = unx

f.write("- direction: x\n")
f.write("  start: %s\n" % str(d_blx))
f.write("  subDomains:\n")

n = 2
L = (u_blx - d_blx)
h1 = L * (sr - 1)/(sr**n - 1)
while h1>h:
	n = n+1
	h1 = L * (sr - 1)/(sr**n - 1)
n = n-1
nx = nx + n
h1 = L * (sr - 1)/(sr**n - 1)
h2 = h1 * sr**(n-1)
print "-X :  h1: %.6f    h2: %.6f    AR: %.1f" % (h1, h2, h2/h1)

f.write("    - end: %s\n" % str(u_blx))
f.write("      cells: %d\n" % n)
f.write("      stretchRatio: %s\n" % str(1.0/sr))

f.write("    - end: %s\n" % str(u_trx))
f.write("      cells: %d\n" % unx)
f.write("      stretchRatio: 1.0\n")

n = 2
L = (d_trx - u_trx)
h1 = L * (sr - 1)/(sr**n - 1)
while h1>h:
	n = n+1
	h1 = L * (sr - 1)/(sr**n - 1)
n = n-1
nx = nx + n
h1 = L * (sr - 1)/(sr**n - 1)
h2 = h1 * sr**(n-1)
print "+X :  h1: %.6f    h2: %.6f    AR: %.1f" % (h1, h2, h2/h1)

f.write("    - end: %s\n" % str(d_trx))
f.write("      cells: %d\n" % n)
f.write("      stretchRatio: %s\n\n" % str(sr))

# domain description in the y-direction

ny = uny

f.write("- direction: y\n")
f.write("  start: %s\n" % str(d_bly))
f.write("  subDomains:\n")

n = 2
L = (u_bly - d_bly)
h1 = L * (sr - 1)/(sr**n - 1)
while h1>h:
	n = n+1
	h1 = L * (sr - 1)/(sr**n - 1)
n = n-1
ny = ny + n
h1 = L * (sr - 1)/(sr**n - 1)
h2 = h1 * sr**(n-1)
print "-Y :  h1: %.6f    h2: %.6f    AR: %.1f" % (h1, h2, h2/h1)

f.write("    - end: %s\n" % str(u_bly))
f.write("      cells: %d\n" % n)
f.write("      stretchRatio: %s\n" % str(1.0/sr))

f.write("    - end: %s\n" % str(u_try))
f.write("      cells: %d\n" % uny)
f.write("      stretchRatio: 1.0\n")

n = 2
L = (u_bly - d_bly)
h1 = L * (sr - 1)/(sr**n - 1)
while h1>h:
	n = n+1
	h1 = L * (sr - 1)/(sr**n - 1)
n = n-1
ny = ny + n
h1 = L * (sr - 1)/(sr**n - 1)
h2 = h1 * sr**(n-1)
print "+Y :  h1: %.6f    h2: %.6f    AR: %.1f" % (h1, h2, h2/h1)
print ''

f.write("    - end: %s\n" % str(d_try))
f.write("      cells: %d\n" % n)
f.write("      stretchRatio: %s\n" % str(sr))

f.close()

print "h  : " + str(h)
print "sr : " + str(sr) + '\n'
print "Mesh size : " + str(nx) + " x " + str(ny) + '\n'

print "Domain information written to file " + args.output

print '-'*120
