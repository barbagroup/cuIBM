#!/usr/bin/env python
import argparse
import numpy as np

# parse command line options
parser = argparse.ArgumentParser(description="Creates a body with uniform segment lengths from a given body input", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--ds", type=float, dest="ds", help="segment length", default=0.004)
parser.add_argument("--infile", dest="infile", help="name of input file", default="snake_cross_section.txt")
parser.add_argument("--outfile", dest="outfile", help="prefix of output file generated", default="snake")
args = parser.parse_args()

x0 = np.arange(0)
y0 = np.arange(0)
nb = 0
f = open(args.infile, 'r')
a = f.readline().strip().split()
x0 = np.append(x0, float(a[0]))
y0 = np.append(y0, float(a[1]))
xmin = x0[nb]
xmax = x0[nb]
xc = 0.
yc = 0.
xc = x0[nb]
yc = y0[nb]
nb = nb+1
while True:
	a = f.readline().strip().split()
	if a==[]:
		break
	x0 = np.append(x0, float(a[0]))
	y0 = np.append(y0, float(a[1]))
	xc = xc + x0[nb]
	yc = yc + y0[nb]
	if x0[nb] < xmin:
		xmin = x0[nb]
	if x0[nb] > xmax:
		xmax = x0[nb]
	nb = nb+1
xc = xc/nb
yc = yc/nb
c = xmax - xmin
x0 = (x0 - xc)/c
y0 = (y0 - yc)/c
f.close()

ds = args.ds
cur = 0
nxt = 1
x = np.zeros(1)
y = np.zeros(1)
x[cur] = x0[0]
y[cur] = y0[0]
while nxt < nb:
	dist = np.sqrt( (x0[nxt]-x[cur])**2 + (y0[nxt]-y[cur])**2 )
	while dist > ds:
		x = np.append( x, x[cur] + (x0[nxt]-x[cur])*ds/dist )
		y = np.append( y, y[cur] + (y0[nxt]-y[cur])*ds/dist )
		cur = cur+1
		dist = np.sqrt( (x0[nxt]-x[cur])**2 + (y0[nxt]-y[cur])**2 )
	while dist < ds:
		nxt = nxt+1
		if nxt == nb:
			break
		dist = np.sqrt( (x0[nxt]-x[cur])**2 + (y0[nxt]-y[cur])**2 )
	if nxt == nb:
		break
	# solve for the next point
	#'''
	upp = 1.0
	low = 0.0
	while abs(dist-ds)/ds > 1e-5:
		eta = 0.5*(upp+low)
		dist = np.sqrt( ( (1-eta)*x0[nxt-1] + eta*x0[nxt] - x[cur] )**2 + ( (1-eta)*y0[nxt-1] + eta*y0[nxt] - y[cur] )**2 )
		if dist < ds:
			low = eta
		else:
			upp = eta
	'''
	A = (x0[nxt]-x0[nxt-1])**2 + (y0[nxt]-y0[nxt-1])**2
	B = - ( (x[cur]-x0[nxt-1])*(x0[nxt]-x0[nxt-1]) + (y[cur]-y0[nxt-1])*(y0[nxt]-y0[nxt-1]) )
	C = (x0[nxt-1]**2 + x[cur]**2 + y0[nxt-1]**2 + y[cur]**2 - ds**2)
	if B*B-4*A*C < 0:
		print cur, A, B, C, B*B-4*A*C
		break
	eta = (-B + np.sqrt(B*B-4*A*C))/(2*A)
	'''
	x = np.append( x, (1-eta)*x0[nxt-1] + eta*x0[nxt] )
	y = np.append( y, (1-eta)*y0[nxt-1] + eta*y0[nxt] )
	cur = cur + 1
dist = np.sqrt( (x0[0]-x[cur])**2 + (y0[0]-y[cur])**2 )
while dist > ds:
	x = np.append( x, x[cur] + (x0[0]-x[cur])*ds/dist )
	y = np.append( y, y[cur] + (y0[0]-y[cur])*ds/dist )
	cur = cur+1
	dist = np.sqrt( (x0[0]-x[cur])**2 + (y0[0]-y[cur])**2 )
if dist < 0.5*ds:
	x[cur] = 0.5*(x[cur-1] + x0[0])
	y[cur] = 0.5*(y[cur-1] + y0[0])

print ""
print "input element width          : %f" % ds
dist = np.sqrt( (x[0]-x[cur])**2 + (y[0]-y[cur])**2 )
print "width of last element        : %f" % dist
dist = np.sqrt( (x[cur]-x[cur-1])**2 + (y[cur]-y[cur-1])**2 )
print "width of penultimate element : %f" % dist
dist = np.sqrt( (x[cur-1]-x[cur-2])**2 + (y[cur-1]-y[cur-2])**2 )
print "width of third last element  : %f" % dist
print "\ntotal number of boundary points: %d\n" % len(x)

outFile = args.outfile + "_" + str(ds) + ".bdy"
g = open(outFile, 'w')
g.write("%d\n" % len(x))
for i in range(len(x)):
	g.write( "%f\t%f\n" % (x[i], y[i]) )
g.close()
