#!/usr/bin/env python
import os
import sys
import argparse
import numpy as np
import readData as rd
import matplotlib.pyplot as plt

# Generate optionsfile as per command line options
parser = argparse.ArgumentParser()
parser.add_argument("--folder", dest="folder", help="Case folder", default=".")
args = parser.parse_args()

folder = args.folder
nt, startStep, nsave, dt = rd.readSimulationParameters(folder)
nx, ny, dx, dy, _, yu, xv, _ = rd.readGridData(folder)

savePoints = 0
k = startStep+nsave
while k < (nt+1):
	savePoints+=1
	k+=nsave

T = 10.
h = dx[0]
totalVorticity = np.zeros(savePoints)
circulation = np.zeros(savePoints)
time = np.zeros(savePoints)

k = startStep+nsave
idx = 0
while k < (nt+1):
	totalVorticity[idx] = 0.
	time[idx] = k*dt
	
	vortFile = args.folder + ("/%07d" % k) + "/vorticity"
	
	f = open(vortFile, 'r')
	i = 0
	
	for line in f:
		try:
			totalVorticity[idx] += float(line.split()[2])
		except IndexError:
			pass
	
	totalVorticity[idx] *= h*h
	circulation[idx] = -np.sin(np.pi*time[idx]/T)
	
	k +=nsave
	idx += 1

plt.ioff()
tV, = plt.plot(time, totalVorticity, 'D')
tC, = plt.plot(time, circulation)
plt.title("Comparison of Total Vorticity and Theoretical Circulation at Re=10,000")
plt.legend([tV, tC], ['Total Vorticity', 'Theoretical circulation'])
plt.xlabel('Time')
plt.ylabel('Total circulation')
plt.axis([0, 80, -1.5, 1.5])
plt.savefig("out.png")
