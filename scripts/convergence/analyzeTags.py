#!/usr/bin/env python

import argparse
import os
import numpy

epsilon = 1e-8

def equals(a, b):
	if numpy.abs(a[0]-b[0]) < epsilon and numpy.abs(a[1]-b[1]) < epsilon:
		return True
	return False

def read_points(filename):
	point_list = numpy.zeros((0,2))
	prev = ""
	for l in open(filename):
		if prev == "" or prev=="\n" or prev=="\r":
			a = l.split()
			point = numpy.array([[float(a[0]), float(a[1])]])
			point_list = numpy.append(point_list, point, axis=0)
		prev = l
	#numpy.savetxt(folder+".txt", point_list)
	return point_list

def print_common_points(point_lists, name):
	print "\n{} grids:".format(name)
	for i in range(len(point_lists)-1):
		total = len(point_lists[i])
		common = 0
		for p0 in point_lists[i]:
			for p1 in point_lists[i+1]:
				if equals(p0, p1):
					common+=1
		print "Common: {}/{}, {:.1f}%".format(common, total, common*100.0/total)

def main():
	# Command line options
	parser = argparse.ArgumentParser(description="Calculates the common points on successive grids.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-folder", dest="caseDir", help="folder in which the cases for different mesh sizes are present", default=".")
	args = parser.parse_args()
	
	# list of folders from which velocity data is to be obtained
	folders = sorted(os.walk(args.caseDir).next()[1], key=int)
	num_grids = len(folders)
	
	# list of arrays that contain the tagged points for each grid
	point_lists_u = []
	point_lists_v = []

	for folder in folders:
		p_list = read_points(args.caseDir + "/" + folder + "/tagx.txt")
		point_lists_u.append(p_list)
		p_list = read_points(args.caseDir + "/" + folder + "/tagy.txt")
		point_lists_v.append(p_list)

	print_common_points(point_lists_u, "u")
	print_common_points(point_lists_v, "v")

if __name__=="__main__":
	main()