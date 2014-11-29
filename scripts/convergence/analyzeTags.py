#!/usr/bin/env python

import argparse
import os
import numpy
import numpy.linalg as la
import sys
sys.path.insert(0, os.path.expandvars("${CUIBM_DIR}/scripts/python"))
from readData import readEta

epsilon = 1e-8

def alpha_predicted(x):
	if x>=8/9. and x<9/9.:
		return 1+numpy.log((-18.*x-6)/(-18.*x+10))/numpy.log(3)
	elif x>=7/9. and x<8/9.:
		return 1+numpy.log((-18.*x-6)/(-18.*x+11))/numpy.log(3)
	elif x>=6/9. and x<7/9.:
		return 1+numpy.log((-18.*x-6)/(-18.*x+12))/numpy.log(3)
	elif x>=5/9. and x<6/9.:
		return 1+numpy.log((-18.*x-3)/(-18.*x+4))/numpy.log(3)
	elif x>=4/9. and x<5/9.:
		return 1+numpy.log((-18.*x-3)/(-18.*x+5))/numpy.log(3)
	elif x>=3/9. and x<4/9.:
		return 1+numpy.log((-18.*x-3)/(-18.*x+6))/numpy.log(3)
	elif x>=2/9. and x<3/9.:
		return 1+numpy.log((-18.*x)/(-18.*x-2))/numpy.log(3)
	elif x>=1/9. and x<2/9.:
		return 1+numpy.log((-18.*x)/(-18.*x-1))/numpy.log(3)
	elif x>=0/9. and x<1/9.:
		return 1.

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
	return point_list

def print_common_points(point_lists, name):
	print "\n{} grids:".format(name)
	for i, point_list in enumerate(point_lists):
		print "{}: {} boundary points".format(i, len(point_list))
	
	for i in range(len(point_lists)-1):
		total = len(point_lists[i])
		common = 0
		for p0 in point_lists[i]:
			for p1 in point_lists[i+1]:
				if equals(p0, p1):
					common+=1
		print "Boundary points unchanged between grids {} and {}: {}/{}, {:.1f}%".format(i, i+1, common, total, common*100.0/total)

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

	alpha_vectorized = numpy.vectorize(alpha_predicted)
	print "\nu: ",
	for folder in folders[:2]:
		eta_u, _ = readEta(args.caseDir + "/" + folder)
		print "{:.3f}".format(la.norm(alpha_vectorized(eta_u),1)/len(eta_u)),
		#print "{:.3f}".format(la.norm(alpha_vectorized(eta_u),2)/len(eta_u)),
		#print "{:.3f}".format(alpha_predicted(la.norm(eta_u,1)/numpy.sqrt(len(eta_u)))),
		#print "{:.3f}".format(alpha_predicted(la.norm(eta_u,2)/numpy.sqrt(len(eta_u)))),
	print "\nv: ",
	for folder in folders[:2]:
		_, eta_v = readEta(args.caseDir + "/" + folder)
		print "{:.3f}".format(la.norm(alpha_vectorized(eta_v),1)/len(eta_v)),
		#print "{:.3f}".format(la.norm(alpha_vectorized(eta_v),2)/len(eta_v)),
		#print "{:.3f}".format(alpha_predicted(la.norm(eta_v,1)/numpy.sqrt(len(eta_v)))),
		#print "{:.3f}".format(alpha_predicted(la.norm(eta_v,2)/numpy.sqrt(len(eta_v)))),
	print " "

if __name__=="__main__":
	main()