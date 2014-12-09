#!/usr/bin/env python

import os
import argparse
import subprocess

def line():
	print '-'*80

def main():
	# Command line options
	parser = argparse.ArgumentParser(description="Calculates the order of convergence.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-folder", dest="case_dir", help="folder in which the cases for different mesh sizes are present", default=os.path.expandvars("${CUIBM_DIR}/cases/convergence/cavityRe100/DirectForcing"))
	parser.add_argument("-interpolation_type", dest="interpolation_type", help="the type of interpolation used in the direct forcing method", default="linear")
	args = parser.parse_args()

	# list of folders where the convergence test cases are present
	folders = sorted(os.walk(args.case_dir).next()[1], key=str)

	for folder in folders:
		runCommand = [os.path.expandvars("${CUIBM_DIR}/scripts/convergence/cavityRe100.py"),
							'-folder', "{}/{}/{}".format(args.case_dir, folder, args.interpolation_type)]
		print " ".join(runCommand)
		subprocess.call(runCommand)

if __name__ == "__main__":
	main()
