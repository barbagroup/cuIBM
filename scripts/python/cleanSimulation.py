#!/usr/bin/env/ python

# file: cleanSimulation.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# description: Clean a cuIBM simulation.


import os
import argparse


def read_inputs():
  """Parses the command-line."""
  # create parser
  parser = argparse.ArgumentParser(description='Clean PetIBM case',
                        formatter_class= argparse.ArgumentDefaultsHelpFormatter)
  # fill parser with arguments
  parser.add_argument('--case', dest='case_directory', type=str, 
                      default=os.getcwd(),
                      help='directory of the PetIBM simulation')
  parser.add_argument('--no-images', dest='images', action='store_false',
                      help='does not remove the images folder')
  parser.add_argument('--no-grid', dest='grid', action='store_false',
                      help='does not remove the grid file')
  parser.add_argument('--no-solutions', dest='solutions', action='store_false',
                      help='does not remove the numrical solution folders')
  parser.add_argument('--no-forces', dest='forces', action='store_false',
                      help='does not remove the forces data file')
  parser.add_argument('--no-logs', dest='logs', action='store_false',
                      help='does not remove log files '
                           '(iterations, run.info)')
  parser.set_defaults(images=True, grid=True, solutions=True, 
                      forces=True, logs=True)
  return parser.parse_args()


def main():
  """Cleans a cuIBM simulation."""
  # parser command-line
  args = read_inputs()

  # get different parts to clean
  parts = {}
  if args.images:
    parts['images'] = '%s/images' % args.case_directory
  if args.grid:
    parts['grid'] = '%s/grid' % args.case_directory
  if args.solutions:
    parts['solutions'] = '%s/0*' % args.case_directory
  if args.forces:
    parts['forces'] = '%s/forces' % args.case_directory
  if args.logs:
    parts['logs'] = ('%s/iterations %s/run.info' 
                     % (args.case_directory, args.case_directory))

  # remove appropriate files/folders
  print '[case-directory] %s' % args.case_directory
  for key, part in parts.iteritems():
    print '\t--> removing %s ...' % key
    os.system('rm -rf %s' % part)


if __name__ == '__main__':
  main()