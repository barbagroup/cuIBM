#!/usr/bin/env python

import os
import os.path

print "\n"+"-"*120
print "Running the validation case for flow impulsively started uniform  over a stationary cylinder at Reynolds number 40"
print "-"*120+"\n"
runCommand = os.path.expandvars("${CUIBM_DIR}/bin/cuIBM \
-flowFile ${CUIBM_DIR}/flows/openFlow.yaml \
-domainFile ${CUIBM_DIR}/domains/openFlow0.025.yaml \
-bodyFile ${CUIBM_DIR}/bodies/cylinder0.025.yaml \
-simulationFile ${CUIBM_DIR}/simParams/cylinder0040.yaml \
-folderName ${CUIBM_DIR}/cylinderRe40 \
-nu 0.025");
print runCommand+"\n"
os.system(runCommand)

print "-"*120
print "Plotting the drag coefficient for impulsively started uniform flow over a stationary cylinder at Reynolds number 40"
print "-"*120
runCommand = os.path.expandvars("gnuplot ${CUIBM_DIR}/scripts/validation/cylinderRe40.plt")
os.system(runCommand)
filePath = os.path.expandvars("${CUIBM_DIR}/cylinderRe40/cylRe40Drag.pdf")
print "\nDONE! Plot saved to file "+filePath+"\n"
