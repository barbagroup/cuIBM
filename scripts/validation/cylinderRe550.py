#!/usr/bin/env python

import os
import os.path

print "\n"+"-"*120
print "Running the validation case for flow impulsively started uniform  over a stationary cylinder at Reynolds number 550"
print "-"*120+"\n"
runCommand = os.path.expandvars("${CUIBM_DIR}/bin/cuIBM \
-flowFile ${CUIBM_DIR}/flows/openFlow.yaml \
-domainFile ${CUIBM_DIR}/domains/openFlow0.010.yaml \
-bodyFile ${CUIBM_DIR}/bodies/cylinder0.010.yaml \
-simulationFile ${CUIBM_DIR}/simParams/cylinder0550.yaml \
-folderName ${CUIBM_DIR}/cylinderRe550 \
-nu 0.00181818181818");
print runCommand+"\n"
os.system(runCommand)

print "-"*120
print "Plotting the drag coefficient for impulsively started uniform flow over a stationary cylinder at Reynolds number 550"
print "-"*120
runCommand = os.path.expandvars("gnuplot ${CUIBM_DIR}/scripts/validation/cylinderRe550.plt")
os.system(runCommand)
filePath = os.path.expandvars("${CUIBM_DIR}/cylinderRe550/cylRe550Drag.pdf")
print "\nDONE! Plot saved to file "+filePath+"\n"
