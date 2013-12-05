#!/usr/bin/env python

import os
import os.path

print "\n"+"-"*120
print "Running the validation case for flow impulsively started uniform  over a stationary cylinder at Reynolds number 3000"
print "-"*120+"\n"
runCommand = os.path.expandvars("${CUIBM_DIR}/bin/cuIBM \
-flowFile ${CUIBM_DIR}/flows/openFlow.yaml \
-domainFile ${CUIBM_DIR}/domains/openFlow0.004.yaml \
-bodyFile ${CUIBM_DIR}/bodies/cylinder0.004.yaml \
-simulationFile ${CUIBM_DIR}/simParams/cylinder3000.yaml \
-folderName ${CUIBM_DIR}/cylinderRe3000 \
-nu 0.00033333333333");
print runCommand+"\n"
os.system(runCommand)

print "-"*120
print "Plotting the drag coefficient for impulsively started uniform flow over a stationary cylinder at Reynolds number 3000"
print "-"*120
runCommand = os.path.expandvars("gnuplot ${CUIBM_DIR}/scripts/validation/cylinderRe3000.plt")
os.system(runCommand)
filePath = os.path.expandvars("${CUIBM_DIR}/cylinderRe3000/cylRe3000Drag.pdf")
print "\nDONE! Plot saved to file "+filePath+"\n"
