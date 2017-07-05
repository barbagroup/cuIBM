# Flapping wing at Reynolds number 75

This directory contains the input files to run a simulation of the flow around a flapping wing at Reynolds number 75.

**WARNING**: this is a long run; for example, it takes about 10 hours to complete the 16,000 time steps on a K20 GPU device.
If you wish to test this example on a smaller number of time steps, you can change the value of the parameter `nt` in the input file `simParams.yaml`.

The post-processing script `scripts/plotVorticity.py` plots the vorticity field at saved time steps (controlled by the parameter `nsave` in the file `simParams.yaml`).
