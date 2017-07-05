# Flow around a flying-snake cross-section (Re=1000, AoA=30deg)

This directory contains the input files to run a simulation of the flow around a flying-snake cross-section at Reynolds number 1000 and a 30-degree angle of attack (with the freestream flow).

This example is based on the simulation reported in Krishnan et al. (2014).
It computes the numerical solution over 200,000 time steps with a time-step size of 0.0004 (for a total of 80 non-dimensional time units).

**WARNING**: this is a long run; for example, it takes about 65 hours to complete the 200,000 time steps on a K20 GPU device.
If you wish to test this example on a smaller number of time steps, you can change the value of the parameter `nt` in the input file `simParams.yaml`.

The post-processing script `scripts/plotForceCoefficients.py` plots the instantaneous force coefficients and computes the mean force coefficients and Strouhal number between 32 and 64 time units (range used in Krishnan et al., 2014); you can adapt the time range to your convenience.

The post-processing script `scripts/plotVorticity.py` plots the vorticity field at saved time steps (controlled by the parameter `nsave` in the file `simParams.yaml`).


_References:_
* Krishnan, A., Socha, J. J., Vlachos, P. V., Barba, L. A. (2014). "Lift and wakes of flying snakes," _Physics of Fluids_, __26__:031901, [10.1063/1.4866444](https://dx.doi.org/10.1063/1.4866444). Preprint on arXiv: https://arxiv.org/abs/1309.2969