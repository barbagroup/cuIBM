---
title: 'cuIBM: a GPU-based immersed boundary method code'
tags:
  - GPU
  - Computational Fluid Dynamics
  - Immersed-Boundary Method
  - CUSP library
authors:
 - name: Anush Krishnan
   affiliation: nuTonomy Inc. (previously at Boston University)
 - name: Olivier Mesnard
   orcid:  0000-0001-5335-7853
   affiliation: 1
 - name: Lorena A. Barba
   orcid: 0000-0001-5812-2711
   affiliation: 1
affiliations:
 - name: The George Washington University
   index: 1
date: 14 July 2016
bibliography: paper.bib
---

# Summary

cuIBM solves the two-dimensional Navier-Stokes equations with an immersed-boundary method on structured Cartesian grids. 
The equations are spatially discretized with a finite-difference technique and temporally integrated via a projection approach seen as an approximate block-LU decomposition (@Perot1993).
cuIBM implements various immersed-boundary techniques that fit into the framework of Perot's projection method.
Among them, the immersed-boundary projection approach from @TairaColonius2007, the direct-forcing method from @FadlunEtAl2000, and a second-order accurate direct-forcing method (@Krishnan2015).

cuIBM exploits NVIDIA GPU hardware by solving the linear systems and applying stencil operations on a single GPU device.
For that purpose, we use [CUSP](https://github.com/cusplibrary/cusplibrary), an open-source C++ library for sparse linear algebra on CUDA-capable GPUs.
cuIBM is written in C++ with CUDA kernels.

As an example, we used cuIBM to investigate the aerodynamics of an anatomically-accurate cross-section of the gliding snake, *Chrysopelea Paradisi* (@KrishnanEtAl2014).

# References
