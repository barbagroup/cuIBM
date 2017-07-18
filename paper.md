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
date: 9 June 2017
bibliography: paper.bib
---

# Summary

cuIBM solves the two-dimensional Navier-Stokes equations with an immersed-boundary method on structured Cartesian grids.
With this solution approach, we release the constraint for the mesh-grid to fit to the surface of a body immersed in the fluid and use simple and easy-to-generate fixed Cartesian grid.
(Thus, we do not have to re-generate the grid when the body is moving.)
The governing equations are modified to take into account the presence of the immersed boundary and solved over the extended domain.
The equations are spatially discretized with a finite-difference technique and temporally integrated via a projection approach seen as an approximate block-LU decomposition (@Perot1993).
cuIBM implements various immersed-boundary techniques that fit into the framework of Perot's projection method.
Among them are the immersed-boundary projection approach from @TairaColonius2007, the direct-forcing method from @FadlunEtAl2000, and a second-order accurate direct-forcing method (@Krishnan2015).

cuIBM is written in C++ and exploits NVIDIA GPU hardware using CUDA and [CUSP](https://github.com/cusplibrary/cusplibrary), 
an open-source C++ library for sparse linear algebra on CUDA-capable GPUs.
cuIBM solves the linear systems of equations and applies stencil operations on a single GPU device.

cuIBM generated the results published in @KrishnanEtAl2014, a study of gliding-snake aerodynamics using an anatomically accurate cross-section of the snake *Chrysopelea Paradisi*.

# References
