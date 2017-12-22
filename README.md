# TTK4550 Engineering Cybernetics Specialization Project

## General center

This repository contains a reference implementation of the algorithm presented in the project «Optimization-Based Control for
Multi-Agent Deployment - Remedying convergence issues by computation of a general center». The code is not optimized (eg. uses SVD in places where QR-decomposition could be more efficient).

## Usage
Add the folder GeneralCenter to your matlab path. The general center and the corresponding Chebyshev radius may be retrieved with the following function call
```
% Polytope in H-rep, Ax <= b.
% x_g is the general center
% r_opt is the Chebyshev radius
[x_g,r_opt] = GeneralCenter(A,b)
```
