# HIFIR MATLAB Interface #

`hifir4m` is a MATLAB interface for HIFIR preconditioner, plus a right-preconditioned GMRES/FGMRES solvers. HIFIR is a hybrid of *multi-level ILU* (*MLILU*) and *rank-revealing QR* (*RRQR*) factorization preconditioner for ill-conditioned and singular systems, which is developed by NumGeom research group at Stony Brook University.

## Installation ##

Copy this project to your preferred location, then under MATLAB command window, simply type

```matlab
>> startup
```

It will try to fetch the C++ kernel (if needed) from BitBucket repository, then build the *mex* kernels (again, if needed). It is worth noting that HIFIR4M links to MATLAB built-in BLAS/LAPACK, thanks to its flexible template and header-only design.

## Usage ##

Look at the random matrix example under `examples` directory. Also, use `help` command for docstrings. This project allows you (1) solving linear systems (reusing same preconditioner) with initial guess and (2) accessing only the preconditioner part for using in other project (e.g., multigrid smoother). For most application codes, (1) is what you want, as you can see, the design is for nonlinear and/or transient problems.

## Copying ##

This project inherits the GLPv3+ license from HIFIR, for more, check `LICENSE` file.

## Contact(s) ##

Qiao Chen, <qiao.chen@stonybrook.edu>, <benechiao@gmail.com>
