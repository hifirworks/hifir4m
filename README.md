# HIFIR High-Level Interface for MATLAB and GNU Octave #

`hifir4m` allows accessing HIFIR from MATLAB and GNU Octave using pre-compiled high-level or intermediate-level interfaces. Its underlying library, HIFIR, is a modern preconditoining technology that combines efficient *multi-level ILU* (*MLILU*) and robust *rank-revealing QR* (*RRQR*) factorization at its coarsest level. HIFIR and `hifir4m` enables users to precondition ill-conditioned and singular systems for `gmres`. In additoin, `hifir4m` also offers a high-level solver named `gmresHifir`, with a similar interface as MATLAB's built-in `gmres` function but uses right-preconditioned GMRES and FGMRES with `HIFIR` as a right-preconditioner.

## Installation ##

Clone this project to your preferred location. Then start MATLAB or GNU Octave under the directory that contains `hifir4m`, or run the command

```matlab
>> run('/path/to/hifir4m/startup_hifir')
```

(and replace `/path/to/hifir4m/` to the directory that contains `hifir4m`) bo build `hifir4m` and load its path. It will build the *mex* kernels if needed by linking to MATLAB and Octave's built-in BLAS/LAPACK libraries for the low-level QRCP.

Note that for the first time, `hifir4m` will download C++ package [HIFIR](https://github.com/hifirworks/hifir) while building *mex* kernels. If you don't have access to network during building *mex* kernels, then you can obtain HIFIR beforehand and put its source in `hifir4m/hifir-[hifir-version]` folder; for instance, you can run the following command to download the C++ HIFIR package

```console
cd /path/to/hifir4m
wget -qO- https://github.com/hifirworks/hifir/archive/refs/tags/v`cat VERSION`.tar.gz|tar xzf -
```

## Usage ##

The easiest way to use `hifir4m` is to call the `gmresHifir` interface. For example,

```matlab
>> [x, flag, relres, iter, reshis, times] = gmresHifir(A, b);
```

where `A` is a MATLAB's built-in sparse matrix or a `MATLAB` `struct` containing the filds of `row_ptr`, `col_ind`, `vals` of a standard CRS storage format, and `b` is a right-hand-side vector (or RHS vectors with two columns).

To access the intermediate-level interfaces of `hifir4m`, please see `gmresHifir.m` for the calling sequence of `hifCreate`, `hifApply`, and `hifDestroy`.

## Copyright and Licenses ##

`HIFIR`, `hifir4m`, and `hifir4py` are developed by the NumGeom Research Group at Stony Brook University.

The software suite is released under a dual-license model. For academic users, individual users, or open-source software developers, you can use HIFIR under the GPLv3+ license free of charge for research and evaluation purpose. For commercial users, separate commercial licenses are available through the Stony Brook University.  For inqueries regarding commercial licenses, please contact Prof. Xiangmin Jiao <xiangmin.jiao@stonybrook.edu>.

## How to Cite `HIFIR` ##

If you use `HIFIR`,  `hifir4m`, or `hifir4py` in your research for nonsingular systems, please cite the `HILUCSI` paper:

```bibtex
@Article{chen2021hilucsi,
  author  = {Chen, Qiao and Ghai, Aditi and Jiao, Xiangmin},
  title   = {{HILUCSI}: Simple, robust, and fast multilevel {ILU} for 
             large-scale saddle-point problems from {PDE}s},
  journal = {Numer. Linear Algebra Appl.},
  year    = {2021},
  doi     = {10.1002/nla.2400},
```

If you use them to solve highly ill-conditioned of singular systems, please cite the `HIFIR` papers:

```bibtex
@Article{jiao2020approximate,
  author  = {Xiangmin Jiao and Qiao Chen},
  journal = {SIAM J. Matrix Anal. Appl.},
  title   = {Approximate Generalized Inverses with Iterative Refinement for 
             $\epsilon$-Accurate Preconditioning of Singular Systems},
  year    = {2021},
  note    = {To appear},
}

@Article{chen2021hifir,
  author  = {Chen, Qiao and Jiao, Xiangmin},
  title   = {{HIFIR}: Hybrid Incomplete Factorization with Iterative Refinement 
             for Preconditioning Ill-conditioned and Singular Systems},
  journal = {arxiv},
  year    = {2021},
  note    = {arXiv:2106.09877},
}
```

## Contacts ##

- Qiao Chen, <qiao.chen@stonybrook.edu>, <benechiao@gmail.com>
- Xiangmin Jiao, <xiangmin.jiao@stonybrook.edu>, <xmjiao@gmail.com>
