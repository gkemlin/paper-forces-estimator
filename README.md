Julia code for asymptotic estimators for quantites of interest in electronic
structure. It can be used to reproduce the main results from our paper
_Practical error bounds for properties in plane-wave electronic structure
calculations_ by Éric Cancès, Geneviève Dusson, Gaspard Kemlin and Antoine
Levitt, available on [arXiv](https://arxiv.org/abs/2111.01470).

This code is used with the refined reference solutions we introduced in the
paper. A small example of the estimates for the forces is available
[here](https://docs.dftk.org/stable/examples/error_estimates_forces/) in
the DFTK documentation.

# Dependencies:
Julia 1.7 with the following Julia libraries:
- [DFTK.jl](https://dftk.org) v0.5.0 for the simulation environment (this code
  might not work with more recent versions);
- ForwardDiff for the computation of `dF` with automatic differentiation;
- LinearAlgebra, IterativeSolvers, LinearMaps, FFTW for the linear algebra
  related computations;
- HDF5, CSV, Tables for saving and reading results.

# Usage
To perform the computations, first open the Julia shell with `julia --project`
from the location of this repository and then run
```
using Pkg
Pkg.instantiate()
```
to install the required dependencies.
Then, you can run all the computations at once with
`include("run_all_computations.jl")`. This will go in every folder `case/` where
`case` is successively `silicon`, `GaAs` and `TiO2` with the setting described
in our paper. Each folder contains:
- a file `case_forces.jl` which run the computations with the actual Schur
  complement from our paper to improve the estimate on the forces and then save
  the results in the file `forces_case.h5`;
- a file `generate_data_case_forces.jl` which generates various `.csv` files for
  plots;
- different `.tex` files that read the associated `.csv` files to generates
  `.pdf` files:
    - `error_case` files represent the estimation on the error `Π(P-P_*)` in the
      H^1-norm, as in Figure 2 (middle), without the bad error bound;
    - `forces_case` files represent the estimation of the error `F-F*` as in
    Figure 3 (right);
    - `diff_forces_case` files represent the improvement of this estimations, as
      in Figures 7 and 9.

Note that, after running the computations, only `.h5` and `.csv` files have been
generated. You need to compile the `.pdf` files with the associated `.tex` files
and your favorite LaTeX compiler.

# Computational time
As we have chosen a quite precise reference grid, with many coarser grids,
we ran all the computations on 32 Intel Xeon E5-2667 (3.30GHz) cores, which takes
around half a day.

# Contact
This is research code, not necessarily user-friendly, actively maintened or extremely robust.
If you have questions or encounter problems, get in touch!


