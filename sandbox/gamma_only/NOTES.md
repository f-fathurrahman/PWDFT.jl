# Gamma-point Trick Implementation Notes

## Overview

Gamma-point trick is used for real quantities that depend of $\mathbf{G}$-vectors. Using the trick, only (about) half of the $\mathbf{G}$-vectors are stored. The expansion coefficient for other half can be calculated from:
$$
c(\mathbf{G}) = c^{*}(-\mathbf{G})
$$

The basic data structure is now `PWGridGamma` which is mainly similar to `PWGrid`.
The main difference is with the `GVectorsGamma` and `GVectorsWGamma`. Other fields of `PWGridGamma` are similar to the usual `PWGrid`.

## Generating GVectorsGamma

The algorithm for generating half set of $\mathbf{G}$-vectors are adapted from the `ggen` subroutine of PWSCF. I simplify a bit and use adapt it according to the convention that I have used before.

Mapping from 3d index to linear index.

## Operations involving wavefunction

Notes: I have not made any comparison with PWSCF implementation. I come up with my own implementation by studying some operations with special arrays.

Test array:

dot product

overlap matrix

Wavefunction orthonormalization: ortho_GS_gamma and ortho_sqrt_gamma.


FIXME: Need to compare gradient calculation: calc_grad using Gamma-trick vs usual.


## Calculating electron density

Blah:

## Hamiltonian

New struct instead of parametric  struct: `HamiltonianGamma`. Several things that change

Construction of local potential