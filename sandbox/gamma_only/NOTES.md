# Gamma-point Trick Implementation Notes

Gamma-point trick is used for real quantities that depend of $\mathbf{G}$-vectors. Using the trick, only (about) half of the $\mathbf{G}$-vectors are stored. The expansion coefficient for other half can be calculated from:
$$
c(\mathbf{G}) = c^{*}(-\mathbf{G})
$$


The basic data structure is now `PWGridGamma` which is mainly similar to `PWGrid`.
The main difference is with the `GVectorsGamma` and `GVectorsWGamma`.
The algorithm for generating half set of G-vectors are adapted from ggen subroutine of PWSCF.
