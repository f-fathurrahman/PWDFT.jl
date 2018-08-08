`PWDFT.jl` is a package to solve
[electronic structure problems](https://en.wikipedia.org/wiki/Electronic_structure)
based on [density functional theory](https://en.wikipedia.org/wiki/Density_functional_theory)
(DFT) and [Kohn-Sham equations](https://en.wikipedia.org/wiki/Kohn%E2%80%93Sham_equations).
It is written in [Julia programming language](https://julialang.org).

The Kohn-Sham orbitals are expanded using plane wave basis. This basis set is
very popular within solid-state community and is also used in several electronic
structure package such as Quantum ESPRESSO, ABINIT, VASP, etc.

## Requirements

- Julia 0.7: with `FFTW` and `SpecialFunctions` packages installed.
- LibXC 3.0: which needs to be compiled and and installed separately.
- spglib: included in this repository and needs to be compiled manually.

## Installation

Currently, this package is not yet registered. You can install this package by
cloning this repository under `.julia/dev` directory.
To make sure the package is installed correctly, you can load the package
and verify that there are no error messages after the loading.

```julia
using PWDFT
```

## Units

`PWDFT.jl` internally uses Hartree atomic units
(energy in Hartree and length in bohr).

## Quick start

- create `Atoms` object

```julia
atoms = init_atoms_xyz("CH4.xyz")
atoms.LatVecs = 16.0*eye(3)
```

- create `Hamiltonian` object

```julia
ecutwfc = 15.0 # in Hartree
pspfiles = ["../pseudopotentials/pade_gth/C-q4.gth",
            "../pseudopotentials/pade_gth/H-q1.gth"]
Ham = Hamiltonian( atoms, pspfiles, ecutwfc )
```

- solve Kohn-Sham equations using any of the following methods

```julia
KS_solve_SCF!( Ham, betamix=0.2 )  # using SCF (self-consistent field) method
# or
KS_solve_Emin_PCG!( Ham ) # direct minimization using preconditioned conjugate gradient
```

## Band structure calculations

![Band structure of silicon (fcc)](images/bands_Si_fcc.svg)

Please see [bandstructure_Si_fcc.jl](sandbox/bandstructure_Si_fcc.jl) as
an example of how this can be obtained.


## Implementation notes

Please see [this](docs/ImplementationNotes.md).

## Some references

Articles:

- M. Bockstedte, A. Kley, J. Neugebauer and M. Scheffler. Density-functional theory
  calculations for polyatomic systems:Electronic structure, static and elastic properties
  and ab initio molecular dynamics. *Comp. Phys. Commun.* **107**, 187 (1997).

- Sohrab Ismail-Beigi and T.A. Arias. New algebraic formulation of density functional calculation.
  *Comp. Phys. Comm.* **128**, 1-45 (2000)


Books:

- Richard Milton Martin. *Electronic Structure: Basic Theory and Practical Methods*.
  Cambridge University Press, 2004.

- Jorge Kohanoff. *Electronic Structure Calculations for Solids and Molecules:
  Theory and Computational Methods*.
  Cambridge University Press, 2006.

- Dominik Marx and Jürg Hutter. *Ab Initio Molecular Dynamics: Basic Theory and
  Advanced Methods*. Cambridge University Press, 2009.
