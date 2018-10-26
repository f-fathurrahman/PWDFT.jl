# PWDFT.jl

`PWDFT.jl` is a package to solve
[electronic structure problems](https://en.wikipedia.org/wiki/Electronic_structure)
based on [density functional theory](https://en.wikipedia.org/wiki/Density_functional_theory)
(DFT) and [Kohn-Sham equations](https://en.wikipedia.org/wiki/Kohn%E2%80%93Sham_equations).
It is written in [Julia programming language](https://julialang.org).

The Kohn-Sham orbitals are expanded using plane wave basis. This basis set is
very popular within solid-state community and is also used in several electronic
structure package such as Quantum ESPRESSO, ABINIT, VASP, etc.

## Requirements

- [Julia](https://julialang.org/downloads) version > 0.7:
  with `FFTW` and `SpecialFunctions` packages installed.
- [LibXC](https://gitlab.com/libxc/libxc) version > 3.0:
  which needs to be compiled and and installed separately.
- [spglib](https://github.com/atztogo/spglib): which needs to be compiled and installed
  separately.
- A working C compiler to compile LibXC and spglib.

## Installation

- Compile and install LibXC and spglib.

- Install Julia package `FFTW` and `SpecialFunctions`. The following
  command can be run under Julia console.

```Julia
Pkg.add("FFTW")
Pkg.add("SpecialFunctions")
```

- Currently, this package is not yet registered. You can use this package by
  cloning this repository under the `$HOME/.julia/dev` directory.

- Create symlink under `$HOME/.julia/dev` to point to `PWDFT.jl`

```bash
ln -fs PWDFT.jl PWDFT
```

- Open the file `extlibs/extlibs.jl` using text editor. Edit the
  `LIBXC` and `LIBSYMSPG` according to your LibXC and spglib installations.
   These variables should point to the appropriate dynamic libraries
   of LibXC and spglib, respectively. Examples:
```
@checked_lib LIBXC "/usr/local/libxc-3.0.0/lib/libxc.so"
@checked_lib LIBSYMSPG "/usr/local/spglib-1.10.4/lib/libsymspg.so"
```

- To make sure that the package is installed correctly, you can load the package
  and verify that there are no error messages during precompilation step.
  You can do this by typing the following in the Julia console.

```julia
using PWDFT
```

- Change directory to `examples` and run the following in the terminal.

```
julia run.jl "atom_H/main_H_atom.jl"
```
  
	The above command will calculate total energy of hyrogen atom by
  SCF method.


## Units

`PWDFT.jl` internally uses Hartree atomic units,
(energy in Hartree and length in bohr).

## A simple work flow

- create an instance of `Atoms`:

```julia
atoms = Atoms(xyz_file="CH4.xyz", LatVecs=gen_lattice_sc(16.0))
```

- create an instance of `Hamiltonian`:

```julia
ecutwfc = 15.0 # in Hartree
pspfiles = ["../pseudopotentials/pade_gth/C-q4.gth",
            "../pseudopotentials/pade_gth/H-q1.gth"]
Ham = Hamiltonian( atoms, pspfiles, ecutwfc )
```

- solve the Kohn-Sham problem

```julia
KS_solve_SCF!( Ham, betamix=0.2 )  # using SCF (self-consistent field) method
# or
KS_solve_Emin_PCG!( Ham ) # direct minimization using preconditioned conjugate gradient
```

## Band structure calculations

![Band structure of silicon (fcc)](images/bands_Si_fcc.svg)

Please see [bandstructure_Si_fcc.jl](sandbox/bandstructure_Si_fcc.jl) as
an example of how this can be obtained.

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
