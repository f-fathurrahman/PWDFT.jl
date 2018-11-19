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
- [SPGLIB](https://github.com/atztogo/spglib): which needs to be compiled and installed
  separately.
- A working C compiler to compile LibXC and SPGLIB.

## Installation

### Compile and install LibXC and SPGLIB.

Configure and install LibXC

```bash
cd libxc-3.0.0 # please change according to your LibXC version
./configure --prefix=/usr/local/libxc-3.0.0 --disable-fortran --enable-shared
make
make install # may need root privilege
```

Configure and install SPGLIB (using Cmake)

```bash
cd spglib-master
mkdir build
cd build
cmake -D CMAKE_INSTALL_PREFIX=/usr/local/spglib-1.10.4
make
make install  # may need root privilege
```


### Install Julia's packages: `FFTW` and `SpecialFunctions`

This can be done by executing the following commands at Julia console.

```julia
using Pkg
Pkg.add("FFTW")
Pkg.add("SpecialFunctions")
```

### Setup `PWDFT.jl` as Julia package

Currently, this package is not yet registered. So, `Pkg.add("PWDFT")` will not work (yet).

We have two alternatives:

1. Using Julia's package manager to install directly from the repository URL:

```julia
Pkg.add(PackageSpec(url="https://github.com/f-fathurrahman/PWDFT.jl"))
```

2. Using Julia development directory. We will use `$HOME/.julia/dev` for this.
   To enable `$HOME/.julia/dev` directory, we need to modify the Julia's
  `LOAD_PATH` variable. Add the following line in your
  `$HOME/.julia/config/startup.jl`.

```julia
push!(LOAD_PATH, expanduser("~/.julia/dev"))
```

  After this has been set, you can download the the package as zip file (using Github) or
  clone this repository to your computer.

  If you download the zip file, extract the zip file under
  `$HOME/.julia/dev$`. You need to rename the extracted directory
  to `PWDFT` (with no `.jl` extension).

  Alternatively, create symlink under `$HOME/.julia/dev`
  to point to you cloned (or extracted) `PWDFT.jl` directory. The link name should not
  contain the `.jl` part. For example:

```bash
ln -fs /path/to/PWDFT.jl $HOME/.julia/dev/PWDFT
```

### Edit the `extlibs/extlibs.jl` file

Open the file `extlibs/extlibs.jl` (found under PWDFT.jl directory)
using text editor. Edit the `LIBXC` and `LIBSYMSPG`
according to your LibXC and SPGLIB installations.
These variables should point to the appropriate dynamic libraries
of LibXC and spglib, respectively. For example:

```julia
@checked_lib LIBXC "/usr/local/libxc-3.0.0/lib/libxc.so"
@checked_lib LIBSYMSPG "/usr/local/spglib-1.10.4/lib/libsymspg.so"
```

To make sure that the package is installed correctly, you can load the package
and verify that there are no error messages during precompilation step.
You can do this by typing the following in the Julia console.

```julia
using PWDFT
```

Change directory to `examples` and run the following in the terminal.

```
julia run.jl "atom_H/main_H_atom.jl"
```
  
The above command will calculate total energy of hydrogen atom
by SCF method.


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
