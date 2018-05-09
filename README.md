`PWDFT.jl` is a package to solve
[electronic structure problems](https://en.wikipedia.org/wiki/Electronic_structure)
based on [density functional theory](https://en.wikipedia.org/wiki/Density_functional_theory)
(DFT) and [Kohn-Sham equations](https://en.wikipedia.org/wiki/Kohn%E2%80%93Sham_equations).
It is written in [Julia programming language](https://julialang.org).

The Kohn-Sham orbitals are expanded using plane wave basis. This basis set is
very popular within solid-state community and is also used in several electronic
structure package such as Quantum ESPRESSO, ABINIT, VASP, etc.

## Calculation steps

- create `Atoms` object

```julia
atoms = init_atoms_xyz("CH4.xyz")
atoms.LatVecs = 16.0*eye(3)
```

- create `PWHamiltonian` object

```julia
ecutwfc = 15.0 # in Hartree
pspfiles = ["../pseudopotentials/pade_gth/C-q4.gth",
            "../pseudopotentials/pade_gth/H-q1.gth"]
Ham = PWHamiltonian( atoms, pspfiles, ecutwfc )
```

- calculate interaction energy between ions (atomic centers):

```julia
Ham.energies.NN = calc_E_NN( atoms )
```

- solve Kohn-Sham equations using any of the following methods

```julia
KS_solve_SCF!( Ham, β=0.2 )  # using SCF (self-consistent field) method
# or
KS_solve_Emin_PCG!( Ham ) # direct minimization using preconditioned conjugate gradient
```

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


