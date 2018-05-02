`ffr-LFDFT` is an experimental package to solve [electronic structure problems](https://en.wikipedia.org/wiki/Electronic_structure)
based on [density functional theory](https://en.wikipedia.org/wiki/Density_functional_theory)
(DFT)
and [Kohn-Sham equations](https://en.wikipedia.org/wiki/Kohn%E2%80%93Sham_equations).

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
ecutwfc_Ry = 30.0
pspfiles = ["../pseudopotentials/pade_gth/C-q4.gth",
            "../pseudopotentials/pade_gth/H-q1.gth"]
Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5 )
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
