# Implementation Notes

## About Julia programming language

Even though Julia is a dynamic language, I tried to always use type annotations
whenever possible.

## Atoms

`Atoms` can be used to represent molecular or crystalline structures.

It is implemented in file `src/Atoms.jl`.

It has the following fields:

- `Natoms::Int64`

- `Nspecies::Int64`

- `positions::Array{Float64,2}`

- `atm2species::Array{Int64,1}`

- `atsymbs::Array{String,1}`

- `SpeciesSymbols::Array{String,1}`

- `LatVecs::Array{Float64,2}`

- `Zvals::Array{Float64,1}`


## Hamiltonian

An instance `Hamiltonian` is a central object in the calculation.
It is used to stores various instances of other important types
such as atoms, plane wave grids, pseudopotentials, etc.

It is implemented in file `src/Hamiltonian.jl`.

To create an instance of `Hamiltonian`, we need to provide at least
two arguments:

- `atoms::Atoms`: an instance of `Atoms`

- `pspfiles::Array{String,1}`: a list of strings specifying the
  locations of pseudopotentials used in the
  calculations. Note that, the order should be the same as species ordering
  of `atoms`, i.e. `pspfiles[isp]` is the path of
  pseudopotentials of species with symbols `atoms.SpeciesSymbols[isp]`.

## Plane wave basis

G-vectors
G-vectorsW

## Wavefunctions

Using `Array{ComplexF64,2}`.

General wavefunction on kpoints

## Potentials

## Energies

## Solving Kohn-Sham problems

## Self-consistent field

## Eigensolver

## Direct energy minimization via nonlinear conjugate gradient method

## Direct constrained minimization method of Yang

## Chebyshev filtered subspace iteration SCF
