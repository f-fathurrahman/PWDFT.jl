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

Plane wave basis is described by type `PWGrid`. This type is defined in file `PWGrid.jl`.
It has the following fields:

- `ecutwfc::Float64`: cutoff for wave function expansion

- `ecutrho::Float64`: cutoff for electron density expansion, for norm-converving
  pseudopotential: `ecutrho = 4*ecutwfc`.

- `Ns::Tuple{Int64,Int64,Int64}`: parameters defining real-space grid points.

- `LatVecs::Array{Float64,2}`: lattice vectors of unit cell

- `RecVecs::Array{Float64,2}`: reciprocal lattice vectors

- `CellVolume::Float64`: the volume of real-space unit cell

- `r::Array{Float64,2}`: real-space grid points

- `gvec::GVectors`: an instace of `Gvectors`: for potentials and density expansion

- `gvecw::GVectorsW`: an instace of `GvectorsW`, for wave function expansion

- `planfw::FFTW.cFFTWPlan{Complex{Float64},-1,false,3}`: FFTW forward plan

- `planbw::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,false,3- Float64}`: FFTW backward plan

The following constructor can be used to create an instance of `PWGrid`:
```julia
function PWGrid( ecutwfc::Float64,
                 LatVecs::Array{Float64,2};
                 kpoints=nothing )
```

### G-vectors

G-vectors are described by type `GVectors`. It is defined in file `PWGrid.jl`.



### G-vectorsW

G-vectors for wave function expansion is described by type `GVectorsW`.


## Wave functions

Using `Array{ComplexF64,2}`.

General wavefunction on kpoints

## Potentials

## Pseudopotentials

## Energies

## Solving Kohn-Sham problems

## Self-consistent field

## Eigensolver

## Direct energy minimization via nonlinear conjugate gradient method

## Direct constrained minimization method of Yang

## Chebyshev filtered subspace iteration SCF
