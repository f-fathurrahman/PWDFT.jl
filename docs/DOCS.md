## Describing an atomic system: `Atoms`

All atomic systems are assumed to be periodic.

The definition of type `Atoms` is given below.

```julia
mutable struct Atoms
    Natoms::Int64
    Nspecies::Int64
    positions::Array{Float64,2}
    atm2species::Array{Int64,1}
    atsymbs::Array{String,1}
    SpeciesSymbols::Array{String,1}
    LatVecs::Array{Float64,2}
    Zvals::Array{Float64,1}
end
```

Information about `LatVecs` and `Zvals` are also available
from `PWGrid` and `PsPots`. They are included to reduce number of
required arguments to several functions.

Currently, the following functions are provided to initialize an `Atoms`:

```julia
atoms = Atoms() # dummy constructor

atoms = init_atoms_xyz(filexyz; in_bohr=false, verbose=false)

atoms = init_atoms_xyz_string(filexyz; in_bohr=false, verbose=false)
```

Note that, `LatVecs` must be set manually by:
```julia
atoms.LatVecs = 16*eye(3) # for example
```

`Zvals` is set when constructing `PWHamiltonian`.
