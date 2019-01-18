# Initializing `Atoms`

From xyz file: supply the path to xyz file as string and
set the lattice vectors:

```julia
atoms = Atoms(xyz_file="file.xyz", LatVecs=gen_lattice_sc(16.0))
```

In extended xyz file, the lattice vectors information is included
(along with several others information, if any):

```julia
atoms = Atoms(ext_xyz_file="file.xyz")
```

For crystalline systems, using keyword argument `xyz_string_frac`
is sometimes convenient:

```julia
atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true,
        LatVecs=gen_lattice_fcc(10.2631))
```

**IMPORTANT**
We need to be careful to also specify `in_bohr` keyword to get
the correct coordinates in bohr (which is used internally in `PWDFT.jl`).



# Referring or including files in `sandbox` (or other dirs in `PWDFT.jl`)

```julia
using PWDFT
const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT,"pseudopotentials","pade_gth")

pspfiles = [joinpath(DIR_PSP,"Ag-q11.gth")]
```


# Generating lattice vectors

Lattice vectors are simply 3x3 array. We can set it manually or use
one of functions defined in
[gen_lattice_pwscf.jl](../src/gen_lattice_pwscf.jl)
for generating lattice vectors for Bravais lattices that used
in Quantum ESPRESSO's PWSCF.

# Using Babel to generate xyz file from SMILES

```
babel file.smi file.sdf
babel file.sdf file.xyz
```

Use `babel -h` to autogenerate hydrogens.



# Setting up pseudopotentials

One can use the function `get_default_psp(::Atoms)` to get default
pseudopotentials set for a given instance of `Atoms`.

Currently, it is not part of main `PWDFT.jl` package. It is located
under `sandbox` subdirectory of `PWDFT.jl` distribution.

```julia
using PWDFT

DIR_PWDFT = jointpath(dirname(pathof(PWDFT)),"..")
include(jointpath(DIRPWDFT,"sandbox","get_default_psp.jl"))

atoms = Atoms(ext_xyz_file="atoms.xyz")
pspfiles = get_default_psp(atoms)
```

Alternatively, one can set `pspfiles` manually because it is simply
an array of `String`:
```julia
pspfiles = ["Al-q3.gth", "O-q6.gth"]
```

**IMPORTANT** Be careful to set the order of species to be same as
`atoms.SpeciesSymbols`. For example, if
```julia
atoms.SpeciesSymbols = ["Al", "O", "H"]
```
then
```julia
pspfiles = ["Al-q3.gth", "O-q6.gth", "H-q1.gth"]
```




# Calculating electron density

Several ways:

```julia
Rhoe = calc_rhoe( Nelectrons, pw, Focc, psiks, Nspin )
# or
Rhoe = calc_rhoe( Ham, psiks )
# or
calc_rhoe!( Ham, psiks, Rhoe )
```

# Subspace rotation

In case need sorting:

```julia
Hr = psiks[ikspin]' * op_H( Ham, psiks[ikspin] )
evals, evecs = eigen(Hr)
evals = real(evals[:])

# Sort in ascending order based on evals 
idx_sorted = sortperm(evals)

# Copy to Hamiltonian
Ham.electrons.ebands[:,ikspin] = evals[idx_sorted]

# and rotate
psiks[ikspin] = psiks[ikspin]*evecs[:,idx_sorted]
```


