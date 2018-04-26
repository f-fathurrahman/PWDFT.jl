function reduce_atoms( atoms::Atoms; symprec=1e-5 )

    lattice = copy(atoms.LatVecs)
    positions = inv(lattice)*copy(atoms.positions) # convert to fractional coordinates

    num_atom = Base.cconvert( Int32, atoms.Natoms )
    types = Base.cconvert(Array{Int32,1}, atoms.atm2species)

    num_primitive_atom =
    ccall( (:spg_find_primitive,SPGLIB_SO_PATH), Int32,
           ( Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Int32, Float64 ),
           lattice, positions, types, num_atom, symprec )

    # Prepare for reduced Atoms
    Natoms = Base.cconvert( Int64, num_primitive_atom )
    LatVecs = lattice
    positions = lattice*positions[:,1:num_primitive_atom]
    atm2species = Base.cconvert( Array{Int64,1}, types[1:num_primitive_atom] )
    Nspecies = atoms.Nspecies
    SpeciesSymbols = atoms.SpeciesSymbols
    atsymbs = Array{String}(Natoms)
    for ia = 1:Natoms
        isp = atm2species[ia]
        atsymbs[ia] = SpeciesSymbols[isp]
    end
    Zvals = zeros(Nspecies)
    return Atoms( Natoms, Nspecies, positions, atm2species, atsymbs, SpeciesSymbols, LatVecs, Zvals )

end


function spg_find_primitive( atoms::Atoms; symprec=1e-5)

    lattice = copy(atoms.LatVecs)
    positions = inv(lattice)*copy(atoms.positions) # convert to fractional coordinates

    num_atom = Base.cconvert( Int32, atoms.Natoms )
    types = Base.cconvert(Array{Int32,1}, atoms.atm2species)

    num_primitive_atom =
    ccall( (:spg_find_primitive,SPGLIB_SO_PATH), Int32,
           ( Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Int32, Float64 ),
           lattice, positions, types, num_atom, symprec )

    return Base.cconvert(Int64,num_primitive_atom)

end


function spg_get_ir_reciprocal_mesh(
             atoms::Atoms, mesh::Array{Int64,1}, is_shift::Array{Int64};
             is_time_reversal=1, symprec=1e-5
         )

    lattice = copy(atoms.LatVecs)
    positions = inv(lattice)*copy(atoms.positions) # convert to fractional coordinates

    cmesh = Base.cconvert( Array{Cint,1}, mesh )
    cis_shift = Base.cconvert( Array{Cint,1}, is_shift )
    ctypes = Base.cconvert( Array{Cint,1}, atoms.atm2species)
    num_atom = Base.cconvert( Cint, atoms.Natoms )
    is_t_rev = Base.cconvert( Cint, is_time_reversal )
    
    # Prepare for output
    Nkpts = prod(mesh)
    kgrid = zeros(Cint,3,Nkpts)
    mapping = zeros(Cint,Nkpts)
    
    num_ir =
    ccall((:spg_get_ir_reciprocal_mesh, SPGLIB_SO_PATH), Cint,
        (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Ptr{Float64}, Ptr{Float64}, 
        Ptr{Cint}, Cint, Float64),
        kgrid, mapping, cmesh, cis_shift, is_t_rev,
        lattice, positions, ctypes, num_atom, symprec)
    
    return Base.cconvert(Int64, num_ir),
           Base.cconvert(Array{Int64,2}, kgrid),
           Base.cconvert(Array{Int64,1}, mapping)

end

