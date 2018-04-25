function reduce_atoms( atoms::Atoms; symprec=1e-5 )

    lattice = copy(atoms.LatVecs)
    positions = inv(lattice)*copy(atoms.positions) # convert to fractional coordinates

    num_atom = Base.cconvert( Int32, atoms.Natoms )
    types = Base.cconvert(Array{Int32,1}, atoms.atm2species)

    num_primitive_atom =
    ccall( (:spg_find_primitive,SPGLIB_SO_PATH), Int32,
           ( Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Int32, Float64 ),
           lattice, positions, types, num_atom, symprec )

    # Reduced Atoms
    println("num_primitive_atom = ", num_primitive_atom)
    Natoms = Base.cconvert( Int64, num_primitive_atom )
    println("Natoms = ", Natoms)
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

    println("lattice after spg_find_primitive:")
    println(lattice)

    println("positions after spg_find_primitive:")
    println(positions)
    
    println("types after spg_find_primitive:")
    println(types)

    return Base.cconvert(Int64,num_primitive_atom)

end
