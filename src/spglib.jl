function spg_find_primitive( atoms::Atoms; symprec=1e-5)

    lattice = atoms.LatVecs
    positions = inv(lattice)*atoms.positions # convert to fractional coordinates

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
