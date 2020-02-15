using LibSymspg

function spg_get_symmetry( atoms::Atoms; symprec=1e-5 )

    lattice = Matrix(atoms.LatVecs')
    positions = Matrix(inv(atoms.LatVecs))*atoms.positions # convert to fractional coordinates

    rots, trans = LibSymspg.get_symmetry(lattice, positions, atoms.atm2species, symprec)

    return size(trans)[2], rots, trans
end


function spg_get_ir_reciprocal_mesh(
             atoms::Atoms, meshk::Array{Int64,1}, is_shift::Array{Int64,1};
             is_time_reversal=1, symprec=1e-5
    )

    lattice = transpose_m3x3(atoms.LatVecs)
    positions = inv_m3x3(atoms.LatVecs)*atoms.positions # convert to fractional coordinates

    cmeshk = Base.cconvert( Array{Int32,1}, meshk )
    cis_shift = Base.cconvert( Array{Int32,1}, is_shift )
    ctypes = Base.cconvert( Array{Int32,1}, atoms.atm2species)
    num_atom = Base.cconvert( Int32, atoms.Natoms )
    is_t_rev = Base.cconvert( Int32, is_time_reversal )

    # Prepare for output
    Nkpts = prod(meshk)
    kgrid = zeros(Int32,3,Nkpts)
    mapping = zeros(Int32,Nkpts)

    num_ir =
    ccall((:spg_get_ir_reciprocal_mesh, LibSymspg.libsymspg), Int32,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Int32, Ptr{Float64}, Ptr{Float64}, 
           Ptr{Int32}, Int32, Float64),
           kgrid, mapping, cmeshk, cis_shift, is_t_rev,
           lattice, positions, ctypes, num_atom, symprec)
    
    return Base.cconvert(Int64, num_ir),
           Base.cconvert(Array{Int64,2}, kgrid),
           Base.cconvert(Array{Int64,1}, mapping)

end