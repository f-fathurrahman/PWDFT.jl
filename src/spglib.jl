using LibSymspg

function spg_get_symmetry( atoms::Atoms; symprec=1e-5 )

    lattice = Matrix(atoms.LatVecs')
    positions = Matrix(inv(atoms.LatVecs))*atoms.positions # convert to fractional coordinates

    ctypes = Base.cconvert( Array{Int32,1}, atoms.atm2species)
    num_atom = Base.cconvert( Int32, atoms.Natoms )

    max_size = 48*num_atom
    cmax_size = Base.cconvert(Int32, max_size)
    out_rot = zeros(Int32,3,3,max_size)
    out_translations = zeros(Float64,3,max_size)

    Nsyms =
    ccall((:spg_get_symmetry, LibSymspg.libsymspg), Int32,
          (Ptr{Int32}, Ptr{Float64}, Int32,
           Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Int32, Float64),
           out_rot, out_translations, cmax_size,
           lattice, positions, ctypes, num_atom, symprec)

    # return minus of out_translations to match QE convention.
    return Nsyms, out_rot[:,:,1:Nsyms], -out_translations[:,1:Nsyms]

    #rots, trans = LibSymspg.get_symmetry(lattice, positions, atoms.atm2species, symprec)
    #return size(trans)[2], rots, trans
end


function spg_get_ir_reciprocal_mesh(
             atoms::Atoms, meshk::Array{Int64,1}, is_shift::Array{Int64,1};
             is_time_reversal=1, symprec=1e-5
    )

    lattice = Matrix(atoms.LatVecs')
    positions = Matrix(inv(atoms.LatVecs))*atoms.positions # convert to fractional coordinates

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