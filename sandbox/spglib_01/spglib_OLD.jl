
# This function is now included in the KPoints constructor
function gen_kgrid_reduced( atoms::Atoms, mesh::Array{Int64,1}, is_shift::Array{Int64,1};
                            time_reversal=1 )
    
    @printf("\n")
    @printf("Generating kpoints:\n")
    @printf("mesh     = (%d,%d,%d)\n", mesh[1], mesh[2], mesh[3])
    @printf("is_shift = (%d,%d,%d)\n", is_shift[1], is_shift[2], is_shift[3])

    num_ir, kgrid, mapping =
    spg_get_ir_reciprocal_mesh( atoms, mesh, is_shift, is_time_reversal=time_reversal )

    @printf("\n")
    @printf("Number of kpoints = %d\n", num_ir)
    @printf("\n")

    umap = unique(mapping)

    Nkpt = prod(mesh)

    list_ir_k = []
    for ikk = 1:num_ir
        for ik = 1:Nkpt
            if umap[ikk] == mapping[ik]
                append!( list_ir_k, [kgrid[:,ik]] )
                break
            end
        end
    end

    kred = zeros(Float64,3,num_ir)
    for ik = 1:num_ir
        kred[1,ik] = list_ir_k[ik][1] / mesh[1]
        kred[2,ik] = list_ir_k[ik][2] / mesh[2]
        kred[3,ik] = list_ir_k[ik][3] / mesh[3]
    end
    
    # prepare for
    kcount = zeros(Int64,num_ir)
    for ik = 1:num_ir
        kcount[ik] = count( i -> ( i == umap[ik] ), mapping )
    end

    # calculate the weights
    wk = kcount[:]/sum(kcount)

    # convert to cartesian
    RecVecs = 2*pi*invTrans_m3x3(atoms.LatVecs)
    kred = RecVecs*kred
    return kred, wk

end



"""
This function try to reduce number of atoms by exploting crystal symmetry.
"""
function reduce_atoms( atoms::Atoms; symprec=1e-5 )

    #lattice = copy(atoms.LatVecs)'
    lattice = transpose_m3x3(atoms.LatVecs)
    positions = inv_m3x3(atoms.LatVecs)*atoms.positions # convert to fractional coordinates

    num_atom = Base.cconvert( Int32, atoms.Natoms )
    types = Base.cconvert(Array{Int32,1}, atoms.atm2species)

    num_primitive_atom =
    ccall( (:spg_find_primitive, LibSymspg.libsymspg), Int32,
           ( Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Int32, Float64 ),
           lattice, positions, types, num_atom, symprec )

    # Prepare for reduced Atoms
    Natoms = Base.cconvert( Int64, num_primitive_atom )

    # Transpose back
    LatVecs = transpose_m3x3(lattice)
    positions = LatVecs*positions[:,1:num_primitive_atom]
    atm2species = Base.cconvert( Array{Int64,1}, types[1:num_primitive_atom] )
    Nspecies = atoms.Nspecies
    SpeciesSymbols = atoms.SpeciesSymbols
    atsymbs = Array{String}(undef,Natoms)
    for ia = 1:Natoms
        isp = atm2species[ia]
        atsymbs[ia] = SpeciesSymbols[isp]
    end
    Zvals = zeros(Nspecies)
    return Atoms( Natoms, Nspecies, positions, atm2species, atsymbs, SpeciesSymbols, LatVecs, Zvals )

end


function spg_find_primitive( atoms::Atoms; symprec=1e-5)

    # We need to transpose lattice
    # For positions we don't need to transpose it.

    #lattice = copy(atoms.LatVecs)'
    lattice = transpose_m3x3(atoms.LatVecs)
    positions = inv_m3x3(atoms.LatVecs)*atoms.positions # convert to fractional coordinates

    num_atom = Base.cconvert( Int32, atoms.Natoms )
    types = Base.cconvert(Array{Int32,1}, atoms.atm2species)

    num_primitive_atom =
    ccall( (:spg_find_primitive, LibSymspg.libsymspg), Int32,
           ( Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Int32, Float64 ),
           lattice, positions, types, num_atom, symprec )

    return Base.cconvert(Int64,num_primitive_atom)

end


function spg_get_ir_reciprocal_mesh(
             atoms::Atoms, meshk::Array{Int64,1}, is_shift::Array{Int64,1};
             is_time_reversal=1, symprec=1e-5
         )

    #lattice = copy(atoms.LatVecs)'
    lattice = PWDFT.transpose_m3x3(atoms.LatVecs)
    positions = PWDFT.inv_m3x3(atoms.LatVecs)*atoms.positions # convert to fractional coordinates

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

