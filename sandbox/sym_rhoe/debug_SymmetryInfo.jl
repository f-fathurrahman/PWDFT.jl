struct SymmetryInfo
    Nsyms::Int64
    s::Array{Int64,3}
    inv_s::Array{Int64,3}
    ft::Array{Float64,2}
end

function SymmetryInfo( atoms::Atoms )
    
    Nsyms, s, ft = spg_get_symmetry(atoms)
    
    inv_s = zeros(Int64,3,3,Nsyms)

    for isym = 1:Nsyms
        inv_s[:,:,isym] = Base.convert(Array{Int64,2}, inv(s[:,:,isym]))
    end

    return SymmetryInfo(Nsyms, s, inv_s, ft)
end


import Base: println

function println( sym_info::SymmetryInfo )

    s = sym_info.s
    inv_s = sym_info.inv_s
    ft = sym_info.ft
    
    @printf("Nsyms = %d\n", sym_info.Nsyms)
    
    for isym = 1:sym_info.Nsyms
        @printf("\nSymmetry element %2d\n", isym)
        @printf("Rotation matrix\n")
        for i = 1:3
            @printf("%2d %2d %2d\n", s[i,1,isym], s[i,2,isym], s[i,3,isym])
        end
        @printf("Translation: %13.10f %13.10f %13.10f\n", ft[1,isym], ft[2,isym], ft[3,isym])
    end
    return
end


function spg_get_symmetry( atoms::Atoms; symprec=1e-5 )

    lattice = Matrix(atoms.LatVecs')
    positions = Matrix(inv(atoms.LatVecs))*atoms.positions # convert to fractional coordinates
    ctypes = Base.cconvert( Array{Int32,1}, atoms.atm2species)
    num_atom = Base.cconvert( Int32, atoms.Natoms )

    max_size = 50
    cmax_size = Base.cconvert(Int32, max_size)
    out_rot = zeros(Int32,3,3,max_size)
    out_translations = zeros(Float64,3,max_size)

    Nsyms_ =
    ccall((:spg_get_symmetry, PWDFT.LIBSYMSPG), Int32,
          (Ptr{Int32}, Ptr{Float64}, Int32,
           Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Int32, Float64),
           out_rot, out_translations, cmax_size,
           lattice, positions, ctypes, num_atom, symprec)
    
    Nsyms = Base.convert(Int64, Nsyms_)

    s = Base.convert(Array{Int64,3}, out_rot[:,:,1:Nsyms])
    ft = out_translations[:,1:Nsyms]

    return Nsyms, s, ft
end
