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

include("SymmetryInfo_io.jl")
