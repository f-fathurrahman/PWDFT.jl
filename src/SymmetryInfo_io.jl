import Base: show
function show( io::IO, sym_info::SymmetryInfo )
    s = sym_info.s
    inv_s = sym_info.inv_s
    ft = sym_info.ft
    non_symmorphic = sym_info.non_symmorphic
    
    @printf(io, "Nsyms = %d\n", sym_info.Nsyms)
    
    for isym = 1:sym_info.Nsyms
        @printf(io, "\nSymmetry element %2d\n", isym)
        @printf(io, "Rotation matrix\n")
        for i = 1:3
            @printf(io, "%2d %2d %2d\n", s[i,1,isym], s[i,2,isym], s[i,3,isym])
        end
        @printf(io, "Translation: %13.10f %13.10f %13.10f\n", ft[1,isym], ft[2,isym], ft[3,isym])
        println(io, "non_symmorphic:  ", non_symmorphic[isym])
    end
    println(io, "No of non_symmorphic translations = ", count(non_symmorphic))
    return
end
show( sym_info::SymmetryInfo ) = show( stdout, sym_info )