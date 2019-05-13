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