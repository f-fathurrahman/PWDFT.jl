import Base: show
function show( io::IO, sym_info::SymmetryInfo )
    
    @printf(io, "Nsyms = %d\n", sym_info.Nsyms)
    @printf(io, "Nrots = %d\n", sym_info.Nrots)

    s = sym_info.s
    sname = sym_info.sname
    ft = sym_info.ft

    for isym in 1:sym_info.Nsyms
        @printf(io, "\nSymmetry operation #%d: %s\n", isym, sname[isym])
        for i in 1:3
            @printf(io, "%4d %4d %4d\n", s[i,1,isym], s[i,2,isym], s[i,3,isym])
        end
        @printf(io, "Fractional translation: %10.4f %10.4f %10.4f\n",
            sym_info.ft[1,isym], sym_info.ft[2,isym], sym_info.ft[3,isym])
    end
end

show( sym_info::SymmetryInfo ) = show( stdout, sym_info )
