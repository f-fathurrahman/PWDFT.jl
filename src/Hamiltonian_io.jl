import Base: print
function print( io::IO, Ham::Hamiltonian; header=true )
    if header
        @printf("\n")
        @printf("                                  -----------\n")
        @printf("                                  Hamiltonian\n")
        @printf("                                  -----------\n")
        @printf("\n")
    end
    @printf(io, "size (MiB) = %18.5f\n", Base.summarysize(Ham)/1024/1024)
    println(io, "")
    println(io, "xcfunc       = ", Ham.xcfunc)
    # FIXME: this is only for Libxc
    #println(io, "xc_calc.x_id = ", Ham.xc_calc.x_id)
    #println(io, "xc_calc.c_id = ", Ham.xc_calc.c_id)
    println(io, "")
    print(io, Ham.atoms)
    print(io, Ham.pw)
    print(io, Ham.pw.gvecw.kpoints)
    print(io, Ham.electrons)
    for isp = 1:Ham.atoms.Nspecies
        print(io, Ham.pspots[isp])
    end
end
print( Ham::Hamiltonian; header=true ) = print( stdout, Ham, header=header )