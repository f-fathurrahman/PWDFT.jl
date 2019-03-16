import Base: println
function println( Ham::Hamiltonian; header=true )
    if header
        @printf("\n")
        @printf("                                  -----------\n")
        @printf("                                  Hamiltonian\n")
        @printf("                                  -----------\n")
        @printf("\n")
    end
    @printf("size (MiB) = %18.5f\n", Base.summarysize(Ham)/1024/1024)
    println("")
    println("xcfunc     = ", Ham.xcfunc)
    println("")
    println(Ham.atoms)
    println(Ham.pw)
    println(Ham.pw.gvecw.kpoints)
    println(Ham.electrons)
    for isp = 1:Ham.atoms.Nspecies
        println(Ham.pspots[isp])
    end
end