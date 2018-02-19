using PWDFT
using PWDFT.PW

mutable struct Wavefun
    data::Array{Complex128,2}
    desc::Ref{PWGrid}
end

function Wavefun( data, pw::PWGrid )
    desc = Ref(pw)
    return Wavefun( data, desc )
end

function add_Wavefun(a::Wavefun, b::Wavefun)
    return Wavefun( a.data + b.data, a.desc )
end

function get_memory_MiB(obj)
    mem = Base.summarysize(obj)/1024./1024.
    return mem
end

function test_main()
    ecutwfc_Ry = 30.0*0.5
    LatVecs = 16.0*diagm( ones(3) )
    pw = PWGrid( ecutwfc_Ry, LatVecs, verbose=true )

    Nstates = 4
    srand(1234)

    data = rand( Complex128, (pw.gvecw.Ngwx, Nstates) )
    psi1 = Wavefun( data, pw )
    #data = 0

    #data
    #psi2 = Wavefun()

    @printf("pw   = %18.10f MiB\n", get_memory_MiB(pw))
    @printf("data = %18.10f MiB\n", get_memory_MiB(data))
    @printf("psi1 = %18.10f MiB\n", get_memory_MiB(psi1))

    println("psi1: ", psi1.desc.x.Ns)

    println("Pass here")
end

test_main()
