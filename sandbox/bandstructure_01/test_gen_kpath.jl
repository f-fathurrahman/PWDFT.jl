using PWDFT
using LinearAlgebra

function test_main()
    dict_spec_kpts = get_special_kpoints("hexagonal")

    kpath = "GMKG"
    Nkpt_spec = length(kpath)
    kpts_spec = zeros(3,Nkpt_spec)
    kpts_spec_label = []
    for ik = 1:Nkpt_spec
        label = string(kpath[ik])
        kvec = dict_spec_kpts[label]
        println(kpath[ik], " kvec = ", kvec)
        append!( kpts_spec_label, label )
        kpts_spec[:,ik] = kvec
    end

    # distance between two adjacent special points in kpath
    d = zeros(Nkpt_spec-1)
    for ik = 1:Nkpt_spec-1
        d[ik] = norm( kpts_spec[:,ik+1] - kpts_spec[:,ik] )
        println(ik, " ", d[ik])
    end

    # number of kpoints between two adjacent special points in kpath
    Δk = 0.02
    Nk = zeros(Int64,Nkpt_spec-1)
    for ik = 1:Nkpt_spec-1
        Nk[ik] = round(Int64, d[ik]/Δk)
        println(ik, " ", Nk[ik])
    end

    ipk = 0
    Nkpt_on_path = sum(Nk .- 1) + 1
    kpath = zeros(3,Nkpt_on_path)
    for ik = 1:Nkpt_spec-1
        dvec = kpts_spec[:,ik+1] - kpts_spec[:,ik]
        dk = dvec./(Nk[ik]-1)
        println("dk = ", dk, " norm(dk) = ", norm(dk))
        for iik = 1:Nk[ik]-1
            kvec = kpts_spec[:,ik] + (iik-1).*dk
            ipk = ipk + 1
            kpath[:,ipk] = kvec
            println(ipk, " kpath = ", kpath[:,ipk])
        end
    end

    println("sum(Nk) = ", sum(Nk))
    println("Nkpt_on_path = ", Nkpt_on_path)
    
    kpath[:,Nkpt_on_path] = kpts_spec[:,Nkpt_spec]
    
    for ik = 1:Nkpt_on_path
        println("kpath = ", kpath[:,ik])
    end
end

test_main()

