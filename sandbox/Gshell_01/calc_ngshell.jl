using PWDFT

function main()
    ecutwfc = 15.0
    LatVecs = gen_lattice_sc(16.0)
    pw = PWGrid(ecutwfc, LatVecs)

    println(pw)

    init_Gshells(pw.gvec)
end

function init_Gshells(gvec::GVectors)

    Ng = gvec.Ng

    eps8 = 1e-8

    G2_sorted = sort(gvec.G2)
    ngl = 1
    for ig = 2:Ng
        if G2_sorted[ig] > (G2_sorted[ig-1] + eps8)
            ngl = ngl + 1
        end
    end
    println("ngl = ", ngl)
end

main()

