using Printf
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

    G2_shells = zeros(ngl)
    idx_gshells = zeros(Int64,Ng)
    
    G2_shells[1] = G2_sorted[1]
    idx_gshells[1] = 1

    igl = 1
    for ig = 2:Ng
        if G2_sorted[ig] > (G2_sorted[ig-1] + eps8)
            igl = igl + 1
            G2_shells[igl] = G2_sorted[ig]
        end
        idx_gshells[ig] = igl
    end

    for igl = 1:ngl
        @printf("igl = %5d, G2 = %18.10f\n", igl, G2_shells[igl])
    end

end

main()

