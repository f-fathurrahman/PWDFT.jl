using PWDFT

function test_read()
    psp = PsPot_GTH("../pseudopotentials/pade_gth/Pt-q10.gth")
    println(psp)
    psp = PsPot_GTH("../pseudopotentials/pade_gth/Pt-q18.gth")
    println(psp)
end


import PyPlot
const plt = PyPlot

function test_main()
    LatVecs = 16*diagm(ones(3))
    pw = PWGrid(15.0, LatVecs)
    println(pw)

    G2 = sort(pw.gvec.G2)
    Ω = pw.Ω

    psp = PsPot_GTH("../pseudopotentials/pade_gth/Ni-q18.gth")
    println(psp)
    Vg = eval_Vloc_G(psp, G2, Ω)

    plt.clf()
    plt.plot( G2, Vg, marker="o" )
    plt.xlim(0.1,1)
    plt.savefig("Ni_V_ps_loc_G.png", dpi=200)
end

#test_read()
test_main()
