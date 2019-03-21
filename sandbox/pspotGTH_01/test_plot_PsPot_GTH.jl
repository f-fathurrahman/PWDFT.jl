using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

import PyPlot
const plt = PyPlot

function test_main()
    LatVecs = gen_lattice_sc(16.0)
    pw = PWGrid(15.0, LatVecs)
    println(pw)

    G2 = sort(pw.gvec.G2)
    CellVolume = pw.CellVolume

    psp = PsPot_GTH(joinpath(DIR_PSP,"Ni-q18.gth"))
    println(psp)
    Vg = eval_Vloc_G(psp, G2)

    plt.clf()
    plt.plot( G2, Vg, marker="o" )
    plt.xlim(0.1,1)
    plt.savefig("TEMP_Ni_V_ps_loc_G.png", dpi=200)
end

test_main()