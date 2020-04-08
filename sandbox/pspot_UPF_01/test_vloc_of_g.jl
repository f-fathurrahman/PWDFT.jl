import LightXML
using Printf
using SpecialFunctions: erf

import PyPlot
const plt = PyPlot

using PWDFT

include("integ_simpson.jl")
include("PsPot_UPF.jl")
include("vloc_of_g.jl")

function main()
    psp = PsPot_GTH("../../pseudopotentials/pade_gth/Si-q4.gth")
    println(psp)

    psp_upf = PsPot_UPF("Si.pz-hgh.UPF")
    println(psp_upf)

    pw = PWGrid(50.0, gen_lattice_fcc(5.0))

    #G2_shells = pw.gvec.G2_shells
    #Ngl = length(G2_shells)

    Ngl = 101
    G2_shells = range(0.0, stop=10.0, length=Ngl)

    Vgl = zeros(Float64, Ngl)
    for igl in 1:Ngl
        Vgl[igl] = eval_Vloc_G( psp, G2_shells[igl] )
    end

    Vgl_upf = init_Vloc_G(
        psp_upf.Nr, psp_upf.r, psp_upf.rab,
        psp_upf.V_local, psp_upf.zval,
        Ngl, collect(G2_shells),
        pw.CellVolume
    )

    plt.clf()
    plt.plot(sqrt.(G2_shells), Vgl, marker="o", label="Vgl")
    plt.plot(sqrt.(G2_shells), Vgl_upf, marker="o", label="Vgl_upf")
    plt.legend()
    plt.xlim(0.0,sqrt(10.0))
    plt.savefig("IMG_Vgl.pdf")
end

main()