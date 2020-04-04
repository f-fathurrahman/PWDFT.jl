import LightXML
using Printf

using PWDFT

include("PsPot_UPF.jl")

function main()
    psp = PsPot_GTH("../../pseudopotentials/pade_gth/Si-q4.gth")
    println(psp)

    psp_upf = PsPot_UPF("Si.pz-hgh.UPF")
    println(psp_upf)
end

main()