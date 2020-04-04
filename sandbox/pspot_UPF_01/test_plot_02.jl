using DelimitedFiles
import LightXML
using Printf

using PWDFT

import PyPlot
const plt = PyPlot

include("PsPot_UPF.jl")

function test_main()

    psp = PsPot_GTH("../../pseudopotentials/pade_gth/Si-q4.gth")
    println(psp)

    psp_upf = PsPot_UPF("Si.pz-hgh.UPF")
    println(psp_upf)

    Nr = psp_upf.Nr
    r = psp_upf.r

    plt.clf()
    plt.plot(psp_upf.r, psp_upf.V_local, label="V_local")
    plt.xlim(0.0, 5.0)
    plt.legend()
    plt.savefig("IMG_Si_Vloc.pdf")

    Nproj = psp_upf.Nproj
    Vproj = zeros(Float64,Nr,Nproj)
    #
    l = 0; iprj = 1
    for i in 1:Nr
        Vproj[i,1] = PWDFT.eval_proj_R(psp, l, iprj, r[i])
    end
    #
    l = 0; iprj = 2
    for i in 1:Nr
        Vproj[i,2] = PWDFT.eval_proj_R(psp, l, iprj, r[i])
    end
    #
    l = 1; iprj = 1
    for i in 1:Nr
        Vproj[i,3] = PWDFT.eval_proj_R(psp, l, iprj, r[i])
    end


    plt.clf()
    plt.plot(psp_upf.r, psp_upf.proj_func[:,1]./r, label="upf proj1")
    plt.plot(psp_upf.r, Vproj[:,1], label="Vproj")
    plt.xlim(0.0, 4.0)
    plt.legend()
    plt.savefig("IMG_Si_proj_func.pdf")

    plt.clf()
    plt.plot(psp_upf.r, psp_upf.proj_func[:,2]./r, label="upf proj2")
    plt.plot(psp_upf.r, Vproj[:,2], label="Vproj2")
    plt.xlim(0.0, 4.0)
    plt.legend()
    plt.savefig("IMG_Si_proj_func_2.pdf")

    plt.clf()
    plt.plot(psp_upf.r, psp_upf.proj_func[:,3]./r, label="upf proj3")
    plt.plot(psp_upf.r, Vproj[:,3], label="Vproj3")
    plt.xlim(0.0, 4.0)
    plt.legend()
    plt.savefig("IMG_Si_proj_func_3.pdf")


end

test_main()
