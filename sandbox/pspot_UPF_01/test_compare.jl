import LightXML
using Printf

using PWDFT

include("PsPot_UPF.jl")

function main()
    psp = PsPot_GTH("../../pseudopotentials/pade_gth/Si-q4.gth")
    println(psp)

    psp_upf = PsPot_UPF("Si.pz-hgh.UPF")
    println(psp_upf)

    idx_r = 2
    r = psp_upf.r[idx_r]
    @printf("r = %18.10f\n", r)
    @printf("V_local = %18.10f\n", psp_upf.V_local[idx_r])
    @printf("eval_V_local = %18.10f\n", PWDFT.eval_Vloc_R(psp, r))

    proj_func = psp_upf.proj_func
    Nproj = psp_upf.Nproj
    iprj = 1
    #for iprj in 1:Nproj
        l = psp_upf.proj_l[iprj]
        
        upf_proj = proj_func[idx_r,iprj]/r
        println("upf_proj = ", upf_proj)
        
        prj_analytic = PWDFT.eval_proj_R(psp, l, 1, r)
        println("eval_proj_R = ", prj_analytic)

        println("ratio = ", prj_analytic/upf_proj)
    #end
    println("Pass here")
end

main()