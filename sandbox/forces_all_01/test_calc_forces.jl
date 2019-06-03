using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("../forces_ewald_01/calc_forces_NN.jl")
include("../forces_local_01/calc_forces_Ps_loc.jl")
include("../forces_nonlocal_01/calc_forces_Ps_nloc.jl")
include("symmetry_atoms.jl")


function symmetrize_vector!(pw::PWGrid, sym_info::SymmetryInfo, irt, v::Array{Float64,2})

    Nsyms = sym_info.Nsyms
    LatVecs = pw.LatVecs
    RecVecs = pw.RecVecs
    s = sym_info.s

    if Nsyms == 1
        return
    end

    Nvecs = size(v)[2]

    tmp = zeros(3,Nvecs)
    
    # bring vector to crystal axis
    for i = 1:Nvecs
        tmp[:,i] = v[1,i]*LatVecs[1,:] + v[2,i]*LatVecs[2,:] + v[3,i]*LatVecs[3,:]
    end
    
    println(tmp)

    # symmetrize in crystal axis
    v[:,:] .= 0.0
    for i = 1:Nvecs
        for isym = 1:Nsyms
            iar = irt[isym,i]
            v[:,i] = v[:,i] + s[:,1,isym]*tmp[1,iar]
                            + s[:,2,isym]*tmp[2,iar]
                            + s[:,3,isym]*tmp[3,iar]
            println(isym, " ", v[:,i])
        end
        #println(v[:,i])
    end
    
    tmp[:,:] = v[:,:]/Nsyms
    
    # bring vector back to cartesian axis
    for i = 1:Nvecs
        v[:,i] = tmp[1,i]*RecVecs[:,1] + tmp[2,i]*RecVecs[:,2] + tmp[3,i]*RecVecs[:,3]
    end

    v = v/(2*pi)

    return
end


include("create_Ham.jl")


function test_main()

    #Ham = create_Ham_H2()

    Ham = create_Ham_Si_fcc()
    #Ham = create_Ham_GaAs_v1()
    #Ham = create_Ham_GaAs_v2()

    irt = init_irt(Ham.atoms, Ham.sym_info)

    Random.seed!(1234)
    
    #KS_solve_Emin_PCG!(Ham, etot_conv_thr=1e-8, savewfc=true)
    #KS_solve_SCF!(Ham, mix_method="rpulay", etot_conv_thr=1e-8, savewfc=true)

    psiks = read_psiks(Ham)

    for psi in psiks
        ortho_check(psi)
    end

    atoms = Ham.atoms
    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs

    F_NN = calc_forces_NN( Ham.atoms )*2.0
    println("NN forces:")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_NN[1,ia], F_NN[2,ia], F_NN[3,ia] )
    end

    F_Ps_loc = calc_forces_Ps_loc( Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe )*2.0
    println("Ps loc forces:")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_loc[1,ia], F_Ps_loc[2,ia], F_Ps_loc[3,ia] )
    end

    F_Ps_nloc =
    calc_forces_Ps_nloc( Ham.atoms, Ham.pw, Ham.pspots, Ham.electrons, Ham.pspotNL, psiks )*2.0
    println("Ps nloc forces:")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_nloc[1,ia], F_Ps_nloc[2,ia], F_Ps_nloc[3,ia] )
    end


    F_total = F_NN + F_Ps_loc + F_Ps_nloc
    println("Total forces:")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_total[1,ia], F_total[2,ia], F_total[3,ia] )
    end

    @printf("Sum of forces in x-dir: %18.10f\n", sum(F_total[1,:]))
    @printf("Sum of forces in y-dir: %18.10f\n", sum(F_total[2,:]))
    @printf("Sum of forces in z-dir: %18.10f\n", sum(F_total[3,:]))

    symmetrize_vector!( Ham.pw, Ham.sym_info, irt, F_total )
    println("Total forces symmetrized:")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_total[1,ia], F_total[2,ia], F_total[3,ia] )
    end


end


test_main()
