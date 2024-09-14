function calc_stress!(Ham, psiks, stress)

    fill!(stress, 0.0)

    s_tmp = zeros(Float64, 3, 3)

    stress_Ps_loc = zeros(Float64, 3, 3)
    calc_stress_Ps_loc!( Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe, stress_Ps_loc )
    stress .+= stress_Ps_loc
    #
    s_tmp .= 2.0*stress_Ps_loc # convert to Ry
    println("\nstress_Ps_loc (in Ry/bohr^3) = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", s_tmp[i,1], s_tmp[i,2], s_tmp[i,3])
    end

    stress_hartree = zeros(Float64, 3, 3)
    calc_stress_hartree!( Ham.pw, Ham.rhoe, Ham.energies.Hartree, stress_hartree )
    stress .+= stress_hartree
    #
    s_tmp .= 2.0*stress_hartree # convert to Ry
    println("\nstress_hartree (in Ry/bohr^3) = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", s_tmp[i,1], s_tmp[i,2], s_tmp[i,3])
    end

    stress_xc = zeros(Float64, 3, 3)
    calc_stress_xc!( Ham.pw, Ham.potentials, Ham.rhoe, Ham.energies.XC, stress_xc )
    stress .+= stress_xc
    #
    s_tmp .= 2.0*stress_xc # convert to Ry
    println("\nstress_xc (in Ry/bohr^3) = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", s_tmp[i,1], s_tmp[i,2], s_tmp[i,3])
    end

    stress_nlcc = zeros(Float64, 3, 3)
    calc_stress_nlcc!( Ham.atoms, Ham.pspots, Ham.pw, Ham.xc_calc, Ham.xcfunc,
        Ham.rhoe, Ham.rhoe_core, stress_nlcc )
    stress .+= stress_nlcc
    #
    s_tmp .= 2.0*stress_nlcc # convert to Ry
    println("\nstress_nlcc (in Ry/bohr^3) = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", s_tmp[i,1], s_tmp[i,2], s_tmp[i,3])
    end

    stress_NN = zeros(Float64, 3, 3)
    calc_stress_NN!( Ham.atoms, Ham.pw, stress_NN )
    stress .+= stress_NN
    #
    s_tmp .= 2.0*stress_NN # convert to Ry
    println("\nstress_NN (in Ry/bohr^3) = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", s_tmp[i,1], s_tmp[i,2], s_tmp[i,3])
    end

    stress_kin = zeros(Float64, 3, 3)
    calc_stress_kinetic!( Ham.pw, Ham.electrons, psiks, stress_kin )
    #
    s_tmp .= 2.0*stress_kin # convert to Ry
    println("\nstress_kin (in Ry/bohr^3) unsymmetrized = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", s_tmp[i,1], s_tmp[i,2], s_tmp[i,3])
    end
    symmetrize_matrix!(Ham.pw.LatVecs, Ham.sym_info, stress_kin)
    stress .+= stress_kin



    stress_Ps_nloc = zeros(Float64, 3, 3)
    calc_stress_Ps_nloc!( Ham.atoms, Ham.pw, Ham.pspots, Ham.pspotNL, Ham.electrons,
        Ham.potentials, psiks, stress_Ps_nloc )
    #
    s_tmp .= 2.0*stress_Ps_nloc # convert to Ry
    println("\nstress_Ps_nloc (in Ry/bohr^3) unsymmetrized = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", s_tmp[i,1], s_tmp[i,2], s_tmp[i,3])
    end
    #
    symmetrize_matrix!(Ham.pw.LatVecs, Ham.sym_info, stress_Ps_nloc)
    stress .+= stress_Ps_nloc

    # Symmetrize all
    symmetrize_matrix!(Ham.pw.LatVecs, Ham.sym_info, stress)

    return
end

function calc_stress(Ham, psiks)
    stress = zeros(Float64, 3, 3)
    calc_stress!(Ham, psiks, stress)
    return stress
end
