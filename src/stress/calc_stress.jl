function calc_stress!(Ham, psiks, stress)

    fill!(stress, 0.0)

    stress_hartree = zeros(Float64, 3, 3)
    calc_stress_hartree!( Ham.pw, Ham.rhoe, Ham.energies.Hartree, stress_hartree )
    stress .+= stress_hartree

    stress_Ps_loc = zeros(Float64, 3, 3)
    calc_stress_Ps_loc!( Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe, stress_Ps_loc )
    stress .+= stress_Ps_loc

    stress_xc = zeros(Float64, 3, 3)
    calc_stress_xc!( Ham.pw, Ham.potentials, Ham.rhoe, Ham.energies.XC, stress_xc )
    stress .+= stress_xc

    stress_nlcc = zeros(Float64, 3, 3)
    calc_stress_nlcc!( Ham.atoms, Ham.pspots, Ham.pw, Ham.xc_calc, Ham.xcfunc,
        Ham.rhoe, Ham.rhoe_core, stress_nlcc )
    stress .+= stress_nlcc

    stress_kin = zeros(Float64, 3, 3)
    calc_stress_kinetic!( Ham.pw, Ham.electrons, psiks, stress_kin )
    symmetrize_matrix!(Ham.pw.LatVecs, Ham.sym_info, stress_kin)
    stress .+= stress_kin

    stress_Ps_nloc = zeros(Float64, 3, 3)
    calc_stress_Ps_nloc!( Ham.atoms, Ham.pw, Ham.pspots, Ham.pspotNL, Ham.electrons,
        Ham.potentials, psiks, stress_Ps_nloc )
    symmetrize_matrix!(Ham.pw.LatVecs, Ham.sym_info, stress_Ps_nloc)
    stress .+= stress_Ps_nloc
    
    return
end

function calc_stress(Ham, psiks)
    stress = zeros(Float64, 3, 3)
    calc_stress!(Ham, psiks, stress)
    return stress
end
