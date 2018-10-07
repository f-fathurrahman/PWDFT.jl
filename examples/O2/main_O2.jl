function main( ; method="SCF" )

    # Atoms
    atoms = Atoms( xyz_file="../structures/O2.xyz",
                   LatVecs = gen_lattice_cubic(16.0) )

    # Initialize Hamiltonian
    ecutwfc_Ry = 30.0
    pspfiles = ["../pseudopotentials/pade_gth/O-q6.gth"]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5 )
    println(Ham)

    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="anderson", betamix=0.2 )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        error( @sprintf("Unknown method = %s", method) )
    end

end
