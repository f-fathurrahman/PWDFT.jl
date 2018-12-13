function main( ; method="SCF" )

    # Atoms
    atoms = Atoms(xyz_string="""
            1
    
            C  0.0  0.0  0.0
            """, LatVecs=gen_lattice_sc(16.0))

    # Initialize Hamiltonian
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, ecutwfc_Ry*0.5, extra_states=3 )

    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="rpulay", betamix=0.1, use_smearing=true )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        println("ERROR: unknown method = ", method)

    end

    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands
    
    println("\nBand energies:")
    for ist = 1:Nstates
        @printf("%8d  %18.10f = %18.10f eV\n", ist, ebands[ist], ebands[ist]*Ry2eV*2)
    end
    
    println("\nTotal energy components")
    println(Ham.energies)

end
