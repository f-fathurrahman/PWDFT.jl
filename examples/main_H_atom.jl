using PWDFT

function test_main( ; method="SCF" )

    #
    # Atoms
    #
    atoms = init_atoms_xyz("H.xyz")
    println(atoms)

    #
    # Initialize Hamiltonian
    #
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 30.0
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, LatVecs )

    #
    # calculate E_NN
    #
    Ham.energies.NN = calc_E_NN( Ham.pw, atoms, [1.0] )

    #
    # Solve the KS problem
    #
    if method == "SCF"
        KS_solve_SCF!( Ham )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham, verbose=true )

    elseif method == "DCM"
        KS_solve_DCM!( Ham )

    else
        println("ERROR: unknow method = ", method)
    end

    #Nstates = Ham.electrons.Nstates
    #println("\nEigenvalues")
    #for ist = 1:Nstates
    #    @printf("%8d  %18.10f = %18.10f eV\n", ist, λ[ist], λ[ist]*Ry2eV*2)
    #end
    #println("\nTotal energy components")
    #println(Ham.energies)

end

val, t, bytes, gctime, memallocs = @timed test_main(method="Emin")
#val, t, bytes, gctime, memallocs = @timed test_main(method="SCF")
#val, t, bytes, gctime, memallocs = @timed test_main(method="DCM")

GiB = 1024.0*1024.0*1024.0
@printf("Allocated memory  = %f GiB\n",bytes/GiB)
@printf("Elapsed time = %f s\n",t)
println("malloc = ",memallocs.malloc)
println("poolalloc = ",memallocs.poolalloc)
