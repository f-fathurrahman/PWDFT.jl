using PWDFT

include("init_V_coulomb_G.jl")
include("calc_strfact.jl")
include("PWHamiltonian.jl")
include("ortho_gram_schmidt.jl")
include("calc_rhoe.jl")
include("calc_energies.jl")
include("calc_E_NN.jl")
include("calc_grad.jl")
include("KS_solve_Emin_PCG.jl")

function test_main()
    const LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 40.0*0.5
    pw = PWGrid( ecutwfc_Ry, LatVecs )
    println(pw)

    #
    # Atoms
    #
    atoms = init_atoms_xyz("H.xyz")
    println(atoms)

    #
    # Structure factor
    #
    strf = calc_strfact( atoms, pw )

    #
    # Initialize Hamiltonian
    #
    Ham = PWHamiltonian(pw)
    Ham.potentials.Ps_loc = init_V_coulomb_G( pw, strf, [1.0] )
    println("sum V Ps loc = ", sum(Ham.potentials.Ps_loc))

    #
    # calculate E_NN
    #
    Ham.energies.NN = calc_E_NN( pw, strf, atoms.positions, atoms.Nspecies, atoms.atm2species, [1.0])

    println("\nAfter calculating E_NN")
    println(Ham.energies)

    # states
    Nstates = 1
    Ham.focc = [1.0]

    #
    KS_solve_Emin_PCG!( Ham, Nstates )
    println("\nAfter calling KS_solve_Emin_PCG:")
    println(Ham.energies)

end

@time test_main()
