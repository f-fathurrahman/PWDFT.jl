using PWDFT

function test_main()

    atoms = init_atoms_xyz("../structures/H.xyz")
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 30.0
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, LatVecs )

    Ham.energies.NN = calc_E_NN( Ham.pw, atoms, [1.0] )

    λ, v = KS_solve_Emin_PCG!( Ham, NiterMax=5, verbose=true )

    f = open("WFC.data","w")
    write(f,v)
    close(f)

    Ngwx = Ham.pw.gvecw.Ngwx
    Nstates = Ham.electrons.Nstates
    f = open("WFC.data","r")
    psi0 = read(f,Complex128,(Ngwx,Nstates))

    λ2, v2 = KS_solve_DCM!( Ham, startingwfc=psi0 )

end

test_main()