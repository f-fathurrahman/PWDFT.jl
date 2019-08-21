include("eval_L_tilde.jl")

function test_SD()
    Random.seed!(1234)

    #Ham = create_Ham_atom_Pt_smearing()
    Ham = create_Ham_Al_fcc_smearing()
    #Ham = create_Ham_Pt_fcc_smearing()

    evars = rand_ElectronicVars(Ham)

    g_evars = copy(evars)

    Ham.energies.NN = calc_E_NN( Ham.atoms )
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    # prepare guess wavefunc
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Rhoe = zeros(Float64,Npoints,Nspin)
    @assert Nspin == 1
    Rhoe[:,1] = guess_rhoe( Ham )
    update!(Ham, Rhoe)
    # eigenvalues are not needed for this case
    _ = diag_LOBPCG!( Ham, evars.psiks, verbose=false, verbose_last=false, NiterMax=20 )

    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin -1)*Nkpt
        evars.Haux[i] = Hermitian( evars.psiks[i]' * ( Ham * evars.psiks[i] ) )
    end


    rotate_evars!( Ham, evars )
    Etot_old = eval_L_tilde!( Ham, evars )

    α_t = 1e-5
    β_t = 1e-1  # not good

    Nconverges = 0
    for iter = 1:200
        
        rotate_evars!( Ham, evars )
        grad_eval_L_tilde!( Ham, evars, g_evars )

        axpy!( -α_t, -α_t, evars, g_evars )

        rotate_evars!( Ham, evars )
        Etot = eval_L_tilde!( Ham, evars )
        print_Haux(evars, "evars after eval_L_tilde!")

        @printf("Iteration %8d %18.10f %18.10e\n", iter, Etot, Etot_old - Etot)
        if abs(Etot_old - Etot) < 1e-6
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end
        if Nconverges >= 2
            @printf("\nEmin_Haux_SD is converged in iter: %d\n", iter)
            break
        end

        Etot_old = Etot
    end

    print_ebands( Ham )
    println( Ham.energies )
end
@time test_SD()