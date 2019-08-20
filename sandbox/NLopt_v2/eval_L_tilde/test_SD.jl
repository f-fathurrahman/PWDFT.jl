include("eval_L_tilde.jl")

function test_SD()
    Random.seed!(1234)

    #Ham = create_Ham_atom_Pt_smearing()
    Ham = create_Ham_Al_fcc_smearing()
    evars = rand_ElectronicVars(Ham)

    g_evars = copy(evars)

    Ham.energies.NN = calc_E_NN( Ham.atoms )
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    Etot_old = eval_L_tilde!( Ham, evars )

    α_t = 1e-5
    β_t = 1e-1

    for iter = 1:50
        
        rotate_evars!( Ham, evars )
        grad_eval_L_tilde!( Ham, evars, g_evars )

        axpy!( -α_t, -β_t, evars, g_evars )

        rotate_evars!( Ham, evars )
        Etot = eval_L_tilde!( Ham, evars )
        print_Haux(evars, "evars after eval_L_tilde!")

        @printf("Iteration %8d %18.10f %18.10e\n", iter, Etot, Etot_old - Etot)
        #print_ebands( Ham )

        Etot_old = Etot
    end
end
@time test_SD()