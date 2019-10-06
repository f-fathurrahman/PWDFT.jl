

function precond_grad!( Ham, g, Kg; Kscalar=1.0 )
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        Kg.psiks[ikspin] = Kprec( ik, Ham.pw, g.psiks[ikspin] )
        Kg.Haux[ikspin] = Kscalar*g.Haux[ikspin]
    end
    return
end

function calc_beta_CG!( g, g_old, Kg, Kg_old, β, β_Haux )

    Nkspin = length(g.psiks)

    for i in 1:Nkspin
        ss = real(sum(conj(g_old.psiks[i]).*Kg_old.psiks[i]))
        if abs(ss) >= 1e-10
            #β[i] = real(sum(conj(g.psiks[i]-g_old.psiks[i]).*Kg.psiks[i]))/ss
            β[i] = real(sum(conj(g.psiks[i]).*Kg.psiks[i]))/ss
        else
            β[i] = 0.0
        end
        if abs(β[i] - 1.0) < 0.2
            β[i] = 0.0
        end
        if β[i] < 0.0
            β[i] = 0.0
        end

        ss = real(sum(conj(g_old.Haux[i]).*Kg_old.Haux[i]))
        if abs(ss) >= 1e-10
            #β_Haux[i] = real(sum(conj(g.Haux[i]-g_old.Haux[i]).*Kg.Haux[i]))/ss
            β_Haux[i] = real(sum(conj(g.Haux[i]).*Kg.Haux[i]))/ss
        else
            β_Haux[i] = 0.0
        end
        if abs(β_Haux[i] - 1.0) < 0.2
            β_Haux[i] = 0.0
        end
        if β_Haux[i] < 0.0
            β_Haux[i] = 0.0
        end
    end
    return
end


function calc_search_dirs!( d, Kg, d_old, β, β_Haux )
    Nkspin = length(d.psiks)
    for i in 1:Nkspin
        d.psiks[i] = Kg.psiks[i] + β[i] * d_old.psiks[i]
        d.Haux[i] = Kg.Haux[i] + β_Haux[i] * d_old.Haux[i]
    end
    #rotate!(d)
    return
end

function trial_evars!( ec, e, d, α_t, α_t_aux )
    Nkspin = length(e.psiks)
    for i in 1:Nkspin
        ec.psiks[i] = e.psiks[i] + α_t*d.psiks[i]
        ec.Haux[i] = e.Haux[i] + α_t_aux*d.Haux[i]
        ec.Haux[i] = 0.5*( ec.Haux[i] + ec.Haux[i]' )
    end
    return
end

function calc_alpha_CG!(
    g::ElectronicVars, gt::ElectronicVars, d::ElectronicVars,
    α_t::Float64, α_t_aux::Float64,
    α, α_aux
)

    Nkspin = length(g.psiks)

    for i in 1:Nkspin

        denum = real(sum(conj(g.psiks[i]-gt.psiks[i]).*d.psiks[i]))
        #println("denum = ", denum)
        #if abs(denum) <= 1e-6
        if denum != 0.0
            α[i] = abs( α_t*real(sum(conj(g.psiks[i]).*d.psiks[i]))/denum )
        else
            α[i] = 0.0
        end

        denum_aux = real(sum(conj(g.Haux[i]-gt.Haux[i]).*d.Haux[i]))
        #println("denum_aux = ", denum_aux)
        #if abs(denum) <= 1e-6
        if denum_aux != 0.0
            α_aux[i] = abs( α_t_aux*real(sum(conj(g.Haux[i]).*d.Haux[i]))/denum_aux )
        else
            α_aux[i] = 0.0
        end

    end
    return
end

function test_CG()
    Random.seed!(12345)

    #Ham = create_Ham_atom_Pt_smearing()
    #Ham = create_Ham_Pt_fcc_smearing()
    Ham = create_Ham_Al_fcc_smearing()

    evars = rand_ElectronicVars(Ham)

    # prepare guess wavefunc
    #Npoints = prod(Ham.pw.Ns)
    #Nspin = Ham.electrons.Nspin
    #Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    #Rhoe = zeros(Float64,Npoints,Nspin)
    #@assert Nspin == 1
    #Rhoe[:,1] = guess_rhoe( Ham )
    #update!(Ham, Rhoe)
    ## eigenvalues are not needed for this case
    #Ham.electrons.ebands[:,:] = diag_LOBPCG!( Ham, evars.psiks, verbose=false, verbose_last=true, NiterMax=10 )
    ##print(evals)
#
#    #for ispin in 1:Nspin, ik in 1:Nkpt
#    #    Ham.ispin = ispin
#    #    Ham.ik = ik
#    #    i = ik + (ispin -1)*Nkpt
#    #    #evars.Haux[i] = Hermitian( evars.psiks[i]' * ( Ham * evars.psiks[i] ) )
#    #    evars.Haux[i] = diagm( 0 => Ham.electrons.ebands[:,i] )
    #end


    evarsc = copy(evars)
    g_evars = copy(evars)
    Kg_evars = copy(evars)
    gt_evars = copy(evars)
    g_old_evars = copy(evars)
    Kg_old_evars = copy(evars)
    d_evars = copy(evars)
    d_old_evars = zeros_ElectronicVars(Ham)

    Nkspin = length(evars.psiks)

    β = zeros(Float64,Nkspin)
    β_Haux = zeros(Float64,Nkspin)

    α = zeros(Float64,Nkspin)
    α_Haux = zeros(Float64,Nkspin)

    Ham.energies.NN = calc_E_NN( Ham.atoms )
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    rotate_evars!( Ham, evars )
    Etot_old = eval_L_tilde!( Ham, evars )

    println("Etot_old = ", Etot_old)

    α_t = 1e-5

    Nconverges = 0

    for iter = 1:100

        rotate_evars!( Ham, evars )
        grad_eval_L_tilde!( Ham, evars, g_evars )
        #print_Haux( evars, "evars after grad_eval_L_tilde")
        #print_Haux( g_evars, "g_evars after grad_eval_L_tilde")

        calc_primary_search_dirs!( Ham, evars, Kg_evars )
        #print_Haux( evars, "evars after calc_primary_search_dirs!")
        #print_Haux( Kg_evars, "Kg_evars after calc_primary_search_dirs!")

        if iter > 1
            calc_beta_CG!( g_evars, g_old_evars, Kg_evars, Kg_old_evars, β, β_Haux )
        end
        #println("β = ", β)
        #println("β_Haux = ", β_Haux)

        calc_search_dirs!( d_evars, Kg_evars, d_old_evars, β, β_Haux )
        #print_Haux( d_evars, "d_evars after calc_search_dirs!")

        trial_evars!( evarsc, evars, d_evars, α_t, α_t )
        #print_Haux( evarsc, "evarsc after trial_evars! and before grad_eval_L_tilde")

        rotate_evars!( Ham, evarsc )
        grad_eval_L_tilde!( Ham, evarsc, gt_evars )
        #print_Haux( evarsc, "evarsc after grad_eval_L_tilde")

        #print_Haux( g_evars, "g_evars after grad_eval_L_tilde")
        #print_Haux( gt_evars, "gt_evars after grad_eval_L_tilde")
        calc_alpha_CG!( g_evars, gt_evars, d_evars, α_t, α_t, α, α_Haux )

        println("α      = ", α)
        println("α_Haux = ", α_Haux)

        # update evars
        axpy!( α, α_Haux, evars, d_evars )

        #print_Haux( evars, "evars before eval_L_tilde")
        rotate_evars!( Ham, evars )
        Etot = eval_L_tilde!( Ham, evars )
        #Etot = eval_L_tilde_and_rotate_dirs!( Ham, evars, d_evars )
        #print_Haux( evars, "evars after eval_L_tilde")

        @printf("Iteration %8d %18.10f %18.10e\n", iter, Etot, Etot_old - Etot)
        #print_ebands( Ham )
        if abs(Etot_old - Etot) < 1e-6
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end
        if Nconverges >= 2
            @printf("\nEmin_Haux_CG is converged in iter: %d\n", iter)
            break
        end



        Etot_old = Etot
        d_old_evars = copy(d_evars)
        g_old_evars = copy(g_evars)
        Kg_old_evars = copy(Kg_evars)

    end

    print_ebands( Ham )


end
@time test_CG()
