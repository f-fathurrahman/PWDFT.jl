include("eval_L_tilde.jl")

function calc_alpha_CG!(
    g::ElectronicVars, gt::ElectronicVars, d::ElectronicVars,
    α_t::Float64,
    α::Vector{Float64},
    α_Haux::Vector{Float64}
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
            α_Haux[i] = abs( α_t*real(sum(conj(g.Haux[i]).*d.Haux[i]))/denum_aux )
        else
            α_Haux[i] = 0.0
        end

    end
    return
end

function trial_evars!( ec, e, d, α_t )
    Nkspin = length(e.psiks)
    for i in 1:Nkspin
        ec.psiks[i] = e.psiks[i] + α_t*d.psiks[i]
        ec.Haux[i] = e.Haux[i] + α_t*d.Haux[i]
        ec.Haux[i] = 0.5*( ec.Haux[i] + ec.Haux[i]' )
    end
    return
end

function test_linmin()
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

    α_t = 3e-5

    Nkspin = length(evars.psiks)
    α = zeros(Nkspin)
    α_Haux = zeros(Nkspin)

    evarsc = copy(evars)
    d_evars = copy(evars)
    gt_evars = copy(evars)
    Kg_evars = copy(evars)
    d_old_evars = copy(evars)

    v1_psiks = zeros(Nkspin)
    v2_psiks = zeros(Nkspin)
    v3_psiks = zeros(Nkspin)

    v1_Haux = zeros(Nkspin)
    v2_Haux = zeros(Nkspin)
    v3_Haux = zeros(Nkspin)

    Nconverges = 0
    for iter = 1:100
        
        #rotate_evars!( Ham, evars )
        grad_eval_L_tilde!( Ham, evars, g_evars )

        if iter > 1
            v1_psiks, v1_Haux = dot( g_evars, d_evars )
            v2_psiks, v2_Haux = dot( g_evars, g_evars )
            v3_psiks, v3_Haux = dot( d_evars, d_evars )
            for i in 1:Nkspin
                @printf("cosine angle psiks = %18.10f, ", v1_psiks[i]/sqrt(v2_psiks[i]*v3_psiks[i]))
                @printf("cosine angle Haux = %18.10f\n", v1_Haux[i]/sqrt(v2_Haux[i]*v3_Haux[i]))
            end
        end
        #print_Haux( evars, "evars after eval_L_tilde")
        #print_Haux( g_evars, "g_evars after eval_L_tilde")

        #calc_primary_search_dirs!( Ham, evars, Kg_evars )
        #d_evars = copy( Kg_evars )

        d_evars = ElectronicVars( -g_evars.psiks, -g_evars.Haux )

        trial_evars!( evarsc, evars, d_evars, α_t )

        rotate_evars!( Ham, evarsc )
        grad_eval_L_tilde!( Ham, evarsc, gt_evars )

        calc_alpha_CG!( g_evars, gt_evars, d_evars, α_t, α, α_Haux )

        #println("α      = ", α)
        #println("α_Haux = ", α_Haux)

        # update evars
        axpy!( α, α_Haux, evars, d_evars )

        rotate_evars!( Ham, evars )
        Etot = eval_L_tilde!( Ham, evars )
        #print_Haux(evars, "evars after eval_L_tilde!")

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
        d_old_evars = copy(d_evars)
    end

    print_ebands( Ham )
    println( Ham.energies )
end
@time test_linmin()