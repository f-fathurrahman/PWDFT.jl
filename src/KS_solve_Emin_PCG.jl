function KS_solve_Emin_PCG!( Ham::Hamiltonian; kwargs... )
    KS_solve_Emin_PCG!( Ham, rand_BlochWavefunc(Ham); kwargs... )
    return
end

# Not yet used
function constrain_search_dir!( d::BlochWavefunc, psiks::BlochWavefunc )
    Nkspin = length(psiks)
    for i in 1:Nkspin
        d[i] = d[i] - psiks[i] * ( psiks[i]' * d[i] )
    end
    return
end

function constrain_search_dir!( d::Array{ComplexF64,2}, psi::Array{ComplexF64,2} )
    d[:] = d - psi * ( psi' * d )
    return
end


"""
    KS_solve_Emin_PCG!( Ham, kwargs... )

Solve Kohn-Sham problem using direct energy minimization as described
by Ismail-Beigi and Arias.
"""
function KS_solve_Emin_PCG!(
    Ham::Hamiltonian, psiks::BlochWavefunc;
    startingrhoe::Symbol=:gaussian,
    skip_initial_diag::Bool=false,
    α_t::Float64=3e-5,
    NiterMax::Int64=200,
    verbose::Bool=true,
    print_final_ebands::Bool=false,
    print_final_energies::Bool=true,
    i_cg_beta::Int64=2,
    etot_conv_thr::Float64=1e-6,
    α_max::Float64=2.1,
    restrict_linmin::Bool=false
)

    pw = Ham.pw
    electrons = Ham.electrons
    
    Focc = electrons.Focc
    Nstates = electrons.Nstates
    Nelectrons = electrons.Nelectrons
    
    Ns = pw.Ns
    Npoints = prod(Ns)
    CellVolume = pw.CellVolume
    Ngw = pw.gvecw.Ngw
    Ngwx = pw.gvecw.Ngwx
    Nkpt = pw.gvecw.kpoints.Nkpt
    
    Nspin = electrons.Nspin_wf
    Nkspin = Nkpt*Nspin

    #
    # Calculated electron density from this wave function and
    # update Hamiltonian (calculate Hartree and XC potential).
    #
    Rhoe = zeros(Float64,Npoints,Nspin)

    if startingrhoe == :gaussian
        @assert Nspin == 1
        Rhoe[:,1] = guess_rhoe( Ham )
        #Rhoe = guess_rhoe_atomic( Ham ) # when smearing is ready this should be preferred
    else
        calc_rhoe!( Ham, psiks, Rhoe )
    end

    update!(Ham, psiks, Rhoe)

    evals = zeros(Nstates,Nkspin)
    if !skip_initial_diag
        evals =
        diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false, NiterMax=10 )
    end

    #
    # Variables for PCG
    #
    g      = zeros_BlochWavefunc( Ham )
    d      = deepcopy(g)
    g_old  = deepcopy(g)
    d_old  = deepcopy(g)
    Kg     = deepcopy(g)
    Kg_old = deepcopy(g)
    psic   = deepcopy(g)
    gt     = deepcopy(g)
    
    β = zeros(Nkspin)
    α = zeros(Nkspin)

    Etot_old = 0.0

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    # Calculate energy at this psi
    energies = calc_energies(Ham, psiks)
    Ham.energies = energies
    Etot = sum(energies)

    Nconverges = 0

    if verbose
        @printf("\n")
        @printf("Minimizing Kohn-Sham energy using PCG\n")
        @printf("-------------------------------------\n")
        @printf("NiterMax  = %d\n", NiterMax)
        @printf("α_t       = %e\n", α_t)
        @printf("conv_thr  = %e\n", etot_conv_thr)
        if i_cg_beta == 1
            @printf("Using Fletcher-Reeves formula for CG_BETA\n")
        elseif i_cg_beta == 2
            @printf("Using Polak-Ribiere formula for CG_BETA\n")
        elseif i_cg_beta == 3
            @printf("Using Hestenes-Stiefeld formula for CG_BETA\n")
        else
            @printf("Using Dai-Yuan formula for CG_BETA\n")
        end
        @printf("\n")
    end


    for iter in 1:NiterMax

        for ispin in 1:Nspin, ik in 1:Nkpt

            Ham.ik = ik
            Ham.ispin = ispin
            i = ik + (ispin - 1)*Nkpt

            g[i] = calc_grad( Ham, psiks[i] )
            Kg[i] = Kprec( Ham.ik, pw, g[i] )

            # XXX: define function trace for real(sum(conj(...)))
            if iter != 1
                if i_cg_beta == 1
                    β[i] =
                    real(sum(conj(g[i]).*Kg[i]))/real(sum(conj(g_old[i]).*Kg_old[i]))
                    if abs(β[i] - 1.0) < 0.2
                        β[i] = 0.0
                    end
                elseif i_cg_beta == 2
                    β[i] = real( dot(g[i]-g_old[i], Kg[i]) )/real( dot(g_old[i],Kg_old[i]) )
                elseif i_cg_beta == 3
                    β[i] =
                    real(sum(conj(g[i]-g_old[ik]).*Kg[i]))/real(sum(conj(g[i]-g_old[i]).*d[i]))
                else
                    β[i] =
                    real(sum(conj(g[i]).*Kg[i]))/real(sum((g[i]-g_old[i]).*conj(d_old[i])))
                    if abs(β[i] - 1.0) < 0.2
                        β[i] = 0.0
                    end                    
                end
            end
            if β[i] < 0.0 || isnan(β[i])
                #println("Resetting β")
                β[i] = 0.0
            end

            d[i] = -Kg[i] + β[i] * d_old[i]

            psic[i] = ortho_sqrt(psiks[i] + α_t*d[i])
        end # ik, ispin
        
        calc_rhoe!( Ham, psic, Rhoe )

        update!(Ham, psic, Rhoe)

        for ispin in 1:Nspin, ik in 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            i = ik + (ispin - 1)*Nkpt
            gt[i] = calc_grad(Ham, psic[i])

            denum = real(sum(conj(g[i]-gt[i]).*d[i]))
            if denum != 0.0
                α[i] = abs( α_t*real(sum(conj(g[i]).*d[i]))/denum )
            else
                α[i] = 0.0
            end

            if restrict_linmin
                if α[i] > α_max
                    @printf("α for ikspin #%d is too large: %f, restrict it to %f\n", i, α[i], α_max)
                    α[i] = α_max
                end
            end

            # Update wavefunction
            psiks[i] = psiks[i] + α[i]*d[i]
            psiks[i] = ortho_sqrt(psiks[i])
        end

        # Update Rhoe and potentials
        calc_rhoe!( Ham, psiks, Rhoe )

        update!(Ham, psiks, Rhoe)

        Ham.energies = calc_energies( Ham, psiks )
        Etot = sum(Ham.energies)
        diffE = Etot_old - Etot

        if verbose
            @printf("Emin_PCG step %4d = %18.10f %16.6e\n", iter, Etot, diffE)
        end
        
        if verbose && (diffE < 0.0) && (iter > 1)
            println("*** WARNING: Etot is not decreasing")
        end
        
        if abs(diffE) < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        if Nconverges >= 2
            if verbose
                @printf("\nEmin_PCG is converged in iter: %d\n", iter)
            end
            break
        end

        g_old = copy(g)
        d_old = copy(d)
        Kg_old = copy(Kg)
        Etot_old = Etot

        flush(stdout)
    end

    # Calculate eigenvalues
    evecs = zeros(ComplexF64,Nstates,Nstates)
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        i = ik + (ispin - 1)*Nkpt
        psiks[i] = ortho_sqrt(psiks[i])
        Hr = Hermitian(psiks[i]' * op_H(Ham, psiks[i]))
        evals, evecs = eigen(Hr)
        Ham.electrons.ebands[:,i] = evals
        psiks[i] = psiks[i]*evecs
    end

    if verbose && print_final_ebands
        @printf("\n")
        @printf("----------------------------\n")
        @printf("Final Kohn-Sham eigenvalues:\n")
        @printf("----------------------------\n")
        @printf("\n")
        print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints, unit="eV")
    end

    if verbose && print_final_energies
        @printf("\n")
        @printf("-------------------------\n")
        @printf("Final Kohn-Sham energies:\n")
        @printf("-------------------------\n")
        @printf("\n")
        println(Ham.energies)
    end

    return

end
