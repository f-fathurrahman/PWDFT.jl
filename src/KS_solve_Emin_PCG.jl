"""
Solves Kohn-Sham problem using direct energy minimization as described
by Ismail-Beigi and Arias.
"""
function KS_solve_Emin_PCG!( Ham::Hamiltonian;
                             startingwfc=nothing, savewfc=false,
                             α_t=3e-5, NiterMax=200, verbose=true,
                             I_CG_BETA=2, ETOT_CONV_THR=1e-6 )

    pw = Ham.pw
    electrons = Ham.electrons
    Focc = electrons.Focc
    Nstates = electrons.Nstates
    Ns = pw.Ns
    Npoints = prod(Ns)
    CellVolume = pw.CellVolume
    Ngw = pw.gvecw.Ngw
    Ngwx = pw.gvecw.Ngwx
    Nkpt = pw.gvecw.kpoints.Nkpt
    Nspin = electrons.Nspin
    Nkspin = Nkpt*Nspin

    #
    # Initial wave function
    #
    if startingwfc == nothing
        psiks = rand_BlochWavefunc(pw, electrons)
    else
        psiks = startingwfc
    end

    #
    # Calculated electron density from this wave function and
    # update Hamiltonian (calculate Hartree and XC potential).
    #
    Rhoe = zeros(Float64,Npoints,Nspin)
    for ispin = 1:Nspin
        idxset = (Nkpt*(ispin-1)+1):(Nkpt*ispin)
        Rhoe[:,ispin] = calc_rhoe( pw, Focc[:,idxset], psiks[idxset] )
    end
    update!(Ham, Rhoe)


    #
    # Variables for PCG
    #
    g = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    d = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    g_old = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    d_old = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    Kg = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    Kg_old = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    psic = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    gt = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    #
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        g[ikspin] = zeros(ComplexF64, Ngw[ik], Nstates)
        d[ikspin] = zeros(ComplexF64, Ngw[ik], Nstates)
        g_old[ikspin] = zeros(ComplexF64, Ngw[ik], Nstates)
        d_old[ikspin] = zeros(ComplexF64, Ngw[ik], Nstates)
        Kg[ikspin] = zeros(ComplexF64, Ngw[ik], Nstates)
        Kg_old[ikspin] = zeros(ComplexF64, Ngw[ik], Nstates)
        psic[ikspin] = zeros(ComplexF64, Ngw[ik], Nstates)
        gt[ikspin] = zeros(ComplexF64, Ngw[ik], Nstates)
    end
    end
    
    β = zeros(Nkspin)
    α = zeros(Nkspin)

    Etot_old = 0.0

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    # calculate PspCore energy
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots, CellVolume )

    # Calculate energy at this psi
    energies = calc_energies(Ham, psiks)
    Ham.energies = energies
    Etot = sum(energies)

    CONVERGED = 0

    if verbose
        @printf("\n")
        @printf("Minimizing Kohn-Sham energy using PCG\n")
        @printf("-------------------------------------\n")
        @printf("NiterMax  = %d\n", NiterMax)
        @printf("α_t       = %e\n", α_t)
        @printf("conv_trh  = %e\n", ETOT_CONV_THR)
        if I_CG_BETA == 1
            @printf("Using Fletcher-Reeves formula for CG_BETA\n")
        elseif I_CG_BETA == 2
            @printf("Using Polak-Ribiere formula for CG_BETA\n")
        elseif I_CG_BETA == 3
            @printf("Using Hestenes-Stiefeld formula for CG_BETA\n")
        else
            @printf("Using Dai-Yuan formula for CG_BETA\n")
        end
        @printf("\n")
    end


    for iter = 1:NiterMax

        for ispin = 1:Nspin
        for ik = 1:Nkpt

            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin - 1)*Nkpt

            g[ikspin] = calc_grad( Ham, psiks[ikspin] )
            Kg[ikspin] = Kprec( Ham.ik, pw, g[ikspin] )

            # XXX: define function trace for real(sum(conj(...)))
            if iter != 1
                if I_CG_BETA == 1
                    β[ikspin] =
                    real(sum(conj(g[ikspin]).*Kg[ikspin]))/real(sum(conj(g_old[ikspin]).*Kg_old[ikspin]))
                elseif I_CG_BETA == 2
                    β[ikspin] =
                    real(sum(conj(g[ikspin]-g_old[ikspin]).*Kg[ikspin]))/real(sum(conj(g_old[ikspin]).*Kg_old[ikspin]))
                elseif I_CG_BETA == 3
                    β[ikspin] =
                    real(sum(conj(g[ikspin]-g_old[ikspin]).*Kg[ikspin]))/real(sum(conj(g[ikspin]-g_old[ikspin]).*d[ikspin]))
                else
                    β[ikspin] =
                    real(sum(conj(g[ikspin]).*Kg[ikspin]))/real(sum((g[ikspin]-g_old[ikspin]).*conj(d_old[ikspin])))
                end
            end
            if β[ikspin] < 0.0
                β[ikspin] = 0.0
            end

            d[ikspin] = -Kg[ikspin] + β[ikspin] * d_old[ikspin]

            psic[ikspin] = ortho_sqrt(psiks[ikspin] + α_t*d[ikspin])
        end # ik
        end # ispin
        #
        for ispin = 1:Nspin
            idxset = (Nkpt*(ispin-1)+1):(Nkpt*ispin)
            Rhoe[:,ispin] = calc_rhoe( pw, Focc[:,idxset], psiks[idxset] )
        end
        #
        update!(Ham, Rhoe)

        for ispin = 1:Nspin
        for ik = 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin - 1)*Nkpt
            gt[ikspin] = calc_grad(Ham, psic[ikspin])

            denum = real(sum(conj(g[ikspin]-gt[ikspin]).*d[ikspin]))
            if denum != 0.0
                α[ikspin] = abs( α_t*real(sum(conj(g[ikspin]).*d[ikspin]))/denum )
            else
                α[ikspin] = 0.0
            end

            # Update wavefunction
            psiks[ikspin] = psiks[ikspin] + α[ikspin]*d[ikspin]

            # Update potentials
            psiks[ikspin] = ortho_sqrt(psiks[ikspin])
        end
        end

        for ispin = 1:Nspin
            idxset = (Nkpt*(ispin-1)+1):(Nkpt*ispin)
            Rhoe[:,ispin] = calc_rhoe( pw, Focc[:,idxset], psiks[idxset] )
        end

        update!(Ham, Rhoe)

        Ham.energies = calc_energies( Ham, psiks )
        Etot = sum(Ham.energies)
        diff = abs(Etot-Etot_old)

        if verbose
            @printf("CG step %8d = %18.10f %10.7e\n", iter, Etot, diff)
        end
        
        if diff < ETOT_CONV_THR
            CONVERGED = CONVERGED + 1
        else
            CONVERGED = 0
        end

        if CONVERGED >= 2
            if verbose
                @printf("CONVERGENCE ACHIEVED\n")
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
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        psiks[ikspin] = ortho_sqrt(psiks[ikspin])
        Hr = psiks[ikspin]' * op_H( Ham, psiks[ikspin] )
        evals, evecs = eigen(Hr)
        evals = real(evals[:])
        # We need to sort this
        idx_sorted = sortperm(evals)
        Ham.electrons.ebands[:,ikspin] = evals[idx_sorted]
        psiks[ikspin] = psiks[ikspin]*evecs[:,idx_sorted]
    end
    end

    if savewfc
        for ikspin = 1:Nkpt*Nspin
            wfc_file = open("WFC_ikspin_"*string(ikspin)*".data","w")
            write( wfc_file, psiks[ikspin] )
            close( wfc_file )
        end
    end


    return

end
