mutable struct RotationsCache
    rotPrev::Vector{Matrix{ComplexF64}}
    rotPrevC::Vector{Matrix{ComplexF64}}
    rotPrevCinv::Vector{Matrix{ComplexF64}}
    Urot::Vector{Matrix{ComplexF64}}
    UrotC::Vector{Matrix{ComplexF64}}
end

function RotationsCache(Nkspin, Nstates)
    rotPrev = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    rotPrevC = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    rotPrevCinv = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    Urot = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    UrotC = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        rotPrev[ikspin] = Matrix(1.0*I(Nstates))
        rotPrevC[ikspin] = Matrix(1.0*I(Nstates))
        rotPrevCinv[ikspin] = Matrix(1.0*I(Nstates))
        Urot[ikspin] = Matrix(1.0*I(Nstates))
        UrotC[ikspin] = Matrix(1.0*I(Nstates))
    end
    return RotationsCache(
        rotPrev,
        rotPrevC,
        rotPrevCinv,
        Urot,
        UrotC,
    )
end


function reset_rotations!(rots_cache)
    Nkspin = size(rots_cache.rotPrev, 1)
    Nstates = size(rots_cache.rotPrev[1], 1)
    for ikspin in 1:Nkspin
        rots_cache.rotPrev[ikspin] = Matrix(1.0*I(Nstates))
        rots_cache.rotPrevC[ikspin] = Matrix(1.0*I(Nstates))
        rots_cache.rotPrevCinv[ikspin] = Matrix(1.0*I(Nstates))
        rots_cache.Urot[ikspin] = Matrix(1.0*I(Nstates))
        rots_cache.UrotC[ikspin] = Matrix(1.0*I(Nstates))
    end
    return
end


function rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)
    Nkspin = length(g)
    rotPrevCinv = rots_cache.rotPrevCinv
    rotPrev = rots_cache.rotPrev
    for ikspin in 1:Nkspin
        g[ikspin][:,:] = g[ikspin][:,:] * rotPrevCinv[ikspin]
        Kg[ikspin][:,:] = Kg[ikspin][:,:] * rotPrevCinv[ikspin]
        g_Haux[ikspin][:,:] = rotPrev[ikspin] * g_Haux[ikspin][:,:] * rotPrev[ikspin]'
        Kg_Haux[ikspin][:,:] = rotPrev[ikspin] * Kg_Haux[ikspin][:,:] * rotPrev[ikspin]'
    end
    return
end



#
# Various functions to update Hamiltonian
#

# Ham.electrons.ebands are already updated elsewhere
function update_from_ebands!(Ham)
    update_from_ebands!(Ham, Ham.electrons.ebands)
    return
end


# Input: ebands
# Modifies: Focc, E_fermi, mTS
# Ham.electrons.ebands are not modified here
function update_from_ebands!(Ham, ebands)

    if !Ham.electrons.use_smearing
        return
    end

    # NOTE: ebands are assumed to be updated outside this function

    # Calculate Kohn-Sham eigenvalues and occupation numbers
    Focc = Ham.electrons.Focc
    Nelectrons = Ham.electrons.Nelectrons
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    wk = Ham.pw.gvecw.kpoints.wk
    kT = Ham.electrons.kT
    @assert kT > 1e-10

    E_fermi, mTS = update_Focc!(
        Focc, smear_fermi, smear_fermi_entropy,
        ebands, Float64(Nelectrons), kT,
        Nkpt, wk
    )
    # Set some output
    Ham.electrons.E_fermi = E_fermi
    Ham.energies.mTS = mTS

    return
end


# Input: psiks
# Modifies: Ham.rhoe, potentials
function update_from_wavefunc!(Ham, psiks)    
    # Compute electron density from psiks
    # Use Ham.rhoe
    calc_rhoe!(Ham, psiks, Ham.rhoe)
    # Update the potentials
    update_from_rhoe!(Ham, psiks, Ham.rhoe)
    # XXX: update_from_rhoe! will not overwrite update Ham.rhoe
    return
end

function get_diag_Haux_from_ebands( Ham )
    Nspin = Ham.electrons.Nspin_channel
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Haux[ikspin] = diagm( 0 => Ham.electrons.ebands[:,ikspin] )
    end
    return Haux
end

# Modifies Ham.electrons.ebands, save rotations in rots_cache
# psiks is assumed to be
function transform_psiks_Haux_update_ebands!(
    Ham, psiks, Haux, rots_cache;
    overwrite_Haux=true,
    do_ortho_psi=true
)
    Nspin = Ham.electrons.Nspin_channel
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    ebands = Ham.electrons.ebands
    #
    Urot = rots_cache.Urot
    UrotC = rots_cache.UrotC
    rotPrev = rots_cache.rotPrev
    rotPrevC = rots_cache.rotPrevC
    rotPrevCinv = rots_cache.rotPrevCinv
    #
    for ispin in 1:Nspin, ik in 1:Nkpt
        # Don't forget to set current index for Hamiltonian
        # Need this in case we need to apply op_S
        Ham.ispin = ispin
        Ham.ik = ik
        ikspin = ik + (ispin-1)*Nkpt
        #
        ebands[:,ikspin], Urot[ikspin] = eigen(Hermitian(Haux[ikspin]))
        if overwrite_Haux
            Haux[ikspin] = diagm( 0 => Ham.electrons.ebands[:,ikspin] )
        end
        #
        # XXX Check this again, probably do_ortho_psi should always be true
        if do_ortho_psi
            if Ham.need_overlap
                UrotC[ikspin] = inv(sqrt(psiks[ikspin]' * op_S(Ham, psiks[ikspin])))
            else
                UrotC[ikspin] = inv(sqrt(psiks[ikspin]' * psiks[ikspin]))
            end
        end
        UrotC[ikspin] = UrotC[ikspin]*Urot[ikspin] # extra rotation
        psiks[ikspin] = psiks[ikspin]*UrotC[ikspin]
    end
    #
    for ikspin in 1:Nkspin
        # Save previous
        rotPrev[ikspin] = rotPrev[ikspin] * Urot[ikspin]
        rotPrevC[ikspin] = rotPrevC[ikspin] * UrotC[ikspin]
        rotPrevCinv[ikspin] = inv(UrotC[ikspin]) * rotPrevCinv[ikspin]
    end
    return
end



# Haux (in diagonal form) is stored in Ham.electrons.ebands
# psiks is already orthonormalized and rotated according to Urot
# that makes Haux diagonal
# 
# The following must be 
# update_from_ebands!( Ham, ebands )
# update_from_wavefunc!( Ham, psiks )
#
function calc_Lfunc(
    Ham::Hamiltonian,
    psiks::BlochWavefunc
)
    calc_energies!(Ham, psiks)
    # get entropy
    # Ham.energies.mTS is computed in update_from_ebands!
    return sum(Ham.energies)
end

function do_step_psiks_Haux!(α::Float64, Ham, psiks, Haux, d, d_Haux, rots_cache)
    Nkspin = length(psiks)
    rotPrev = rots_cache.rotPrev
    rotPrevC = rots_cache.rotPrevC
    # Step
    for ikspin in 1:Nkspin
        psiks[ikspin] += α * d[ikspin] * rotPrevC[ikspin]
        Haux[ikspin]  += α * rotPrev[ikspin]' * d_Haux[ikspin] * rotPrev[ikspin]
    end

    transform_psiks_Haux_update_ebands!( Ham, psiks, Haux, rots_cache, do_ortho_psi=true )
    return
end



function do_compute_energy(Ham, psiks)
    # Update Hamiltonian terms
    update_from_ebands!( Ham )
    update_from_wavefunc!( Ham, psiks )
    # Now, we are ready to evaluate
    E = calc_Lfunc( Ham, psiks )
    return E
end



function linmin_quad_v01!(
    α_t,
    Ham, psiks, Haux, Hsub, g, g_Haux,
    Kg, Kg_Haux,
    d, d_Haux, rots_cache,
    E_old
)

    gd = 2*real(dot(g,d)) + real(dot(g_Haux, d_Haux))
    println("gd = $(gd)")
    α_prev = 0.0
    if gd > 0
        println("ERROR: Bad step direction")
        return E_old, false, α_prev
    end

    NtryMax = 5

    α_safe =  1e-5 # safe step size
    α = α_safe # declare this outside for loop, set to a "safe" value
    α_t_ReduceFactor = 0.1
    α_t_IncreaseFactor = 3.0
    is_success = false
    for itry in 1:NtryMax
        println("--- Begin itry linmin trial step = $(itry) using α_t=$(α_t)")
        #
        do_step_psiks_Haux!(α_t - α_prev, Ham, psiks, Haux, d, d_Haux, rots_cache)
        # make explicit calls to update_* functions
        #
        α_prev = α_t
        E_t = do_compute_energy(Ham, psiks) # this will update ebands, Focc, and Rhoe
        #
        if !isfinite(E_t)
            α_t *= α_t_ReduceFactor
            println("α_t is reduced to=$(α_t)")
            continue # continue
        end
        # prediciton of step size
        c = ( E_t - (E_old + α_t*gd) ) / α_t^2
        α = -gd/(2*c)
        if α < 0
            println("Wrong curvature, α is negative: E_t=$(E_t), E_old=$(E_old)")
            α_t *= α_t_IncreaseFactor
            println("Trial step will become true step. α_t will be set to $(α_t)")
            # calculate gradients
            calc_grad_psiks!(Ham, psiks, g, Kg, Hsub)
            calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
            #@infiltrate
            # return trial energy and status
            is_success = true
            return E_t, is_success, α_t
        end
        break
    end
    println("Find α = $(α)")
    
    # actual step
    for itry in 1:NtryMax
        #
        println("--- Begin itry linmin actual step = $(itry) using α=$(α)")
        #
        do_step_psiks_Haux!(α - α_prev, Ham, psiks, Haux, d, d_Haux, rots_cache)
        α_prev = α
        # calculate energy and gradients
        E_t2 = do_compute_energy(Ham, psiks)
        calc_grad_psiks!(Ham, psiks, g, Kg, Hsub)
        calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
        #
        println("Actual step energy 2: E_t2 = $(E_t2)")
        #
        if !isfinite(E_t2)
            α *= α_t_ReduceFactor
            println("α is reduced to=$(α)")
        end
        # trial energy is higher, reduce α
        if E_t2 > E_old
            α *= α_t_ReduceFactor
            println("Energy is not decreasing, try do decrease α to $(α)")
            continue # continue iteration
        else
            println("Actual step is successful")
            is_success = true
            return E_t2, is_success, α
        end
    end

    # default is unsuccessful try
    return Inf, false, α

end


function constrain_search_dir!( Ham, d::BlochWavefunc, psiks::BlochWavefunc )
    Nspin = Ham.electrons.Nspin_channel
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    #XXX need overlap ?
    for ispin in 1:Nspin, ik in 1:Nkpt
        # Don't forget to set current index for Hamiltonian
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin-1)*Nkpt
        if Ham.need_overlap
            d[i][:,:] = d[i] - psiks[i] * ( psiks[i]' * op_S(Ham, d[i]) )
        else
            d[i][:,:] = d[i] - psiks[i] * ( psiks[i]' * d[i] )
        end
    end
    return
end


function my_Kprec!(Ham, g, Kg)
    Nspin = Ham.electrons.Nspin_channel
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        Kprec!( ik, Ham.pw, g[ikspin], Kg[ikspin] )
    end
    return
end


# for psiks, renamed to calc_grad_psiks! to avoid name clashing
# input: Ham, psiks
# output: g, Hsub
function calc_grad_psiks!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    g::BlochWavefunc, Kg::BlochWavefunc,
    Hsub::Vector{Matrix{ComplexF64}}
)
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    Nspin = Ham.electrons.Nspin_channel
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Focc = Ham.electrons.Focc
    wk = Ham.pw.gvecw.kpoints.wk

    #XXX Preallocate Hpsi ? Using Ngwx?
    # Ngwx = maximum(Ham.pw.gvecw.Ngw)

    for ispin in 1:Nspin, ik in 1:Nkpt
        # Don't forget to set current index for Hamiltonian
        Ham.ispin = ispin
        Ham.ik = ik
        ikspin = ik + (ispin-1)*Nkpt
        #
        Hpsi = op_H( Ham, psiks[ikspin] )
        Hsub[ikspin][:,:] = psiks[ikspin]' * Hpsi
        if Ham.need_overlap
            Spsi = op_S(Ham, psiks[ikspin])
            Hpsi[:,:] -= Spsi * Hsub[ikspin]
        else
            Hpsi[:,:] -= psiks[ikspin] * Hsub[ikspin] # op_S(psiks[iskspin])
        end
        for ist in 1:Nstates
            # dont forget Focc and wk factor 
            g[ikspin][:,ist] .= Focc[ist,ikspin] .* Hpsi[:,ist] * wk[ik]
            # FIXME: for nonspin-pol?
        end
        Kprec!( ik, Ham.pw, Hpsi, Kg[ikspin] )
    end
    return
end


# Gradient for Haux
# The real input is actually stored in Ham.electrons.ebands which
# is calculated from diagonalizing Haux
# Haux need to be diagonal here
# Input: Ham, Hsub
# Output: g_Haux, Kg_Haux
function calc_grad_Haux!(
    Ham, Hsub, g_Haux, Kg_Haux;
    κ=1.0
)
    # κ is a scalar multiplier for Kg_Haux. There is some heuristics mentioned
    # in the original paper about how to tune this, however it is not implemented
    # yet

    Nspin = Ham.electrons.Nspin_channel
    Nstates = Ham.electrons.Nstates
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    fprime = zeros(Float64, Nstates)
    fprimeNum = zeros(Float64, Nstates)
    dmuNum = zeros(Float64, Nspin)
    dmuDen = zeros(Float64, Nspin)

    # These variables are not updated or calculated here
    # They are assumed to be calculated elsewhere
    ebands = Ham.electrons.ebands
    kT = Ham.electrons.kT
    E_fermi = Ham.electrons.E_fermi
    wk = Ham.pw.gvecw.kpoints.wk

    if Nspin == 1
        w_spin = 2.0
    else
        w_spin = 1.0 
    end
    for ispin in 1:Nspin
        dmuNum[ispin] = 0.0
        dmuDen[ispin] = 0.0
        for ik in 1:Nkpt
            ikspin = ik + (ispin-1)*Nkpt
            # accumulate with Nkpt? Add wk ?
            for ist in 1:Nstates
                fprime[ist] = smear_fermi_prime( ebands[ist,ikspin], E_fermi, kT )
                fprimeNum[ist] = fprime[ist] * ( real(Hsub[ikspin][ist,ist]) - ebands[ist,ikspin] )
            end
            # smear_fermi_prime might return NaN if E_fermi is not set properly
            dmuNum[ispin] += wk[ik] * sum(fprimeNum) * w_spin
            dmuDen[ispin] += wk[ik] * sum(fprime) * w_spin
        end
        #println("ispin=$(ispin) dmu = $(dmuNum[ispin]) $(dmuDen[ispin])")
    end

    dmuContrib = sum(dmuNum)/sum(dmuDen)
    if isnan(dmuContrib)
        @warn "dmuContrib is problematic"
        sign_frac = sign(sum(dmuNum))*sign(sum(dmuDen))
        if sign_frac == 0
            dmuContrib = 0.0 #1
        else
            dmuContrib = 0.0 #1*sign_frac
        end
        println("dmuContrib = $(dmuContrib)")
    end
    #dBzContrib = 0.0 # not used

    gradF0 = zeros(ComplexF64, Nstates, Nstates)
    gradF = zeros(ComplexF64, Nstates, Nstates)
    g_tmp = zeros(ComplexF64, Nstates, Nstates)

    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        gradF0[:,:] = Hsub[ikspin] - diagm( 0 => ebands[:,ikspin] )
        gradF[:,:] = copy(gradF0)
        for ist in 1:Nstates
            gradF[ist,ist] = gradF0[ist,ist] - dmuContrib # FIXME: not tested for spinpol
        end
        g_tmp[:,:] = grad_smear( smear_fermi, smear_fermi_prime, ebands[:,ikspin], E_fermi, kT, gradF )
        g_Haux[ikspin][:,:] = wk[ik] * 0.5 * (g_tmp' + g_tmp) * w_spin
        Kg_Haux[ikspin][:,:] = -κ * gradF0[:,:] # preconditioning here?
    end

    return

end


# Probably useful for refactoring

#=
# Eq (24)
function calc_dFdmu(
    Hsub::Matrix{Float64},
    ebands,
    Focc_in,
    kT
)
    Nstates = size(ebands, 1)
    dFdmu = zeros(Float64, Nstates)
    if Nspin == 1
        Focc = 0.5*Focc_in
    else
        Focc = Focc_in
    end
    for ist in 1:Nstates
        dFdmu[ist] = (Hsub[ist,ist] - ebands[ist,1])*Focc[ist,1]*(1 - Focc[ist,1])
    end
    return dFdmu/kT
end


function offdiag_elements( Hsub, ebands, E_f::Float64, kT::Float64 )
    Nstates = size(evals, 1)
    mat = zeros(Nstates,Nstates)
    for j in 1:Nstates, i in 1:Nstates
        de = ebands[i] - ebands[j]
        if abs(de) > 1e-6
            mat[i,j] = Hsub[i,j] * ( smear_fermi(ebands[i], E_f, kT) - smear_fermi(ebands[j], E_f, kT) ) / de
        end
    end
    return mat
end
=#


# Like `electrons_scf` but using direct minimization algorithm for metals
# TODO: Add original references (including JDFTx)
function electrons_Emin_Haux!(Ham; NiterMax=100, psiks=nothing, Haux=nothing)

    Nspin = Ham.electrons.Nspin_channel
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin
    Nstates = Ham.electrons.Nstates
    Rhoe = Ham.rhoe
    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns)

    # Initialize electronic variables: `psiks` and `Haux`:
    Random.seed!(1234)
    if isnothing(psiks)
        psiks = rand_BlochWavefunc(Ham)
    end

    Hsub = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        Hsub[ikspin] = zeros(ComplexF64, Nstates, Nstates)
    end
    # Calculate Hsub
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        ikspin = ik + (ispin-1)*Nkpt
        Hsub[ikspin][:,:] = psiks[ikspin]' * (Ham * psiks[ikspin])
    end

    # Prepare Haux (random numbers)
    #
    if isnothing(Haux)
        Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
        for ikspin in 1:Nkspin
            Haux[ikspin] = randn(ComplexF64, Nstates, Nstates)
            # the same as Hsub
            Haux[ikspin][:,:] = Hsub[ikspin][:,:]
            Haux[ikspin][:,:] = 0.5*( Haux[ikspin] + Haux[ikspin]' ) # make symmetric
            #Haux[ikspin] = diagm(0 => sort(randn(Float64, Nstates)))
        end
    end

    # Gradients, subspace Hamiltonian
    g = zeros_BlochWavefunc(Ham)
    Kg = zeros_BlochWavefunc(Ham)
    d = zeros_BlochWavefunc(Ham)
    gPrev = zeros_BlochWavefunc(Ham)
    #
    g_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    Kg_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    d_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    gPrev_Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    for ikspin in 1:Nkspin
        g_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        Kg_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        d_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        gPrev_Haux[ikspin] = zeros(ComplexF64, Nstates, Nstates)
    end

    rots_cache = RotationsCache(Nkspin, Nstates);

    # psiks is already orthonormal
    # Make Haux diagonal and rotate psiks
    # Ham.electrons.ebands are updated here
    transform_psiks_Haux_update_ebands!( Ham, psiks, Haux, rots_cache )

    # Update Hamiltonian, compute energy and gradients at current psiks and Haux:

    # Update Hamiltonian before evaluating free energy
    update_from_ebands!( Ham )
    #update_from_wavefunc!( Ham, psiks )
    E1 = calc_Lfunc( Ham, psiks )
    println("E1 = $(E1)")
    #
    # Calculate gradients
    calc_grad_psiks!(Ham, psiks, g, Kg, Hsub)
    calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
    rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)

    println("Initial Focc = ")
    display(Ham.electrons.Focc); println()
    println("Initial ebands (w.r.t Fermi) = ")
    display(Ham.electrons.ebands .- Ham.electrons.E_fermi); println()

    α_t_start = 1.0
    α_t_min = 1e-5
    α_t = α_t_start

    do_force_grad_dir = true
    gKNorm = 0.0
    gKNormPrev = 0.0
    # current and previous norms of the preconditioned gradient

    ok_paw = any(Ham.pspotNL.are_paw) # only used for printing?

    for iterCG in 1:NiterMax

        println("\nStart iterCG = ", iterCG)

        gKNorm = 2*real(dot(g, Kg)) + real(dot(g_Haux, Kg_Haux))

        β = 0.0
        if !do_force_grad_dir
            gPrevKg = 2*real(dot(gPrev, Kg)) + real(dot(gPrev_Haux, Kg_Haux))
            gd = 2*real(dot(g, d)) + real(dot(g_Haux, d_Haux))
            gg = 2*real(dot(g, g)) + real(dot(g_Haux, g_Haux))
            dd = 2*real(dot(d, d)) + real(dot(d_Haux, d_Haux))
            if gg*dd > 0
                @printf("linmin: %10.3le\n", gd/sqrt(gg*dd))
            else
                @warn "Negative gg*dd encountered"
            end
            if gKNorm*gKNormPrev > 0
                @printf("cgtest: %10.3le\n", gPrevKg/sqrt(gKNorm*gKNormPrev))
            else
                @warn "Negative gKNorm*gKNormPrev encountered"
            end
            # Update beta:
            println("gKNorm = $(gKNorm), gPrevKg = $(gPrevKg)")
            β = (gKNorm - gPrevKg)/gKNormPrev
            println("β = ", β)              
            if β < 0.0
                println("!!!! Resetting CG because β is negative")
                β = 0.0
            end
        end

        do_force_grad_dir = false

        # XXX TODO Check convergence here?

        # Save previous gradient
        gKNormPrev = gKNorm
        for ikspin in 1:Nkspin
            gPrev[ikspin][:] = g[ikspin][:]
            gPrev_Haux[ikspin][:] = g_Haux[ikspin][:]
        end

        # Set direction
        for ikspin in 1:Nkspin
            d[ikspin] = -Kg[ikspin] + β*d[ikspin]
            d_Haux[ikspin] = -Kg_Haux[ikspin] + β*d_Haux[ikspin]
        end
        constrain_search_dir!(Ham, d, psiks)

        #
        # Do line minimization:
        # XXX: Probably check gd: if it is already small then no need to do line minimization?
        E_new, is_success, α = linmin_quad_v01!(
            α_t,
            Ham, psiks, Haux, Hsub, g, g_Haux, Kg, Kg_Haux, d, d_Haux, rots_cache, E1
        )
        #println("Test grad psiks before rotate: $(2*dot(g, psiks))")
        #println("Test grad Haux before rotate: $(dot(Haux, g_Haux))")
        rotate_gradients!(g, Kg, g_Haux, Kg_Haux, rots_cache)
 
        #
        if is_success
            α_t = α
            println("linminQuad is successful. α_t is updated to α = $α")
            if α_t < α_t_min
                # bad step size: make sure next test step size is not too bad
                α_t = α_t_start 
                println("Bad step size is encountered, α_t is set to α_t_start = $(α_t_start)")
            end
        else
            println("WARN: Line minimization is not successful")
            #
            do_step_psiks_Haux!(-α, Ham, psiks, Haux, d, d_Haux, rots_cache)
            # calculate energy and gradients
            update_from_ebands!( Ham )
            update_from_wavefunc!( Ham, psiks )
            # Now, we are ready to evaluate
            E_new = Inf # calc_Lfunc( Ham, psiks )
            calc_grad_psiks!(Ham, psiks, g, Kg, Hsub)
            calc_grad_Haux!(Ham, Hsub, g_Haux, Kg_Haux)
            #
            if β > 0.0
                # Failed, but not along the gradient direction:
                println("Forcing gradient direction")
                do_force_grad_dir = true
            else
                # Failed along the gradient direction
                println("Probably round off error")
                break
            end
        end

        if Nspin == 2
            magn = sum(Rhoe[:,1] - Rhoe[:,2])*dVol
            integRhoe = sum(Rhoe)*dVol
            println("integRhoe = $integRhoe integ magn = $magn")
        else
            integRhoe = sum(Rhoe)*dVol
            println("integRhoe = $integRhoe")
        end
        ΔE = abs(E_new - E1)
        println("\niterCG: $(iterCG) E_new = $(E_new) ΔE = $(ΔE)\n")
        #println("Focc = ")
        #display(Ham.electrons.Focc); println()
        #println("ebands (w.r.t) Fermi energy = ")
        #display(Ham.electrons.ebands .- Ham.electrons.E_fermi); println()
        #println("Energies:")
        #println(Ham.energies, use_smearing=true, is_paw=ok_paw)

        if ΔE < 1e-8
            println("\nConverged !!!")
            break
        end 
        
        # New iterations, variables are updated in linmin_quad_v01
        E1 = E_new
    end

    println("Final Focc = ")
    display(Ham.electrons.Focc); println()
    println("Final ebands (w.r.t) Fermi energy = ")
    display(Ham.electrons.ebands .- Ham.electrons.E_fermi); println()
    println("Final Energies:")
    println(Ham.energies, use_smearing=true, is_paw=ok_paw)

    #@infiltrate

    return
end

