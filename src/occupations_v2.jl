struct OccupationUpdater
    smear_func::Function
    smear_func_entropy::Function
    smear_func_prime::Function
    kT::Float64
end

function OccupationUpdater(; smear_type=:FermiDirac, kT::Float64=0.01)
    if smear_type == :FermiDirac
        return OccupationUpdater(
            smear_fermi, smear_fermi_entropy, smear_fermi_prime, kT
        )
    elseif smear_type == :Gaussian
        return OccupationUpdater(
            smear_gauss, smear_gauss_entropy, smear_gauss_prime, kT
        )
    elseif (smear_type == :MarzariVanderbilt) || (smear_type == :cold)
        return OccupationUpdater(
            smear_cold, smear_cold_entropy, smear_cold_prime, kT
        )
    else
        println("Unknown smear_type: ", smear_type)
        error()
    end
end


function update_Focc!(
    occ_updater::OccupationUpdater,
    Focc::Matrix{Float64},
    ebands::Matrix{Float64},
    Nelectrons, Nkpt::Int64, wk
)
    return update_Focc!(
        Focc, occ_updater.smear_func, occ_updater.smear_func_entropy,
        ebands,
        Nelectrons, occ_updater.kT, Nkpt, wk
    )
end

# Update occupation number according to `smear_func`
# with given energy eigenvalues `ebands`
function update_Focc!(
    Focc::Array{Float64,2},
    smear_func, smear_func_entropy,
    ebands::Array{Float64,2},
    Nelectrons,
    kT::Float64,
    Nkpt::Int64,
    wk::Vector{Float64}
)

    E_f = find_E_fermi( smear_func, ebands, Nelectrons, kT, Nkpt, wk )
    #@info "Found E_fermi = $(E_f)"

    Nstates = size(ebands, 1)
    Nkspin = size(ebands, 2)
    Nspin = Int64(Nkspin/Nkpt)
    
    if Nspin == 1
        w_spin = 2.0 # weight factor
    else
        w_spin = 1.0
    end

    # Occupation numbers
    for ispin in 1:Nspin, ik in 1:Nkpt, ist in 1:Nstates
        ikspin = ik + (ispin-1)*Nkpt
        Focc[ist,ikspin] = w_spin * smear_func( ebands[ist,ikspin], E_f, kT )
    end
    
    # Entropy
    mTS = 0.0
    for ispin in 1:Nspin, ik in 1:Nkpt, ist in 1:Nstates
        ikspin = ik + (ispin-1)*Nkpt
        mTS = mTS - w_spin * wk[ik] * kT * smear_func_entropy( ebands[ist,ikspin], E_f, kT )
    end
    return E_f, mTS  
end


# FIXME: include in update_Focc! ?
function calc_electronic_entropy(
    smear_func_entropy, ebands, E_f, kT, Nkpt::Int64, wk
)
    Nstates = size(ebands, 1)
    Nkspin = size(ebands, 2)
    Nspin = Int64(Nkspin/Nkpt)
    mTS = 0.0
    if Nspin == 1
        w_spin = 2.0 # weight factor
    else
        w_spin = 1.0
    end
    for ispin in 1:Nspin, ik in 1:Nkpt, ist in 1:Nstates
        ikspin = ik + (ispin-1)*Nkpt
        mTS = mTS - w_spin * wk[ik] * kT * smear_func_entropy( ebands[ist,ikspin], E_f, kT )
    end
    return mTS
end


function sum_Focc(
    smear_func,
    ebands::Array{Float64,2},
    E_f::Float64,
    kT::Float64,
    Nkpt::Int64,
    wk
)
    Nstates = size(ebands, 1)
    Nkspin = size(ebands, 2)
    Nspin = Int64(Nkspin/Nkpt)
    #
    if Nspin == 1
        w_spin = 2.0 # weight factor
    else
        w_spin = 1.0
    end
    #
    ss = 0.0
    for ispin in 1:Nspin, ik in 1:Nkpt, ist in 1:Nstates
        ikspin = ik + (ispin-1)*Nkpt
        ss = ss + w_spin * wk[ik] * smear_func( ebands[ist,ikspin], E_f, kT )
    end
    return ss
end


function find_E_fermi(
    smear_func,
    ebands::Array{Float64,2},
    Nelectrons,
    kT::Float64,
    Nkpt::Int64,
    wk;
    NiterMax=300, verbose=false
)

    Nstates = size(ebands,1)
    Nkspin = size(ebands,2)
    Nspin = Int64(Nkspin/Nkpt)

    # determine lower and upper bound for bisection
    Elw = minimum(ebands[1,:]) # minimum for all kspin
    Eup = maximum(ebands[Nstates,:]) # maximum for all kspin

    Elw = Elw - 2*kT
    Eup = Eup + 2*kT

    verbose && println("Elw = ", Elw)
    verbose && println("Eup = ", Eup)

    sumlw = sum_Focc( smear_func, ebands, Elw, kT, Nkpt, wk )
    sumup = sum_Focc( smear_func, ebands, Eup, kT, Nkpt, wk )

    SMALL = 1e-10

    if ( (sumup - Nelectrons) < -eps() ) ||
       ( (sumlw - Nelectrons) >  eps() )
        @printf("sumup = %18.10f\n", sumup)
        @printf("sumlw = %18.10f\n", sumlw)
        @printf("Nelectrons = %18.10f\n", Nelectrons)
        error("Bounds for E_fermi is not found")
    end

    #
    # Start bisection
    #
    Ef = 0.5*(Eup + Elw)
    Ef_old = Ef
    for iter in 1:NiterMax
        sum_mid = sum_Focc( smear_func, ebands, Ef, kT, Nkpt, wk )
        if abs(sum_mid-Nelectrons) < SMALL
            verbose && println("Fermi converged: diff_Ef = $(abs(sum_mid-Nelectrons)), E_fermi=$(Ef)")
            return Ef
        elseif (sum_mid-Nelectrons) < -SMALL
            Elw = Ef
        else
            Eup = Ef
        end
        Ef = 0.5*(Eup + Elw)
        diff_Ef = abs(Ef-Ef_old)
        if verbose
            @printf("find_E_fermi: %3d %18.10f %18.10f %18.10e\n", iter, Ef, sum_mid, diff_Ef)
        end
        if diff_Ef < SMALL
            verbose && println("Fermi converged: diff_Ef = $(diff_Ef), E_fermi=$(Ef)")
            return Ef
        end
        Ef_old = Ef
    end

    @warn "WARNING: Ef is not found after $(NiterMax) iterations"
    return Ef
    
end