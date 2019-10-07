"""
Similar to `calc_primary_search_dirs`, but using the gradients explicitly.
FIXME: If this is ever used, no need to call `calc_grad_Haux` again.
Simply use the one that has been calculated previously.
"""
function calc_primary_search_dirs_v1!(
    Ham::Hamiltonian,
    evars::ElectronicVars,
    Δ_evars::ElectronicVars;
    kT=1e-3,
    κ=0.5
)

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    Δ_ψ = Δ_evars.ψ
    Δ_η = Δ_evars.η
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin - 1)*Nkpt
        Δ_ψ[i], Δ_η[i] = calc_grad_Haux(Ham, evars.ψ[i], kT)
        Δ_η[i] = -κ*Δ_η[i]  # XXX reverse the sign of Δ_η ?
    end

    return

end
