function calc_forces_scf_corr( Ham::Hamiltonian )
    F_scf_corr = zeros(Float64, 3, Ham.atoms.Natoms)
    calc_forces_scf_corr!(Ham, F_scf_corr)
    return F_scf_corr
end

# Set to zeros in case of PsPot_GTH
function calc_forces_scf_corr!( Ham::Hamiltonian{PsPot_GTH}, F_scf_corr )
    fill!(F_scf_corr, 0.0)
    return
end

function calc_forces_scf_corr!( Ham::Hamiltonian{PsPot_UPF}, F_scf_corr )
    calc_forces_scf_corr!(
        Ham.atoms, Ham.pw, Ham.pspots, Ham.potentials, F_scf_corr
    )
    return
end


# XXX: This currently only defined for PsPot_UPF because PsPot_GTH does not
# XXX: provide rhoatom information.
# XXX: This term should be small for converged SCF.
function calc_forces_scf_corr!(
    atoms::Atoms,
    pw::PWGrid,
    pspots::Vector{PsPot_UPF},
    potentials::Potentials,
    F_scf_corr::Matrix{Float64}
)

    # Calculate difference in potential
    VtotOld = potentials.TotalOld # actually Hartree and XC
    #!!! Check whether VtotOld is properly saved in update_from_rhoe
    
    eps8 = 1e-8
    G2_shells = pw.gvec.G2_shells
    idx_g2r = pw.gvec.idx_g2r
    G = pw.gvec.G
    idx_g2shells = pw.gvec.idx_g2shells
    Ngl = length(G2_shells)
    Ng = pw.gvec.Ng

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species
    X = atoms.positions

    Nspin = size(VtotOld,2)
    Npoints = size(VtotOld,1)

    # Calculate difference between new and old potential
    ctmp = zeros(ComplexF64, Npoints)
    for ispin in 1:Nspin, ip in 1:Npoints
        ctmp[ip] += ( potentials.Hartree[ip] + potentials.XC[ip,ispin] - VtotOld[ip,ispin] )
    end
    # FIXME: probably need to calculate this only just after convergence is reached

    R_to_G!(pw, ctmp)
    @views ctmp[:] /= Npoints # rescale

    rhocgnt = zeros(Float64, Ngl)
    fill!(F_scf_corr, 0.0)

    for isp in 1:Nspecies
        psp = pspots[isp]
        aux = zeros(Float64, psp.Nr)
        # G != 0 terms
        for igl in 2:Ngl
            Gm = sqrt(G2_shells[igl])
            for ir in 1:psp.Nr
                if psp.r[ir] < eps8
                   aux[ir] = psp.rhoatom[ir]
                else
                   aux[ir] = psp.rhoatom[ir]*sin(Gm*psp.r[ir])/(psp.r[ir]*Gm)
                end
            end
            rhocgnt[igl] = integ_simpson( psp.Nr, aux, psp.rab )
        end
        #
        # sum over atoms
        for ia in 1:Natoms
            if isp != atm2species[ia]
                continue
            end
            #
            for ig in 2:Ng
                igl = idx_g2shells[ig]
                ip = idx_g2r[ig]
                GX = G[1,ig]*X[1,ia] + G[2,ig]*X[2,ia] + G[3,ig]*X[3,ia]
                Sf = sin(GX) + im*cos(GX)
                Sfpsi = real(Sf*conj(ctmp[ip]))
                @views F_scf_corr[:,ia] .+= rhocgnt[igl] * G[:,ig] * Sfpsi
            end
        end

    end

    return

end