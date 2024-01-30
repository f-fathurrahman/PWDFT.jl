function calc_forces_scf_corr!( Ham::Hamiltonian{PsPot_UPF}, F_scf_corr )
    calc_forces_scf_corr!(
        Ham.atoms, Ham.pw, Ham.pspots, Ham.potentials, F_scf_corr
    )
    return
end



function calc_forces_scf_corr!(
    atoms::Atoms,
    pw::PWGrid,
    pspots::Vector{PsPot_UPF},
    potentials::Potentials,
    F_scf_corr::Matrix{Float64}
)

    println("\nBegin calculating F_scf_corr")

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
    
    println("Ng = ", Ng)
    println("Npoints = ", Npoints)
    println("pw.Ns = ", pw.Ns)

    println("sum abs VHartree = ", sum(abs.(potentials.Hartree)))
    println("sum XC = ", sum(potentials.XC))

    # Calculate difference between new and old potential
    ctmp = zeros(ComplexF64, Npoints)
    for ispin in 1:Nspin, ip in 1:Npoints
        ctmp[ip] += ( potentials.Hartree[ip] + potentials.XC[ip,ispin] - VtotOld[ip,ispin] )
    end

    # FIXME: probably need to calculate this only just after convergence is reached
    println("sum ctmp before fwfft = ", sum(ctmp))

    R_to_G!(pw, ctmp)
    #ff = reshape(ctmp, pw.Ns)
    ## Use FFT again
    #planfw = plan_fft!( zeros(ComplexF64,pw.Ns) ) # using default plan
    #planfw*ff
    @views ctmp[:] /= Npoints # rescale

    println("sum ctmp after forward FFT = ", sum(ctmp))

    rhocgnt = zeros(Float64, Ngl)
    fill!(F_scf_corr, 0.0)

    for isp in 1:Nspecies
        psp = pspots[isp]
        aux = zeros(Float64, psp.Nr)
        # G != 0 terms
        for igl in 2:Ngl
            gx = sqrt(G2_shells[igl])
            for ir in 1:psp.Nr
                if psp.r[ir] < eps8
                   aux[ir] = psp.rhoatom[ir]
                else
                   aux[ir] = psp.rhoatom[ir]*sin(gx*psp.r[ir])/(psp.r[ir]*gx)
                end
            end
            rhocgnt[igl] = PWDFT.integ_simpson( psp.Nr, aux, psp.rab )
        end
        #
        println("sum rhocgnt = ", sum(rhocgnt))
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