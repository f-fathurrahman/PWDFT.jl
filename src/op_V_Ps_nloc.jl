function op_V_Ps_nloc( Ham::Hamiltonian, psiks::BlochWavefunc )
    Nstates = size(psiks[1])[2] # Nstates should be similar for all Bloch states
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    out = zeros_BlochWavefunc(Ham)
    
    for ispin = 1:Nspin
    for ik=1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        out[ikspin] = op_V_Ps_nloc( Ham, psiks[ikspin] )
    end
    end
    return out
end

function op_V_Ps_nloc( Ham::Hamiltonian, psi::Array{ComplexF64,2} )

    ik = Ham.ik

    # Take `Nstates` to be the size of psi and not from `Ham.electrons.Nstates`.
    Nstates = size(psi)[2]

    # first dimension of psi should be Ngw[ik]

    atoms = Ham.atoms
    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    Pspots = Ham.pspots
    prj2beta = Ham.pspotNL.prj2beta
    betaNL = Ham.pspotNL.betaNL

    betaNL_psi = calc_betaNL_psi( ik, Ham.pspotNL.betaNL, psi )
    
    Ngw = Ham.pw.gvecw.Ngw
    Vpsi = zeros( ComplexF64, Ngw[ik], Nstates )

    for ist = 1:Nstates
        for ia = 1:Natoms
            isp = atm2species[ia]
            psp = Pspots[isp]
            for l = 0:psp.lmax
            for m = -l:l
                for iprj = 1:psp.Nproj_l[l+1]
                for jprj = 1:psp.Nproj_l[l+1]
                    ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                    jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                    hij = psp.h[l+1,iprj,jprj]
                    for ig = 1:Ngw[ik]
                        Vpsi[ig,ist] = Vpsi[ig,ist] + hij*betaNL[ig,ibeta,ik]*betaNL_psi[ist,jbeta]
                    end
                end # iprj
                end # jprj
            end # m
            end # l
        end
    end
    return Vpsi
end

