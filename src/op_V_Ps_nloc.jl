function op_V_Ps_nloc( Ham::PWHamiltonian, psi::Array{Complex128,2} )

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
    Vpsi = zeros( Complex128, Ngw[ik], Nstates )

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
                    Vpsi[:,ist] = Vpsi[:,ist] + hij * betaNL[1:Ngw[ik],ibeta,ik] * betaNL_psi[ist,jbeta]
                end # iprj
                end # jprj
            end # m
            end # l
        end
    end
    return Vpsi
end

