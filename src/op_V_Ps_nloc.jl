function op_V_Ps_nloc( Ham::PWHamiltonian, psi::Array{Complex128,2} )
    
    Nstates = size(psi)[2]
    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx

    atoms = Ham.atoms
    atm2species = atoms.atm2species

    Pspots = Ham.pspots

    # calculate op_V_NL | psi >
    Vpsi = zeros( Complex128, Ngwx, Nstates )

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
                    Vpsi[:,ist] = Vpsi[:,ist] + hij * betaNL[:,ibeta] * betaNL_psi[ia,ist,jbeta]
                end # iprj
                end # jprj
            end # m
            end # l
        end
    end
end

