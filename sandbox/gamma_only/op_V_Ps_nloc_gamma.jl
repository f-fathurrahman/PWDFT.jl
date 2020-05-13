import PWDFT: op_V_Ps_nloc

function op_V_Ps_nloc( Ham::HamiltonianGamma, psis::BlochWavefuncGamma )
    
    Nstates = size(psis.data[1],2) # Nstates should be similar for all Bloch states
    Nspin = Ham.electrons.Nspin
    out = zeros_BlochWavefuncGamma(Ham)
    
    for ispin in 1:Nspin
        out.data[ispin] = op_V_Ps_nloc( Ham, psis.data[ispin] )
    end

    return out

end


function op_V_Ps_nloc( Ham::HamiltonianGamma, psi::Array{ComplexF64,2} )

    Nstates = size(psi,2)

    atoms = Ham.atoms
    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    
    pspots = Ham.pspots
    
    prj2beta = Ham.pspotNL.prj2beta
    betaNL = Ham.pspotNL.betaNL

    betaNL_psi = calc_betaNL_psi( Ham.pspotNL, psi )
    
    Ngw = Ham.pw.gvecw.Ngw
    Vpsi = zeros( ComplexF64, Ngw, Nstates )

    for ist = 1:Nstates
        for ia = 1:Natoms
            isp = atm2species[ia]
            psp = pspots[isp]
            for l = 0:psp.lmax
            for m = -l:l
                for iprj = 1:psp.Nproj_l[l+1]
                for jprj = 1:psp.Nproj_l[l+1]
                    ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                    jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                    hij = psp.h[l+1,iprj,jprj]
                    for igw in 1:Ngw
                        Vpsi[igw,ist] = Vpsi[igw,ist] + hij*betaNL[igw,ibeta]*betaNL_psi[ist,jbeta]
                    end
                end # iprj
                end # jprj
            end # m
            end # l
        end
    end

    #println("betaNL_psi = ", )
    #display(betaNL_psi); println()

    #println(Vpsi[1,1:Nstates])
    
    return Vpsi
end


# FIXME: remove redundant code
function op_V_Ps_nloc( Ham::HamiltonianGamma, psi::Array{ComplexF64,1} )

    atoms = Ham.atoms
    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    Pspots = Ham.pspots
    prj2beta = Ham.pspotNL.prj2beta
    betaNL = Ham.pspotNL.betaNL

    betaNL_psi = calc_betaNL_psi( Ham.pspotNL, psi )
    
    Ngw = Ham.pw.gvecw.Ngw
    Vpsi = zeros( ComplexF64, Ngw )

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
                for igw in 1:Ngw[ik]
                    Vpsi[igw] = Vpsi[igw] + hij*betaNL[igw,ibeta]*betaNL_psi[jbeta]
                end
            end # iprj
            end # jprj
        end # m
        end # l
    end

    return Vpsi
end
