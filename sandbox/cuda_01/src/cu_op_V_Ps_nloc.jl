import PWDFT: op_V_Ps_nloc

function op_V_Ps_nloc( Ham::CuHamiltonian, psiks::CuBlochWavefunc )
    
    Nstates = size(psiks[1])[2] # Nstates should be similar for all Bloch states
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    out = zeros_CuBlochWavefunc(Ham)
    
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        out[ikspin] = op_V_Ps_nloc( Ham, psiks[ikspin] )
    end
    return out
end


function kernel_op_V_Ps_nloc!( ist, ibeta, jbeta, hij, betaNL_ik, betaNL_psi, Vpsi )
    
    igk = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    Ngw_k = size(Vpsi,1)
    
    if igk <= Ngw_k
        Vpsi[igk,ist] = Vpsi[igk,ist] + hij*betaNL_ik[igk,ibeta]*betaNL_psi[ist,jbeta]
    end

    return
end


function kernel_op_V_Ps_nloc_1state!( ibeta, jbeta, hij, betaNL_ik, betaNL_psi, Vpsi )
    
    igk = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    Ngw_k = size(Vpsi,1)
    
    if igk <= Ngw_k
        Vpsi[igk] = Vpsi[igk] + hij*betaNL_ik[igk,ibeta]*betaNL_psi[jbeta]
    end

    return
end


function op_V_Ps_nloc( Ham::CuHamiltonian, psi::CuArray{ComplexF64,2} )

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
    Vpsi = CuArrays.zeros( ComplexF64, Ngw[ik], Nstates )

    Nthreads = 256
    Nblocks = ceil(Int64, Ngw[ik]/Nthreads)

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
                    #Vpsi[:,ist] = Vpsi[:,ist] + hij*betaNL[ik][:,ibeta]*betaNL_psi[ist,jbeta]
                    @cuda threads=Nthreads blocks=Nblocks kernel_op_V_Ps_nloc!( ist, ibeta, jbeta, hij, betaNL[ik], betaNL_psi, Vpsi )
                end # iprj
                end # jprj
            end # m
            end # l
        end
    end
    return Vpsi
end

# FIXME: remove redundant code
function op_V_Ps_nloc( Ham::CuHamiltonian, psi::CuArray{ComplexF64,1} )

    ik = Ham.ik

    # first dimension of psi should be Ngw[ik]

    atoms = Ham.atoms
    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    Pspots = Ham.pspots
    prj2beta = Ham.pspotNL.prj2beta
    betaNL = Ham.pspotNL.betaNL

    betaNL_psi = calc_betaNL_psi( ik, Ham.pspotNL.betaNL, psi )
    
    Ngw = Ham.pw.gvecw.Ngw
    Vpsi = CuArrays.zeros( ComplexF64, Ngw[ik] )

    Nthreads = 256
    Nblocks = ceil(Int64, Ngw[ik]/Nthreads)

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
                #Vpsi[:] = Vpsi[:] + hij*betaNL[ik][:,ibeta]*betaNL_psi[jbeta]
                @cuda threads=Nthreads blocks=Nblocks kernel_op_V_Ps_nloc( ist, ibeta, jbeta, hij, betaNL[ik], betaNL_psi, Vpsi )
            end # iprj
            end # jprj
        end # m
        end # l
    end

    return Vpsi
end
