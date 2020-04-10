import PWDFT: calc_energies, calc_E_kin, calc_E_local, calc_E_Ps_nloc

function kernel_calc_E_Ps_nloc!(
    Natoms::Int64,
    atm2species,
    psp_lmax,
    psp_Nproj_l,
    psp_h,
    ik::Int64,
    ikspin::Int64,
    wk,
    Focc,
    prj2beta,
    betaNL_psi,
    e_ps_nloc_k
)
    
    ist = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x

    Nstates = length(e_ps_nloc_k)

    if ist <= Nstates

        enl1 = 0.0
        for ia = 1:Natoms
            
            isp = atm2species[ia]
            
            for l = 0:psp_lmax[isp]
            for m = -l:l
            for iprj = 1:psp_Nproj_l[l+1,isp]
            for jprj = 1:psp_Nproj_l[l+1,isp]
                
                ibeta = prj2beta[iprj, ia, l+1, m+psp_lmax[isp]+1]
                
                jbeta = prj2beta[jprj, ia, l+1, m+psp_lmax[isp]+1]
                
                hij = psp_h[l+1, iprj, jprj, isp]

                enl1 = enl1 + hij*real( conj(betaNL_psi[ist,ibeta]) * betaNL_psi[ist,jbeta] )

            end
            end
            end # m
            end # l
        end
        
        e_ps_nloc_k[ist] = wk[ik]*Focc[ist,ikspin]*enl1
    
    end

    return
end


function calc_E_Ps_nloc( Ham::CuHamiltonian, psiks::CuBlochWavefunc )

    Nstates = Ham.electrons.Nstates
    Natoms = Ham.atoms.Natoms
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    NbetaNL = Ham.pspotNL.NbetaNL
    Nspin = Ham.electrons.Nspin

    atm2species = Ham.pspotNL.atm2species
    prj2beta = Ham.pspotNL.prj2beta_gpu
    Focc = Ham.pspotNL.Focc
    wk = Ham.pspotNL.wk
    psp_lmax = Ham.pspotNL.psp_lmax
    psp_Nproj_l = Ham.pspotNL.psp_Nproj_l
    psp_h = Ham.pspotNL.psp_h

    E_Ps_nloc = 0.0

    betaNL_psi = CuArrays.zeros(ComplexF64,Nstates,NbetaNL)

    e_ps_nloc_k = CuArrays.zeros(Float64, Nstates)

    Nthreads = min( 256, Nstates )
    Nblocks = ceil(Int64, Nstates/Nthreads)

    for ispin in 1:Nspin, ik in 1:Nkpt
        
        ikspin = ik + (ispin - 1)*Nkpt
        
        betaNL_psi[:] = calc_betaNL_psi( ik, Ham.pspotNL.betaNL, psiks[ikspin] )

        @cuda threads=Nthreads blocks=Nblocks kernel_calc_E_Ps_nloc!(
                Natoms, atm2species,
                psp_lmax, psp_Nproj_l, psp_h,
                ik, ikspin,
                wk, Focc,
                prj2beta, betaNL_psi,
                e_ps_nloc_k )

        E_Ps_nloc = E_Ps_nloc + sum( e_ps_nloc_k )
    end

    return E_Ps_nloc

end

# Calculate Ps loc, Hartree, and XC components of total energy
function calc_E_local( Ham::CuHamiltonian )

    Npoints = prod(Ham.pw.Ns)
    dVol = Ham.pw.CellVolume/Npoints
    Nspin = Ham.electrons.Nspin
    potentials = Ham.potentials

    Rhoe_tot = CuArrays.zeros(Float64, Npoints)
    for ispin = 1:Nspin
        Rhoe_tot[:] = Rhoe_tot[:] + Ham.rhoe[:,ispin]
    end

    E_Hartree = 0.5*dot( potentials.Hartree, Rhoe_tot ) * dVol
    E_Ps_loc = dot( potentials.Ps_loc, Rhoe_tot ) * dVol

    if Ham.xcfunc == "PBE"
        epsxc = calc_epsxc_PBE( Ham.xc_calc, Ham.pw, Ham.rhoe )
    else
        epsxc = calc_epsxc_VWN( Ham.xc_calc, Ham.rhoe )
    end
    E_xc = dot( epsxc, Rhoe_tot ) * dVol

    return E_Ps_loc, E_Hartree, E_xc
end




function kernel_calc_E_kin!( ist, idx_gw2g_ik, G, k1, k2, k3, psi, psiKpsi )
    
    igk = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x

    Ngw_ik = length( idx_gw2g_ik )

    if igk <= Ngw_ik
        
        ig = idx_gw2g_ik[igk]
        
        Gw2 = ( G[1,ig] + k1 )^2 + ( G[2,ig] + k2 )^2 + ( G[3,ig] + k3 )^2
        
        psiKpsi[igk] = CUDAnative.abs( psi[igk,ist] )^2*Gw2
    
    end
    
    return
end


function calc_E_kin( Ham::CuHamiltonian, psiks::CuBlochWavefunc )

    Focc = Ham.electrons.Focc
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nstates = Ham.electrons.Nstates
    wk = Ham.pw.gvecw.kpoints.wk
    Nspin = Ham.electrons.Nspin

    Ngw = Ham.pw.gvecw.Ngw
    idx_gw2g = Ham.pw.gvecw.idx_gw2g
    G = Ham.pw.gvec.G
    k = Ham.pw.gvecw.kpoints.k

    Ngwx = Ham.pw.gvecw.Ngwx
    psiKpsi = CuArrays.zeros(Float64, Ngwx)

    Nthreads = 256

    E_kin = 0.0
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        psi = psiks[ikspin]
        k1 = k[1,ik]
        k2 = k[2,ik]
        k3 = k[3,ik]
        for ist = 1:Nstates
            
            psiKpsi[:] .= 0.0 # need this ?
            
            Nblocks = ceil(Int64, Ngw[ik]/Nthreads)

            @cuda threads=Nthreads blocks=Nblocks kernel_calc_E_kin!( ist, idx_gw2g[ik], G, k1, k2, k3, psi, psiKpsi )

            E_kin = E_kin + wk[ik]*Focc[ist,ikspin]*sum( psiKpsi[1:Ngw[ik]] )

        end
    end
    end
    return 0.5*E_kin

end

#
# psi is assumed to be already orthonormalized elsewhere
# `potentials` and `Rhoe` are not updated
# Ham is assumed to be already updated at input psi
#
# Ham.energies.NN should be calculated outside this function
function calc_energies( Ham::CuHamiltonian, psiks::CuBlochWavefunc )
    
    E_kin = calc_E_kin( Ham, psiks )

    E_Ps_loc, E_Hartree, E_xc = calc_E_local( Ham )

    if Ham.pspotNL.NbetaNL > 0
        E_Ps_nloc = calc_E_Ps_nloc( Ham, psiks )
    else
        E_Ps_nloc = 0.0
    end

    energies = Energies()
    energies.Kinetic = E_kin
    energies.Ps_loc  = E_Ps_loc
    energies.Ps_nloc = E_Ps_nloc
    energies.Hartree = E_Hartree
    energies.XC      = E_xc
    energies.NN      = Ham.energies.NN

    return energies
end

