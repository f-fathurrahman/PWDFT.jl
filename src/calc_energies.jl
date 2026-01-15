"""
    calc_E_Ps_nloc(Ham, psiks)

Compute and return non-local pseudopotential energy for a given Hamiltonian `Ham` (with properly updated
electron density) and BlochWavefunc `psiks`.
"""
function calc_E_Ps_nloc(
    Ham::Hamiltonian{PsPot_GTH},
    psiks::BlochWavefunc
)

    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    Natoms = Ham.atoms.Natoms
    atm2species = Ham.atoms.atm2species
    Pspots = Ham.pspots
    prj2beta = Ham.pspotNL.prj2beta
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    wk = Ham.pw.gvecw.kpoints.wk
    NbetaNL = Ham.pspotNL.NbetaNL
    Nspin = Ham.electrons.Nspin_channel

    # calculate E_NL
    E_Ps_nloc = 0.0

    betaNL_psi = zeros(ComplexF64,Nstates,NbetaNL)

    # The loops can be simplified, especially to conform with the usual
    # multiple projectors per angular momentum channel.
    # This is the form that I found closely resemble
    # the equations found in the papers/books

    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        psi = psiks[ikspin]
        @views betaNL_psi[:,:] = calc_betaNL_psi( ik, Ham.pspotNL.betaNL, psi )
        for ist in 1:Nstates
            enl1 = 0.0
            for ia in 1:Natoms
                isp = atm2species[ia]
                psp = Pspots[isp]
                for l in 0:psp.lmax, m in -l:l
                    for iprj in 1:psp.Nproj_l[l+1], jprj in 1:psp.Nproj_l[l+1]
                        ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                        jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                        hij = psp.h[l+1,iprj,jprj]
                        enl1 += hij*real(conj(betaNL_psi[ist,ibeta])*betaNL_psi[ist,jbeta])
                    end
                end
            end
            # accumulate
            E_Ps_nloc += wk[ik]*Focc[ist,ikspin]*enl1
        end
    end

    return E_Ps_nloc

end


function calc_E_Ps_nloc(
    Ham::Hamiltonian{PsPot_UPF},
    psiks::BlochWavefunc
)
    if Ham.options.noncollinear
        return calc_E_Ps_nloc_noncollinear(Ham, psiks)
    end

    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    Nspin = Ham.electrons.Nspin_channel

    Natoms = Ham.atoms.Natoms
    Nspecies = Ham.atoms.Nspecies
    atm2species = Ham.atoms.atm2species

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    wk = Ham.pw.gvecw.kpoints.wk
    Ngw = Ham.pw.gvecw.Ngw
    Ngwx = maximum(Ngw)

    NbetaNL = Ham.pspotNL.NbetaNL
    Dvan = Ham.pspotNL.Dvan # XXX Note that we are using Dvan instead of Deeq here
    nh = Ham.pspotNL.nh
    indv_ijkb0 = Ham.pspotNL.indv_ijkb0

    betaNL_psi = zeros(ComplexF64, NbetaNL, Nstates)
    ps = zeros(ComplexF64, NbetaNL, Nstates)
    Vpsi = zeros(ComplexF64, Ngwx, Nstates)

    E_Ps_nloc = 0.0

    for ispin in 1:Nspin, ik in 1:Nkpt
        
        ikspin = ik + (ispin - 1)*Nkpt
        psi = psiks[ikspin]
        
        #
        # Here is the operation of op_V_Ps_nloc
        #
        betaNL_k = Ham.pspotNL.betaNL[ik]

        mul!(betaNL_psi, betaNL_k', psi)

        for isp in 1:Nspecies
            if nh[isp] == 0
                continue
            end
            for ia in 1:Natoms
                if atm2species[ia] != isp
                    continue
                end
                idx1 = (indv_ijkb0[ia]+1):(indv_ijkb0[ia]+nh[isp])
                @views ps[idx1,:] = Dvan[1:nh[isp],1:nh[isp],isp] * betaNL_psi[idx1,:]
            end
        end

        # Accumulate
        @views Vpsi[1:Ngw[ik],:] = betaNL_k * ps
        for ist in 1:Nstates
            @views ss = real(dot(psi[:,ist], Vpsi[1:Ngw[ik],ist]))
            E_Ps_nloc += wk[ik]*Focc[ist,ikspin]*ss
        end

    end

    return E_Ps_nloc

end


function calc_E_Ps_nloc_noncollinear(
    Ham::Hamiltonian{PsPot_UPF},
    psiks::BlochWavefunc
)

    Npol = 2 # HARDCODED
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc

    Natoms = Ham.atoms.Natoms
    atm2species = Ham.atoms.atm2species

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    wk = Ham.pw.gvecw.kpoints.wk
    Ngw = Ham.pw.gvecw.Ngw
    Ngwx = maximum(Ngw)

    # XXX Note that we are using Dvan instead of Deeq here
    if Ham.options.lspinorb
        Dvan = Ham.pspotNL.Dvan_so
    else
        Dvan = Ham.pspotNL.Dvan
    end

    NbetaNL = Ham.pspotNL.NbetaNL
    nh = Ham.pspotNL.nh
    indv_ijkb0 = Ham.pspotNL.indv_ijkb0

    ps = zeros(ComplexF64, NbetaNL, Npol, Nstates)
    Vpsi = zeros(ComplexF64, Ngwx, Npol, Nstates)
    betaNL_psi = zeros(ComplexF64, NbetaNL, Npol, Nstates)

    E_Ps_nloc = 0.0

    for ik in 1:Nkpt
        
        psir = reshape(psiks[ik], Ngw[ik], Npol, Nstates)
        
        #
        # Here is the operation of op_V_Ps_nloc
        #
        betaNL_k = Ham.pspotNL.betaNL[ik]        
        betaNL_psi[:,1,:] = betaNL_k' * psir[:,1,:]
        betaNL_psi[:,2,:] = betaNL_k' * psir[:,2,:]

        for ia in 1:Natoms
            isp = atm2species[ia]
            if nh[isp] == 0
                continue
            end
            for ist in 1:Nstates
                for jh in 1:nh[isp]
                    jkb = indv_ijkb0[ia] + jh
                    for ih in 1:nh[isp]
                        ikb = indv_ijkb0[ia] + ih
                        ps[ikb,1,ist] += Dvan[ih,jh,1,isp]*betaNL_psi[jkb,1,ist] + 
                                         Dvan[ih,jh,2,isp]*betaNL_psi[jkb,2,ist] 
                        ps[ikb,2,ist] += Dvan[ih,jh,3,isp]*betaNL_psi[jkb,1,ist] +
                                         Dvan[ih,jh,4,isp]*betaNL_psi[jkb,2,ist]
                    end
                end
            end
        end # loop over Natoms

        # Accumulate
        @views Vpsi[1:Ngw[ik],1,:] = betaNL_k * ps[:,1,:]
        @views Vpsi[1:Ngw[ik],2,:] = betaNL_k * ps[:,2,:]
        for ist in 1:Nstates
            # XXX: something missing here?
            @views ss = real(dot(psir[:,1:Npol,ist], Vpsi[1:Ngw[ik],1:Npol,ist]))
            E_Ps_nloc += wk[ik]*Focc[ist,ik]*ss
        end

    end

    return E_Ps_nloc

end





"""
    calc_E_local(Ham)

Compute and return local energy terms (local pseudopotential, Hartree, and XC) for a given
Hamiltonian `Ham` with electron density stored in `Ham.rhoe`.
"""
function calc_E_local( Ham::Hamiltonian, psiks::BlochWavefunc )

    Npoints = prod(Ham.pw.Ns)
    dVol = Ham.pw.CellVolume/Npoints
    Nspin = Ham.electrons.Nspin_comp
    potentials = Ham.potentials

    Rhoe = Ham.rhoe # alias
    Rhoe_tot = zeros(Float64, Npoints)
    if Nspin in [1,2]
        # XXX: use reduce?
        for ispin in 1:Nspin, ip in 1:Npoints
            Rhoe_tot[ip] += Rhoe[ip,ispin]
        end
    else
        # Nspin == 4
        Rhoe_tot[:] = Rhoe[:,1]
    end

    E_Hartree = 0.5*dot( potentials.Hartree, Rhoe_tot ) * dVol
    E_Ps_loc = dot( potentials.Ps_loc, Rhoe_tot ) * dVol
    
    # XC energy is handled by different function
    E_xc = calc_E_xc(Ham, psiks)

    return E_Ps_loc, E_Hartree, E_xc
end

function calc_E_xc(Ham, psiks)

    Rhoe = Ham.rhoe # alias    
    Npoints = prod(Ham.pw.Ns)
    dVol = Ham.pw.CellVolume/Npoints
    Nspin = Ham.electrons.Nspin_comp

    if !isnothing(Ham.rhoe_core)
        if Nspin == 2
            Rhoe[:,1] .+= Ham.rhoe_core*0.5
            Rhoe[:,2] .+= Ham.rhoe_core*0.5
        else # for both Nspin == 1 and Nspin == 4
            Rhoe[:,1] .+= Ham.rhoe_core
        end
    end

    domag = Ham.electrons.domag
    epsxc = zeros(Float64, Npoints)
    if Ham.xcfunc == "VWN"
        # VWN, handle 
        if Nspin <= 2
            calc_epsxc_VWN!( Ham.xc_calc, Rhoe, epsxc )
        elseif Nspin == 4
            if domag
                calc_epsxc_VWN_noncollinear!( Ham.xc_calc, Rhoe, epsxc )
            else
                # XXX Special case for noncollinear, not using magnetism
                @views calc_epsxc_VWN!( Ham.xc_calc, Rhoe[:,1], epsxc )
            end
        else
            @error("Wrong Nspin=$Nspin")
        end 
    elseif Ham.xcfunc == "PBE"
        # check if this is working
        calc_epsxc_PBE!( Ham.xc_calc, Ham.pw, Rhoe, epsxc )
    else
        # FIXME: no spinpol and core-correction yet
        @assert isnothing(Ham.rhoe_core)
        @assert Nspin == 1
        calc_epsxc_SCAN!( Ham, psiks, Rhoe, epsxc )
    end
    if Nspin == 4
        E_xc = sum(epsxc .* Rhoe[:,1])*dVol
    else
        E_xc = sum(epsxc .* Rhoe)*dVol
    end

    if !isnothing(Ham.rhoe_core)
        # Recover
        if Nspin == 2
            Rhoe[:,1] .-= Ham.rhoe_core*0.5
            Rhoe[:,2] .-= Ham.rhoe_core*0.5
        else
            Rhoe[:,1] .-= Ham.rhoe_core
        end
    end

    return E_xc
end



"""
    calc_E_kin(Ham, psiks)

Compute and return kinetic energy term for a given Hamiltonian `Ham` and BlochWavefunc `psiks`.
"""
function calc_E_kin( Ham::Hamiltonian, psiks::BlochWavefunc )

    if Ham.options.noncollinear
        return calc_E_kin_noncollinear(Ham, psiks)
    end

    Focc = Ham.electrons.Focc
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nstates = Ham.electrons.Nstates
    wk = Ham.pw.gvecw.kpoints.wk
    Nspin = Ham.electrons.Nspin_channel

    Ngw = Ham.pw.gvecw.Ngw
    idx_gw2g = Ham.pw.gvecw.idx_gw2g
    G = Ham.pw.gvec.G
    k = Ham.pw.gvecw.kpoints.k

    E_kin = 0.0
    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        psi = psiks[ikspin]
        for ist in 1:Nstates
            psiKpsi = 0.0
            for igk in 1:Ngw[ik]
                ig = idx_gw2g[ik][igk]
                Gw2 = (G[1,ig] + k[1,ik])^2 + (G[2,ig] + k[2,ik])^2 + (G[3,ig] + k[3,ik])^2
                psiKpsi += abs(psi[igk,ist])^2*Gw2
            end
            E_kin += wk[ik]*Focc[ist,ikspin]*psiKpsi
        end
    end
    return 0.5*E_kin

end

function calc_E_kin_noncollinear( Ham::Hamiltonian, psiks::BlochWavefunc )

    Npol = 2 # HARDCODED

    Focc = Ham.electrons.Focc
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nstates = Ham.electrons.Nstates
    wk = Ham.pw.gvecw.kpoints.wk

    Ngw = Ham.pw.gvecw.Ngw
    idx_gw2g = Ham.pw.gvecw.idx_gw2g
    G = Ham.pw.gvec.G
    k = Ham.pw.gvecw.kpoints.k

    E_kin = 0.0
    for ik in 1:Nkpt
        psir = reshape(psiks[ik], (Ngw[ik], 2, Nstates))
        for ist in 1:Nstates
            psiKpsi = 0.0
            for ipol in 1:Npol, igk in 1:Ngw[ik]
                ig = idx_gw2g[ik][igk]
                Gw2 = (G[1,ig] + k[1,ik])^2 + (G[2,ig] + k[2,ik])^2 + (G[3,ig] + k[3,ik])^2
                psiKpsi += abs(psir[igk,ipol,ist])^2*Gw2
            end
            E_kin += wk[ik]*Focc[ist,ik]*psiKpsi
        end
    end
    return 0.5*E_kin

end



"""
    calc_energies(Ham, psiks)

Compute and return total energy components of type `Energies` for a given Hamiltonian `Ham`
and BlochWavefunc `psiks`.

Each `psiks` is assumed to be already orthonormalized elsewhere.

`Ham.potentials` and `Ham.rhoe` are not updated.

`Ham.energies.NN` should be calculated outside this function if needed.
"""
function calc_energies( Ham::Hamiltonian, psiks::BlochWavefunc )
    calc_energies!(Ham, psiks) # call the in-place version
    return Ham.energies
end

# This is the in-place version
function calc_energies!( Ham::Hamiltonian, psiks::BlochWavefunc )
    
    E_kin = calc_E_kin( Ham, psiks )

    E_Ps_loc, E_Hartree, E_xc = calc_E_local( Ham, psiks )

    if Ham.pspotNL.NbetaNL > 0
        E_Ps_nloc = calc_E_Ps_nloc( Ham, psiks )
    else
        E_Ps_nloc = 0.0
    end

    Ham.energies.Kinetic = E_kin
    Ham.energies.Ps_loc  = E_Ps_loc
    Ham.energies.Ps_nloc = E_Ps_nloc
    Ham.energies.Hartree = E_Hartree
    Ham.energies.XC      = E_xc
    Ham.energies.NN      = Ham.energies.NN

    ok_paw = any(Ham.pspotNL.are_paw)
    if ok_paw
        Ham.energies.EHxc_paw = Ham.pspotNL.paw.EHxc_paw
    end

    return
end
