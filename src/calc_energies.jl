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
    Nspin = Ham.electrons.Nspin

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

    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    Nspin = Ham.electrons.Nspin

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




"""
    calc_E_local(Ham)

Compute and return local energy terms (local pseudopotential, Hartree, and XC) for a given
Hamiltonian `Ham` with electron density stored in `Ham.rhoe`.
"""
function calc_E_local( Ham::Hamiltonian, psiks::BlochWavefunc )

    Npoints = prod(Ham.pw.Ns)
    dVol = Ham.pw.CellVolume/Npoints
    Nspin = Ham.electrons.Nspin
    potentials = Ham.potentials

    Rhoe = Ham.rhoe # alias
    Rhoe_tot = zeros(Float64, Npoints)
    for ispin in 1:Nspin, ip in 1:Npoints
        Rhoe_tot[ip] += Rhoe[ip,ispin]
    end
    # XXX: use reduce?

    E_Hartree = 0.5*dot( potentials.Hartree, Rhoe_tot ) * dVol
    E_Ps_loc = dot( potentials.Ps_loc, Rhoe_tot ) * dVol

    if !isnothing(Ham.rhoe_core)
        if Nspin == 2
            Rhoe[:,1] .+= Ham.rhoe_core*0.5
            Rhoe[:,2] .+= Ham.rhoe_core*0.5
        else # should be Nspin == 1
            Rhoe[:,1] .+= Ham.rhoe_core
        end
    end
    #
    epsxc = zeros(Float64, Npoints)
    #
    if Ham.xcfunc == "SCAN"
        # FIXME: no spinpol yet
        calc_epsxc_SCAN!( Ham, psiks, Rhoe, epsxc )
        #
    elseif Ham.xcfunc == "PBE"
        #
        calc_epsxc_PBE!( Ham.xc_calc, Ham.pw, Rhoe, epsxc )
        #
    else
        calc_epsxc_VWN!( Ham.xc_calc, Rhoe, epsxc )
    end
    # Compute the energy
    E_xc = sum( epsxc .* Rhoe ) * dVol
    # Restore
    if !isnothing(Ham.rhoe_core)
        if Nspin == 2
            Rhoe[:,1] .-= Ham.rhoe_core*0.5
            Rhoe[:,2] .-= Ham.rhoe_core*0.5
        else # should be Nspin == 1
            Rhoe[:,1] .-= Ham.rhoe_core
        end
    end
    return E_Ps_loc, E_Hartree, E_xc
end


"""
    calc_E_kin(Ham, psiks)

Compute and return kinetic energy term for a given Hamiltonian `Ham` and BlochWavefunc `psiks`.
"""
function calc_E_kin( Ham::Hamiltonian, psiks::BlochWavefunc )

    Focc = Ham.electrons.Focc
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nstates = Ham.electrons.Nstates
    wk = Ham.pw.gvecw.kpoints.wk
    Nspin = Ham.electrons.Nspin

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
