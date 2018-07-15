# FIXME: psi is not used
function calc_E_xc( Ham::Hamiltonian, psi::Array{ComplexF64,2} )
    Ω = Ham.pw.Ω
    Npoints = prod(Ham.pw.Ns)
    rhoe = Ham.rhoe
    if Ham.xcfunc == "PBE"
        epsxc = calc_epsxc_PBE( Ham.pw, rhoe )
    else
        epsxc = calc_epsxc_VWN( rhoe )
    end
    E_xc = dot( epsxc, rhoe ) * Ω/Npoints
    return E_xc
end


function calc_E_Hartree( Ham::Hamiltonian, psi::Array{ComplexF64,2} )
    potentials = Ham.potentials
    Ω = Ham.pw.Ω
    Npoints = prod(Ham.pw.Ns)
    rhoe = Ham.rhoe
    E_Hartree = 0.5*dot( potentials.Hartree, rhoe ) * Ω/Npoints
    return E_Hartree
end

#
# Use ik and spin
#
function calc_E_Ps_nloc( Ham::Hamiltonian, psiks::Array{Array{ComplexF64,2},1} )

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
    E_ps_NL = 0.0

    betaNL_psi = zeros(ComplexF64,Nstates,NbetaNL)
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        psi = psiks[ikspin]
        betaNL_psi = calc_betaNL_psi( ik, Ham.pspotNL.betaNL, psi )
        for ist = 1:Nstates
            enl1 = 0.0
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
                    enl1 = enl1 + hij*real(conj(betaNL_psi[ist,ibeta])*betaNL_psi[ist,jbeta])
                end
                end
                end # m
                end # l
            end
            E_ps_NL = E_ps_NL + wk[ik]*Focc[ist,ikspin]*enl1
        end
    end
    end

    return E_ps_NL

end


#
# psi is assumed to be already orthonormalized elsewhere
# `potentials` and `Rhoe` are not updated
# Ham is assumed to be already updated at input psi
#
function calc_energies( Ham::Hamiltonian, psiks::Array{Array{ComplexF64,2},1} )

    pw = Ham.pw
    potentials = Ham.potentials
    Focc = Ham.electrons.Focc

    Ω = pw.Ω
    Ns = pw.Ns
    Npoints = prod(Ns)
    dVol = Ω/Npoints
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nstates = Ham.electrons.Nstates
    wk = Ham.pw.gvecw.kpoints.wk
    Nspin = Ham.electrons.Nspin
    
    #
    # Kinetic energy
    #
    E_kin = 0.0
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        psi = psiks[ikspin]
        Kpsi = op_K( Ham, psi )
        for ist = 1:Nstates
            E_kin = E_kin + wk[ik] * Focc[ist,ikspin] * real( dot( psi[:,ist], Kpsi[:,ist] ) )
        end
    end
    end

    Rhoe_total = zeros(Npoints)
    for ispin = 1:Nspin
        Rhoe_total[:] = Rhoe_total[:] + Ham.rhoe[:,ispin]
    end

    E_Hartree = 0.5*dot( potentials.Hartree, Rhoe_total ) * dVol

    E_Ps_loc = dot( potentials.Ps_loc, Rhoe_total ) * dVol

    Rhoe = Ham.rhoe
    if Ham.xcfunc == "PBE"
        epsxc = calc_epsxc_PBE( Ham.pw, Rhoe )
    else
        epsxc = calc_epsxc_VWN( Rhoe )
    end
    E_xc = dot( epsxc, Rhoe_total ) * dVol

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
    energies.Total   = E_kin + E_Ps_loc + E_Ps_nloc + E_Hartree + E_xc + Ham.energies.NN

    return energies
end



# Special use for Nkpt=1 and Nspin=1
# !!! Not extended for general k-dependent wavefunction
#
# psi is assumed to be already orthonormalized elsewhere
# `potentials` and `Rhoe` are not updated
# Ham is assumed to be already updated at input psi
#
function calc_energies( Ham::Hamiltonian, psi::Array{ComplexF64,2} )

    @assert( Ham.pw.gvecw.kpoints.Nkpt == 1 )
    @assert( Ham.electrons.Nspin == 1 )

    pw = Ham.pw
    potentials = Ham.potentials
    Focc = Ham.electrons.Focc

    ik = 1
    ispin = 1

    Ω = pw.Ω
    Ns = pw.Ns
    Npoints = prod(Ns)

    Ngwx = size(psi)[1]  # This should be guaranted by Nkpt = 1
    Nstates = size(psi)[2]

    Kpsi = op_K( Ham, psi )
    E_kin = 0.0
    for ist = 1:Nstates
        E_kin = E_kin + Focc[ist,1] * real( dot( psi[:,ist], Kpsi[:,ist] ) )
    end

    Rhoe = Ham.rhoe

    E_Hartree = 0.5*dot( potentials.Hartree, Rhoe ) * Ω/Npoints

    E_Ps_loc = dot( potentials.Ps_loc, Rhoe ) * Ω/Npoints

    if Ham.xcfunc == "PBE"
        epsxc = calc_epsxc_PBE( Ham.pw, Rhoe )
    else
        epsxc = calc_epsxc_VWN( Rhoe )
    end
    E_xc = dot( epsxc, Rhoe ) * Ω/Npoints

    if Ham.pspotNL.NbetaNL > 0
        E_Ps_nloc = calc_E_Ps_nloc( Ham, psi )
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
    energies.Total   = E_kin + E_Ps_loc + E_Ps_nloc + E_Hartree + E_xc + Ham.energies.NN

    return energies
end
