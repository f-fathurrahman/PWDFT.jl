function calc_E_xc( Ham::PWHamiltonian, psi::Array{Complex128,2} )
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


function calc_E_Hartree( Ham::PWHamiltonian, psi::Array{Complex128,2} )
    potentials = Ham.potentials
    Ω = Ham.pw.Ω
    Npoints = prod(Ham.pw.Ns)
    rhoe = Ham.rhoe
    E_Hartree = 0.5*dot( potentials.Hartree, rhoe ) * Ω/Npoints
    return E_Hartree
end

"""
Use ik
"""
function calc_E_Ps_nloc( Ham::PWHamiltonian, psik::Array{Array{Complex128,2},1} )

    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    Natoms = Ham.atoms.Natoms
    atm2species = Ham.atoms.atm2species
    Pspots = Ham.pspots
    prj2beta = Ham.pspotNL.prj2beta
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    wk = Ham.pw.gvecw.kpoints.wk
    NbetaNL = Ham.pspotNL.NbetaNL

    # calculate E_NL
    E_ps_NL = 0.0

    betaNL_psi = zeros(Complex128,Nstates,NbetaNL)
    for ik = 1:Nkpt
        psi = psik[ik]
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
            E_ps_NL = E_ps_NL + wk[ik]*Focc[ist]*enl1
        end
    end

    return E_ps_NL

end


#
# psi is assumed to be already orthonormalized elsewhere
# `potentials` and `Rhoe` are not updated
# Ham is assumed to be already updated at input psi
#
function calc_energies( Ham::PWHamiltonian, psik::Array{Array{Complex128,2},1} )

    pw = Ham.pw
    potentials = Ham.potentials
    Focc = Ham.electrons.Focc

    Ω = pw.Ω
    Ns = pw.Ns
    Npoints = prod(Ns)
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nstates = Ham.electrons.Nstates
    wk = Ham.pw.gvecw.kpoints.wk
    
    #
    # Kinetic energy
    #
    E_kin = 0.0
    for ik = 1:Nkpt
        Ham.ik = ik
        psi = psik[ik]
        Kpsi = op_K( Ham, psi )
        for ist = 1:Nstates
            E_kin = E_kin + wk[ik] * Focc[ist] * real( dot( psi[:,ist], Kpsi[:,ist] ) )
        end
    end

    rhoe = Ham.rhoe

    E_Hartree = 0.5*dot( potentials.Hartree, rhoe ) * Ω/Npoints

    E_Ps_loc = dot( potentials.Ps_loc, rhoe ) * Ω/Npoints

    if Ham.xcfunc == "PBE"
        epsxc = calc_epsxc_PBE( Ham.pw, rhoe )
    else
        epsxc = calc_epsxc_VWN( rhoe )
    end
    E_xc = dot( epsxc, rhoe ) * Ω/Npoints

    if Ham.pspotNL.NbetaNL > 0
        E_Ps_nloc = calc_E_Ps_nloc( Ham, psik )
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



#
# psi is assumed to be already orthonormalized elsewhere
# `potentials` and `Rhoe` are not updated
# Ham is assumed to be already updated at input psi
#
function calc_energies( Ham::PWHamiltonian, psi::Array{Complex128,2} )

    pw = Ham.pw
    potentials = Ham.potentials
    Focc = Ham.electrons.Focc

    Ω = pw.Ω
    Ns = pw.Ns
    Npoints = prod(Ns)

    Ngwx = size(psi)[1]
    Nstates = size(psi)[2]
    #
    Kpsi = op_K( Ham, psi )
    E_kin = 0.0
    for ist = 1:Nstates
        E_kin = E_kin + Focc[ist] * real( dot( psi[:,ist], Kpsi[:,ist] ) )
    end

    rhoe = Ham.rhoe

    E_Hartree = 0.5*dot( potentials.Hartree, rhoe ) * Ω/Npoints

    E_Ps_loc = dot( potentials.Ps_loc, rhoe ) * Ω/Npoints

    if Ham.xcfunc == "PBE"
        epsxc = calc_epsxc_PBE( Ham.pw, rhoe )
    else
        epsxc = calc_epsxc_VWN( rhoe )
    end
    E_xc = dot( epsxc, rhoe ) * Ω/Npoints

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
