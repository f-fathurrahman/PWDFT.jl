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
