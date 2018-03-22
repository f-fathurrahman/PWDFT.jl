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
    Potentials = Ham.potentials
    Ω = Ham.pw.Ω
    Npoints = prod(Ham.pw.Ns)
    rhoe = Ham.rhoe
    E_Hartree = 0.5*dot( Potentials.Hartree, rhoe ) * Ω/Npoints
    return E_Hartree
end

#
# psi is assumed to be already orthonormalized elsewhere
# Potentials and Rhoe are not updated
# Ham is assumed to be already updated at input psi
#
function calc_energies( Ham::PWHamiltonian, psi::Array{Complex128,2} )

    pw = Ham.pw
    Potentials = Ham.potentials
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

    E_Hartree = 0.5*dot( Potentials.Hartree, rhoe ) * Ω/Npoints

    E_Ps_loc = dot( Potentials.Ps_loc, rhoe ) * Ω/Npoints

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

    Energies = EnergiesT()
    Energies.Kinetic = E_kin
    Energies.Ps_loc  = E_Ps_loc
    Energies.Ps_nloc = E_Ps_nloc
    Energies.Hartree = E_Hartree
    Energies.XC      = E_xc
    Energies.NN      = Ham.energies.NN
    Energies.Total   = E_kin + E_Ps_loc + E_Ps_nloc + E_Hartree + E_xc + Ham.energies.NN

    return Energies
end
