#
# psi is assumed to be already orthonormalized elsewhere
# Potentials and Rhoe are not updated
# Ham is assumed to be already updated at input psi
#
function calc_energies( Ham::PWHamiltonian, psi::Array{Complex128,2} )

    PW = Ham.pw
    Potentials = Ham.potentials
    Focc = Ham.electrons.Focc

    Ω = PW.Ω
    Ns = PW.Ns
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

    E_xc = dot( excVWN(rhoe), rhoe ) * Ω/Npoints

    E_Ps_loc = dot( Potentials.Ps_loc, rhoe ) * Ω/Npoints

    Energies = EnergiesT()
    Energies.Kinetic = E_kin
    Energies.Ps_loc  = E_Ps_loc
    Energies.Hartree = E_Hartree
    Energies.XC      = E_xc
    Energies.NN      = Ham.energies.NN
    Energies.Total   = E_kin + E_Ps_loc + E_Hartree + E_xc + Ham.energies.NN

    return Energies
end
