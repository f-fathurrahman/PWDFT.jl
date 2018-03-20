mutable struct EnergiesT
    Total::Float64
    Kinetic::Float64
    Ps_loc::Float64
    Ps_nloc::Float64
    Hartree::Float64
    XC::Float64
    NN::Float64
end

# Default: all zeroes
function EnergiesT()
    return EnergiesT(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

import Base.println
function println( Energies::EnergiesT )
    @printf("\n")
    @printf("Kinetic    energy: %18.10f\n", Energies.Kinetic )
    @printf("Ps_loc     energy: %18.10f\n", Energies.Ps_loc )
    @printf("Ps_nloc    energy: %18.10f\n", Energies.Ps_nloc )
    @printf("Hartree    energy: %18.10f\n", Energies.Hartree )
    @printf("XC         energy: %18.10f\n", Energies.XC )
    @printf("-------------------------------------\n")
    E_elec = Energies.Kinetic + Energies.Ps_loc + Energies.Ps_nloc +
             Energies.Hartree + Energies.XC
    @printf("Electronic energy: %18.10f\n", E_elec)
    @printf("NN         energy: %18.10f\n", Energies.NN )
    @printf("-------------------------------------\n")
    @printf("Total      energy: %18.10f\n", Energies.Total )
end

mutable struct PotentialsT
    Ps_loc::Array{Float64,1}
    Hartree::Array{Float64,1}
    XC::Array{Float64,1}
end

include("PsPotNL.jl")

mutable struct PWHamiltonian
    pw::PWGrid
    potentials::PotentialsT
    energies::EnergiesT
    rhoe::Array{Float64,1}
    electrons::ElectronsInfo
    atoms::Atoms
    pspots::Array{PsPot_GTH,1}
    pspotNL::PsPotNL
end


function PWHamiltonian( atoms::Atoms, pspfiles::Array{String,1},
                        ecutwfc::Float64, LatVecs::Array{Float64,2} )
    # Initialize plane wave grids
    pw = PWGrid( ecutwfc, LatVecs )
    println(pw)

    Nspecies = atoms.Nspecies
    if Nspecies != size(pspfiles)[1]
        @printf("ERROR length of pspfiles is not equal to %d\n", Nspecies)
        exit()
    end

    Npoints = prod(pw.Ns)
    Ω = pw.Ω
    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r

    strf = calc_strfact( atoms, pw )

    #
    # Initialize pseudopotentials and local potentials
    #
    Vg = zeros(Complex128, Npoints)
    V_Ps_loc = zeros(Float64, Npoints)

    Pspots = Array{PsPot_GTH}(Nspecies)

    for isp = 1:Nspecies
        Pspots[isp] = PsPot_GTH( pspfiles[isp] )
        psp = Pspots[isp]
        println(psp)
        for ig = 1:Ng
            ip = idx_g2r[ig]
            Vg[ip] = strf[ig,isp] * eval_Vloc_G( psp, G2[ig] ) / Ω
        end
        #
        V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(pw, Vg) ) * Npoints
    end

    # other potential terms are set to zero
    V_Hartree = zeros( Float64, Npoints )
    V_XC = zeros( Float64, Npoints )
    potentials = PotentialsT( V_Ps_loc, V_Hartree, V_XC )
    #
    energies = EnergiesT()
    #
    rhoe = zeros( Float64, Npoints )

    electrons = ElectronsInfo( atoms, Pspots )
    println(electrons)

    # NL pseudopotentials
    pspotNL = PsPotNL( pw, atoms, Pspots, check_norm=false )

    return PWHamiltonian( pw, potentials, energies, rhoe, electrons, atoms, Pspots, pspotNL )
end


#
# No pspfiles given. Use Coulomb potential (all electrons)
#
function PWHamiltonian( atoms::Atoms, ecutwfc::Float64, LatVecs::Array{Float64,2} )

    #
    # Initialize plane wave grids
    #
    pw = PWGrid( ecutwfc, LatVecs )
    println(pw)

    Nspecies = atoms.Nspecies

    Npoints = prod(pw.Ns)

    # dummy pspots
    Pspots = Array{PsPot_GTH}(1)
    Pspots[1] = PsPot_GTH()

    #
    # Coulomb potential as local potential
    #
    Zvals = get_Zatoms(atoms)
    strf = calc_strfact( atoms, pw )
    V_Ps_loc = init_V_coulomb_G( pw, strf, Zvals )

    # other potentials are set to zero
    V_Hartree = zeros( Float64, Npoints )
    V_XC = zeros( Float64, Npoints )
    potentials = PotentialsT( V_Ps_loc, V_Hartree, V_XC )
    #
    energies = EnergiesT()
    #
    rhoe = zeros( Float64, Npoints )

    electrons = ElectronsInfo( atoms, Zvals )
    println(electrons)

    rhoe = zeros(Float64,Npoints)

    pspotNL = PsPotNL()

    return PWHamiltonian( pw, potentials, energies, rhoe, electrons, atoms, Pspots, pspotNL )
end



include("op_K.jl")
include("op_V_loc.jl")
include("op_V_Ps_nloc.jl")
include("op_H.jl")

include("Poisson_solve.jl")
include("LDA_VWN.jl")


"""
Given rhoe in real space, update Ham.rhoe, Hartree and XC potentials.
"""
function update!(Ham::PWHamiltonian, rhoe::Array{Float64,1})
    Ham.rhoe = rhoe
    Ham.potentials.Hartree = real( G_to_R( Ham.pw, Poisson_solve(Ham.pw, rhoe) ) )
    #Ham.potentials.XC = excVWN( rhoe ) + rhoe .* excpVWN( rhoe )
    Ham.potentials.XC = calc_Vxc_VWN( rhoe )
end


function update!( Ham::PWHamiltonian, atoms::Atoms, strf::Array{Complex128,2}, pspfiles::Array{String,1} )

    Nspecies = atoms.Nspecies
    if Nspecies != size(pspfiles)[1]
        @printf("ERROR length of pspfiles is not equal to %d\n", Nspecies)
        exit()
    end

    pw = Ham.pw

    Npoints = prod(pw.Ns)
    Ω = pw.Ω
    G2 = pw.gvec.G2

    Vg = zeros(Complex128, Npoints)
    V_Ps_loc = zeros(Float64, Npoints)

    pw = Ham.pw
    for isp = 1:Nspecies
        psp = PsPot_GTH( pspfiles[isp] )
        println(psp)
        for ig = 1:Npoints
            Vg[ig] = strf[ig,isp] * eval_Vloc_G( psp, G2[ig], Ω )
        end
        #
        V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(pw, Vg) ) * Npoints
    end

    # Update
    Ham.potentials.Ps_loc = V_Ps_loc
    return

end
