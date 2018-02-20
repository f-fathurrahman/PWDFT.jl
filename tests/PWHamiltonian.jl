mutable struct EnergiesT
    Total::Float64
    Kinetic::Float64
    Ps_loc::Float64
    Hartree::Float64
    XC::Float64
    NN::Float64
end

# Default: all zeroes
function EnergiesT()
    return EnergiesT(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

import Base.println
function println( Energies::EnergiesT )
    @printf("Kinetic energy: %18.10f\n", Energies.Kinetic )
    @printf("Ps_loc  energy: %18.10f\n", Energies.Ps_loc )
    @printf("Hartree energy: %18.10f\n", Energies.Hartree )
    @printf("XC      energy: %18.10f\n", Energies.XC )
    @printf("NN      energy: %18.10f\n", Energies.NN )
    @printf("----------------------------------\n")
    @printf("Total   energy: %18.10f\n", Energies.Total )
end

mutable struct PotentialsT
    Ps_loc::Array{Float64,1}
    Hartree::Array{Float64,1}
    XC::Array{Float64,1}
end

mutable struct PWHamiltonian
    pw::PWGrid
    potentials::PotentialsT
    energies::EnergiesT
    rhoe::Array{Float64,1}
    focc   # not yet determined
end

include("calc_strfact.jl")
include("init_V_coulomb_G.jl")

function PWHamiltonian( pw::PWGrid, atoms::Atoms )
    Npoints = prod(pw.Ns)
    #
    strf = calc_strfact( atoms, pw )
    V_Ps_loc = init_V_coulomb_G( pw, strf, [1.0] )
    #
    V_Hartree = zeros(Float64,Npoints)
    V_XC = zeros(Float64,Npoints)
    potentials = PotentialsT( V_Ps_loc, V_Hartree, V_XC )
    #
    energies = EnergiesT()
    #
    rhoe = zeros(Float64,Npoints)
    return PWHamiltonian( pw, potentials, energies, rhoe, nothing )
end

include("op_K.jl")
include("op_V_loc.jl")
include("op_H.jl")

include("Poisson_solve.jl")
include("LDA_VWN.jl")

function update!(Ham::PWHamiltonian, rhoe::Array{Float64,1})
    Ham.rhoe = rhoe
    Ham.potentials.Hartree = real( G_to_R( Ham.pw.Ns, Poisson_solve(Ham.pw, rhoe) ) )
    Ham.potentials.XC = excVWN( rhoe ) + rhoe .* excpVWN( rhoe )
end
