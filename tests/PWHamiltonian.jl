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
    Rhoe::Array{Float64,1}
    Focc   # not yet determined
end

function PWHamiltonian( pw::PWGrid, atoms::Atoms )
    Npoints = prod(pw.Ns)
    #
    V_Ps_loc = Array{Float64}(Npoints)
    V_Hartree = Array{Float64}(Npoints)
    V_XC = Array{Float64}(Npoints)
    potentials = PotentialsT( V_Ps_loc, V_Hartree, V_XC )
    #
    energies = EnergiesT()
    #
    Rhoe = Array{Float64}(Npoints)
    Focc = nothing
    return PWHamiltonian( pw, potentials, energies, Rhoe, nothing )
end

include("op_K.jl")
include("op_V_loc.jl")
