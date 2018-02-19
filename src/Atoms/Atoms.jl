mutable struct Atoms
    Natoms::Int64
    Nspecies::Int64
    positions::Array{Float64,2}    
    atm2species::Array{Int64,1}    
    symbols::Array{String,1}
    SpeciesSymbols::Array{String,1}  # unique symbols
end

function init_xyz( filexyz::String )
    f = open(filexyz)
end