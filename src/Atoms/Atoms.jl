mutable struct Atoms
    symbols::Array{String,1}
    positions::Array{Float64,2}
end

function init_xyz( filexyz::String )
    f = open(filexyz)
end