mutable struct Potentials
    Ps_loc::Array{Float64,1}
    Hartree::Array{Float64,1}
    XC::Array{Float64,2}  # spin dependent
end