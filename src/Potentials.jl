"""
The type for representing local part of Kohn-Sham Hamiltonians in real space:
- `Ps_loc`: local pseudopotential
- `Hartree`: Hartree potential
- `XC`: XC potential
"""
mutable struct Potentials
    Ps_loc::Array{Float64,1}
    Hartree::Array{Float64,1}
    XC::Array{Float64,2}  # spin dependent
    Total::Array{Float64,2}
    TotalSmooth::Union{Array{Float64,2},Nothing} # smooth potential
end