mutable struct CuPotentials
    Ps_loc::CuArray{Float64,1}
    Hartree::CuArray{Float64,1}
    XC::CuArray{Float64,2}
    Total::CuArray{Float64,2}
end
