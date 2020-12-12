abstract type AbstractXCCalculator
end

# Use internal XC
struct XCCalculator <: AbstractXCCalculator
end

# Use Libxc
struct LibxcXCCalculator <: AbstractXCCalculator
    Vlapl::Array{Float64,2}
    Vtau::Array{Float64,2}
end

function LibxcXCCalculator( ; is_metagga=false, Npoints=0, Nspin=1 )
    if is_metagga
        @assert Npoints > 0
        Vlapl = zeros(Npoints,Nspin)
        Vtau = zeros(Npoints,Nspin)
        return LibxcXCCalculator(Vlapl, Vtau)
    else
        return LibxcXCCalculator(zeros(1,1),zeros(1,1))
    end
end

# These types should contain xcfunc field
