abstract type AbstractXCCalculator
end

# Use internal XC
struct XCCalculator <: AbstractXCCalculator
end

# Use Libxc
struct LibxcXCCalculator <: AbstractXCCalculator
    x_id::Int64
    c_id::Int64
    Vlapl::Array{Float64,2}
    Vtau::Array{Float64,2}
end

function LibxcXCCalculator( ; x_id=1, c_id=7, is_metagga=false, Npoints=0, Nspin=1 )
    if is_metagga
        @assert Npoints > 0
        Vlapl = zeros(Npoints,Nspin)
        Vtau = zeros(Npoints,Nspin)
        return LibxcXCCalculator(x_id, c_id, Vlapl, Vtau) # XXX use SCAN
    else
        return LibxcXCCalculator(x_id, c_id, zeros(1,1),zeros(1,1))
    end
end

# These types should contain xcfunc field
