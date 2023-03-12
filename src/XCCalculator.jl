abstract type AbstractXCCalculator
end

# Use internal XC
struct XCCalculator <: AbstractXCCalculator
end

# Use Libxc
struct LibxcXCCalculator <: AbstractXCCalculator
    x_id::Int64
    c_id::Int64
    need_gradient::Bool # of electron density
    need_KE_dens::Bool
    need_laplacian::Bool # of electron density
    Vlapl::Union{Nothing,Array{Float64,2}}
    Vtau::Union{Nothing,Array{Float64,2}}
end

"""
For now, is_gga and is_metagga need to be set manually
Ideally it should be inferred from x_id and c_id

The current behavior of this constructor:

- is_gga=false and is_metagga=false, then c_id=1 and x_id=7

- is_gga=true and is_metagga=false, then x_id=101 and c_id=130

- is_metagga=true and is_gga=false, then x_id=263 and c_id=267

Note that we x_id and c_id are not used.
"""
function LibxcXCCalculator( ; x_id=1, c_id=7, is_gga=false, is_metagga=false, Npoints=0, Nspin=1 )
    
    if is_metagga && !is_gga

        @assert Npoints > 0
        Vlapl = nothing # Not needed for SCAN
        Vtau = zeros(Npoints,Nspin)
        need_gradient = true
        need_KE_dens = true
        need_laplacian = false # for SCAN
        x_id = 263
        c_id = 267
        return LibxcXCCalculator(x_id, c_id, need_gradient, need_KE_dens, need_laplacian, Vlapl, Vtau)
    
    elseif is_gga && !is_metagga

        Vlapl = nothing
        Vtau = nothing
        need_gradient = true
        need_KE_dens = false
        need_laplacian = false
        x_id = 101
        c_id = 130
        return LibxcXCCalculator(x_id, c_id, need_gradient, need_KE_dens, need_laplacian, Vlapl, Vtau)

    else
        
        Vlapl = nothing
        Vtau = nothing
        need_gradient = false
        need_KE_dens = false
        need_laplacian = false
        return LibxcXCCalculator(x_id, c_id, need_gradient, need_KE_dens, need_laplacian, Vlapl, Vtau)

    end

end

