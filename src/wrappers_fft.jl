#
# Using plan_fft and plan_ifft
#
function G_to_R( pw::PWGrid, fG::Array{ComplexF64,1} )
    Ns = pw.Ns
    Npoints = prod(Ns)
    plan = pw.planbw
    out = reshape( plan*reshape(fG,Ns), Npoints )
    return out
end

# without reshape
function G_to_R( pw::PWGrid, fG::Array{ComplexF64,3} )
    return pw.planbw*fG
end

function G_to_R( pw::PWGrid, fG::Array{ComplexF64,2} )
    Ns = pw.Ns
    Npoints = prod(Ns)
    plan = pw.planbw
    out = zeros( ComplexF64, size(fG) )
    for ic = 1:size(fG,2)
        out[:,ic] = reshape( plan*reshape(fG[:,ic],Ns), Npoints )
    end
    return out
end

function R_to_G( pw::PWGrid, fR::Array{ComplexF64,1} )
    Ns = pw.Ns
    Npoints = prod(Ns)
    plan = pw.planfw
    out = reshape( plan*reshape(fR,Ns), Npoints )
    return out
end

function R_to_G( pw::PWGrid, fR::Array{ComplexF64,3} )
    return pw.planfw*fR
end

function R_to_G( pw::PWGrid, fR_::Array{Float64,1} )
    Ns = pw.Ns
    Npoints = prod(Ns)
    plan = pw.planfw
    fR = convert(Array{ComplexF64,1},fR_)
    out = reshape( plan*reshape(fR,Ns), Npoints )
    return out
end

function R_to_G( pw::PWGrid, fR::Array{ComplexF64,2} )
    Ns = pw.Ns
    plan = pw.planfw
    Npoints = prod(Ns)
    Ncol = size(fR,2)
    out = zeros( ComplexF64, size(fR) )
    for ic = 1:Ncol
        out[:,ic] = reshape( plan*reshape(fR[:,ic],Ns), Npoints )
    end
    return out
end

#
# Using fft and ifft directly
#

function G_to_R( Ns::Tuple{Int64,Int64,Int64}, fG::Array{ComplexF64,1} )
    out = reshape( ifft( reshape(fG,Ns) ),size(fG) )
end

# without reshape
function G_to_R( Ns::Tuple{Int64,Int64,Int64}, fG::Array{ComplexF64,3} )
    out = ifft(fG)
end

# multicolumn
function G_to_R( Ns::Tuple{Int64,Int64,Int64}, fG::Array{ComplexF64,2} )
    Npoints = prod(Ns)
    out = zeros( ComplexF64, size(fG) ) # Is this safe?
    for ic = 1:size(fG)[2]
        out[:,ic] = reshape( ifft( reshape(fG[:,ic],Ns) ), Npoints )
    end
    return out
end

function R_to_G( Ns::Tuple{Int64,Int64,Int64}, fR::Array{ComplexF64,1} )
    out = reshape( fft( reshape(fR,Ns) ), size(fR) )
end

# without reshape
function R_to_G( Ns::Tuple{Int64,Int64,Int64}, fR::Array{ComplexF64,3} )
    out = fft(fR)
end

# In case we forget to convert the input, we convert it in this version
function R_to_G( Ns::Tuple{Int64,Int64,Int64}, fR_::Array{Float64,1} )
    fR = convert(Array{ComplexF64,1},fR_)
    out = reshape( fft( reshape(fR,Ns) ), size(fR) )
end
    
function R_to_G( Ns::Tuple{Int64,Int64,Int64}, fR::Array{ComplexF64,2} )
    Npoints = prod(Ns)
    Ncol = size(fR)[2]
    out = zeros( ComplexF64, size(fR) )
    for ic = 1:Ncol
        out[:,ic] = reshape( fft( reshape(fR[:,ic],Ns) ), Npoints )
    end
    return out
end