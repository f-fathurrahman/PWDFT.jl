#
# In-place version, input 3d data as 3d array
#
function G_to_R!( pw, fG::Array{ComplexF64,3} )
    pw.planbw*fG
    return
end

function R_to_G!( pw, fR::Array{ComplexF64,3} )
    pw.planfw*fR
    return
end

#
# In-place version, input 3d data as column vector
#
function G_to_R!( pw, fG::Vector{ComplexF64} )
    ff = reshape(fG, pw.Ns)
    pw.planbw*ff
    return
end

function R_to_G!( pw, fR::Vector{ComplexF64} )
    ff = reshape(fR, pw.Ns)
    pw.planfw*ff
    return
end


#
# Return a new array
#
function G_to_R( pw, fG::Vector{ComplexF64} )
    ff = copy(fG)
    ff = reshape(ff, pw.Ns)
    pw.planbw*ff
    return reshape(ff, prod(pw.Ns))
end

function R_to_G( pw, fR::Vector{ComplexF64} )
    ff = copy(fR)
    ff = reshape(fR, pw.Ns)
    pw.planfw*ff
    return reshape(ff, prod(pw.Ns))
end


#
# used in Poisson solver
#
function R_to_G( pw, fR_::Vector{Float64} )
    fR = convert(Array{ComplexF64,1}, fR_) # This will make a copy
    ff = reshape(fR, pw.Ns)
    pw.planfw*ff
    return reshape(ff, prod(pw.Ns))
end


#
# Used in calc_rhoe
#
function G_to_R!( pw, fG::Matrix{ComplexF64} )
    plan = pw.planbw
    for i in 1:size(fG,2)
        @views ff = reshape(fG[:,i], pw.Ns)
        pw.planbw*ff
    end
    return
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
        @views out[:,ic] = reshape( ifft( reshape(fG[:,ic],Ns) ), Npoints )
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