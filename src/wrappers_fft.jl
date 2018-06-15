
const USING_FFT_JULIA = true

# --- Call FFTW3 directly --- #

# XXX Need automatic way to do this
const FFT_SO_PATH = "/home/efefer/WORKS/my_github_repos/PWDFT.jl/src/extlibs/fft3d.so"

function c_G_to_R( Ns::Tuple{Int64,Int64,Int64}, fG::Array{ComplexF64,1} )
    out = zeros(ComplexF64,Ns[1]*Ns[2]*Ns[3])
    ccall( (:fftw_inv_fft3d, FFT_SO_PATH), Nothing,
           (Ptr{ComplexF64}, Ptr{ComplexF64}, Int64, Int64, Int64),
            fG, out, Ns[3],Ns[2],Ns[1] )
    return out
end

"""
Multicolumn version
"""
function c_G_to_R( Ns::Tuple{Int64,Int64,Int64}, fG::Array{ComplexF64,2} )
    Npoints = size(fG)[1]
    Ncol = size(fG)[2]
    out = zeros(ComplexF64,Npoints,Ncol)
    for ic = 1:Ncol
        out[:,ic] = c_G_to_R( Ns, fG[:,ic] )
    end
    return out
end

# --- Call FFTW3 directly
function c_R_to_G( Ns::Tuple{Int64,Int64,Int64}, fR::Array{ComplexF64,1} )
    out = zeros(ComplexF64,Ns[1]*Ns[2]*Ns[3])
    ccall( (:fftw_fw_fft3d, FFT_SO_PATH), Nothing,
           (Ptr{ComplexF64}, Ptr{ComplexF64}, Int64, Int64, Int64),
            fR, out, Ns[3],Ns[2],Ns[1] )
    return out
end

function c_R_to_G( Ns::Tuple{Int64,Int64,Int64}, fR_::Array{Float64,1} )
    out = zeros(ComplexF64,Ns[1]*Ns[2]*Ns[3])
    fR = convert(Array{ComplexF64,1},fR_)
    ccall( (:fftw_fw_fft3d, FFT_SO_PATH), Nothing,
           (Ptr{ComplexF64}, Ptr{ComplexF64}, Int64, Int64, Int64),
            fR, out, Ns[3],Ns[2],Ns[1] )
    return out
end

"""
Multicolumn version
"""
function c_R_to_G( Ns::Tuple{Int64,Int64,Int64}, fR::Array{ComplexF64,2} )
    Npoints = size(fR)[1]
    Ncol = size(fR)[2]
    out = zeros(ComplexF64,Npoints,Ncol)
    for ic = 1:Ncol
        out[:,ic] = c_R_to_G( Ns, fR[:,ic] )
    end
    return out
end



"""Using plan_fft and plan_ifft"""


function G_to_R( pw::PWGrid, fG::Array{ComplexF64,1} )
    Ns = pw.Ns
    Npoints = prod(Ns)
    plan = pw.planbw
    out = reshape( plan*reshape(fG,Ns), Npoints )
    return out
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
# Call FFTW from Julia
#
if USING_FFT_JULIA

function G_to_R( Ns::Tuple{Int64,Int64,Int64}, fG::Array{ComplexF64,1} )
    out = reshape( ifft( reshape(fG,Ns) ),size(fG) )
end

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

else   # not USING_FFT_JULIA

const R_to_G = c_R_to_G
const G_to_R = c_G_to_R

end   # USING_FFT_JULIA
