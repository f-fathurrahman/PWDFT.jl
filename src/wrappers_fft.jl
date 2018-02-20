
const USING_FFT_NAIVE = true

# --- Call FFTW3 directly --- #

# XXX Need automatic way to do this
const FFT_SO_PATH = "/home/efefer/WORKS/my_github_repos/PWDFT.jl/src/extlibs/fft3d.so"

function c_G_to_R( Ns::Array{Int,1}, fG::Array{Complex128,1} )
    out = zeros(Complex128,Ns[1]*Ns[2]*Ns[3])
    ccall( (:fftw_inv_fft3d, FFT_SO_PATH), Void,
           (Ptr{Complex128}, Ptr{Complex128}, Int64, Int64, Int64),
            fG, out, Ns[3],Ns[2],Ns[1] )
    return out
end

"""
Multicolumn version
"""
function c_G_to_R( Ns::Array{Int,1}, fG::Array{Complex128,2} )
    Npoints = size(fG)[1]
    Ncol = size(fG)[2]
    out = zeros(Complex128,Npoints,Ncol)
    for ic = 1:Ncol
        out[:,ic] = c_G_to_R( Ns, fG[:,ic] )
        #ccall( (:fftw_inv_fft3d, FFT_SO_PATH), Void,
        #      (Ptr{Complex128}, Ptr{Complex128}, Int64, Int64, Int64),
        #       fG[:,ic], out[:,ic], Ns[3],Ns[2],Ns[1] )
    end
    return out
end

# --- Call FFTW3 directly
function c_R_to_G( Ns::Array{Int,1}, fR::Array{Complex128,1} )
    out = zeros(Complex128,Ns[1]*Ns[2]*Ns[3])
    ccall( (:fftw_fw_fft3d, FFT_SO_PATH), Void,
           (Ptr{Complex128}, Ptr{Complex128}, Int64, Int64, Int64),
            fR, out, Ns[3],Ns[2],Ns[1] )
    return out
end

function c_R_to_G( Ns::Array{Int,1}, fR_::Array{Float64,1} )
    out = zeros(Complex128,Ns[1]*Ns[2]*Ns[3])
    fR = convert(Array{Complex128,1},fR_)
    ccall( (:fftw_fw_fft3d, FFT_SO_PATH), Void,
           (Ptr{Complex128}, Ptr{Complex128}, Int64, Int64, Int64),
            fR, out, Ns[3],Ns[2],Ns[1] )
    return out
end

"""
Multicolumn version
"""
function c_R_to_G( Ns::Array{Int,1}, fR::Array{Complex128,2} )
    Npoints = size(fR)[1]
    Ncol = size(fR)[2]
    out = zeros(Complex128,Npoints,Ncol)
    for ic = 1:Ncol
        out[:,ic] = c_R_to_G( Ns, fR[:,ic] )
        #ccall( (:fftw_fw_fft3d, FFT_SO_PATH), Void,
        #       (Ptr{Complex128}, Ptr{Complex128}, Int64, Int64, Int64),
        #        fR[:,ic], out[:,ic], Ns[3],Ns[2],Ns[1] )
    end
    return out
end

#
# Naive version of calling FFT
#
if USING_FFT_NAIVE

function G_to_R( Ns::Array{Int,1}, fG::Array{Complex128,1} )
    out = reshape( ifft( reshape(fG,Ns[1],Ns[2],Ns[3]) ),size(fG) )
end

function G_to_R( Ns::Array{Int,1}, fG::Array{Complex128,2} )
    Npoints = prod(Ns)
    out = zeros( Complex128, size(fG) ) # Is this safe?
    for ic = 1:size(fG,2)
        out[:,ic] = reshape( ifft( reshape(fG[:,ic],Ns[1],Ns[2],Ns[3]) ), Npoints )
    end
    return out
end

function R_to_G( Ns::Array{Int,1}, fR::Array{Complex128,1} )
    out = reshape( fft( reshape(fR,Ns[1],Ns[2],Ns[3]) ), size(fR) )
end

# In case we forget to convert the input, we convert it in this version
function R_to_G( Ns::Array{Int,1}, fR_::Array{Float64,1} )
    fR = convert(Array{Complex128,1},fR_)
    out = reshape( fft( reshape(fR,Ns[1],Ns[2],Ns[3]) ), size(fR) )
end

function R_to_G( Ns::Array{Int,1}, fR::Array{Complex128,2} )
    Npoints = prod(Ns)
    Ncol = size(fR,2)
    out = zeros( Complex128, size(fR) )
    for ic = 1:Ncol
        out[:,ic] = reshape( fft( reshape(fR[:,ic],Ns[1],Ns[2],Ns[3]) ), Npoints )
    end
    return out
end

else   # not USING_FFT_NAIVE

const R_to_G = c_R_to_G
const G_to_R = c_G_to_R

end   # USING_FFT_NAIVE
