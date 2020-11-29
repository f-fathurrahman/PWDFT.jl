## TODO: Merge with wrappers_fft

#
# Backward transform
#

# fG is assumed to have size Npoints=Ns[1]*Ns[2]*Ns[3]
function G_to_R( pw::PWGridGamma, fG::Array{ComplexF64,1} )
    Ns = pw.Ns
    Npoints = prod(Ns)
    plan = pw.planbw
    out = reshape( plan*reshape(fG,Ns), Npoints )
    return out
end

# without reshape
function G_to_R( pw::PWGridGamma, fG::Array{ComplexF64,3} )
    return pw.planbw*fG
end

# fG is assumed to have size Npoints=Ns[1]*Ns[2]*Ns[3]
function G_to_R!( pw::PWGridGamma, fG::Array{ComplexF64,1} )
    Ns = pw.Ns
    Npoints = prod(Ns)
    plan = pw.planbw
    @views fG[:] = reshape( plan*reshape(fG[:],Ns), Npoints )
    return
end

# without fG
function G_to_R!( pw::PWGridGamma, fG::Array{ComplexF64,3} )
    @views fG[:,:,:] = pw.planbw*fG[:,:,:]
    return
end

#
# Forward transform
#

function R_to_G( pw::PWGridGamma, fR::Array{ComplexF64,1} )
    Ns = pw.Ns
    Npoints = prod(Ns)
    plan = pw.planfw
    out = reshape( plan*reshape(fR,Ns), Npoints )
    return out
end

function R_to_G( pw::PWGridGamma, fR::Array{ComplexF64,3} )
    return pw.planfw*fR
end

function R_to_G!( pw::PWGridGamma, fR::Array{ComplexF64,1} )
    Ns = pw.Ns
    plan = pw.planfw
    Npoints = prod(Ns)
    @views fR[:] = reshape( plan*reshape(fR[:],Ns), Npoints )
    return
end

function R_to_G!( pw::PWGridGamma, fR::Array{ComplexF64,3} )
    plan = pw.planfw
    @views fR[:,:,:] = pw.planfw*fR[:,:,:]
    return
end