import PWDFT: calc_strfact

# Calculate structure factor
function calc_strfact( atpos::Array{Float64,2}, Nspecies::Int64,
    atm2species::Array{Int64,1}, G::CuArray{Float64,2} )
    #
    Ng = size(G)[2]
    Na = size(atpos)[2]
    Sf = CuArrays.zeros(ComplexF64,Ng,Nspecies)

    Nthreads = 256
    Nblocks = ceil(Int64, Ng/Nthreads)

    for ia = 1:Na
        isp = atm2species[ia]
        x = atpos[1,ia]
        y = atpos[2,ia]
        z = atpos[3,ia]
        @cuda threads=Nthreads blocks=Nblocks kernel_calc_strfact!( isp, Ng, x, y, z, G, Sf )
    end
    return Sf
end

function calc_strfact( atoms::Atoms, pw::CuPWGrid )
    return calc_strfact( atoms.positions, atoms.Nspecies, atoms.atm2species, pw.gvec.G )
end

function kernel_calc_strfact!( isp, Ng, x, y, z, G, Sf )
    ig = ( blockIdx().x - 1 ) * blockDim().x + threadIdx().x
    if ig <= Ng
        GX = x*G[1,ig] + y*G[2,ig] + z*G[3,ig]
        Sf[ig,isp] = Sf[ig,isp] + CUDAnative.cos(GX) - im*CUDAnative.sin(GX)
    end
    return
end
