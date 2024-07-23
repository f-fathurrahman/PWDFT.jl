"""
The type for set of G-vectors for describing density and potentials
"""
struct GVectors
    Ng::Int64
    G::Array{Float64,2}
    G2::Array{Float64,1}
    idx_g2r::Array{Int64,1}
    G2_shells::Array{Float64,1}
    idx_g2shells::Array{Int64,1}
    idx_g2miller::Vector{Tuple{Int64,Int64,Int64}} #::Array{Int64,2}
end

"""
The type for set of G-vectors for describing wave function.
"""
struct GVectorsW
    Ngwx::Int64          # maximum(Ngk)
    Ngw::Array{Int64,1}  # no of GvectorsW for each kpoints
    idx_gw2g::Array{Array{Int64,1},1}
    idx_gw2r::Array{Array{Int64,1},1}
    kpoints::KPoints
end

const PLANFW_TYPE = typeof(plan_fft!(zeros(ComplexF64,(1,1,1)), flags=FFTW.MEASURE))
const PLANBW_TYPE = typeof(plan_ifft!(zeros(ComplexF64,(1,1,1)), flags=FFTW.MEASURE))

"""
The type for describing plane wave basis set for a given periodic unit cell.
"""
struct PWGrid
    using_dual_grid::Bool
    ecutwfc::Float64
    ecutrho::Float64
    Ns::Tuple{Int64,Int64,Int64}
    Nss::Union{Tuple{Int64,Int64,Int64},Nothing}
    LatVecs::Array{Float64,2}
    RecVecs::Array{Float64,2}
    CellVolume::Float64
    gvec::GVectors
    gvecs::Union{GVectors,Nothing}
    gvecw::GVectorsW
    planfw::PLANFW_TYPE
    planbw::PLANBW_TYPE
    planfws::Union{PLANFW_TYPE,Nothing}
    planbws::Union{PLANBW_TYPE,Nothing}
end

# Adding type annotation to planfw and planbw can reduce number of allocations

"""
    pw = PWGrid( ecutwfc, LatVecs [, kpoints] )

Create an instance of `PWGrid` given the following inputs:

- `ecutwfc`: cutoff energy for wave function
- `LatVecs`: unit cell lattice vectors.
  Convention used: The three lattice vectors v1, v2, and v3 are assumed to be
  arranged **by column**.
- `kpoints`: optional, an instance of `KPoints`.
"""
function PWGrid(
    ecutwfc::Float64,
    LatVecs::Array{Float64,2};
    kpoints::Union{Nothing,KPoints}=nothing,
    Ns_::Tuple{Int64,Int64,Int64}=(0,0,0),
    dual::Float64=4.0
)

    @assert dual >= 4.0

    using_dual_grid = dual > 4.0

    ecutrho = dual*ecutwfc
    RecVecs = 2*pi*inv(Matrix(LatVecs'))
    # XXX RecVecs are stored by rows
    CellVolume = abs(det(LatVecs))
    #
    LatVecsLen = Array{Float64}(undef,3)
    LatVecsLen[1] = norm(LatVecs[:,1])
    LatVecsLen[2] = norm(LatVecs[:,2])
    LatVecsLen[3] = norm(LatVecs[:,3])

    Ns1 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[1]/pi ) + 1
    Ns2 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[2]/pi ) + 1
    Ns3 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[3]/pi ) + 1

    if any(Ns_ .== 0)
        Ns1 = good_fft_order(Ns1)
        Ns2 = good_fft_order(Ns2)
        Ns3 = good_fft_order(Ns3)
        Ns = (Ns1,Ns2,Ns3)
    else
        Ns = Ns_[:]
    end

    Npoints = prod(Ns)
    
    gvec = init_gvec( Ns, RecVecs, ecutrho )
    planfw = plan_fft!( zeros(ComplexF64,Ns), flags=FFTW.MEASURE )
    planbw = plan_ifft!( zeros(ComplexF64,Ns), flags=FFTW.MEASURE )


    # If dual is larger than 4 then we need to allocate a smooth grid:
    # gvecs and its grid size, and the FFR plans associated with it
    #
    # "Smooth" here means that the functions described by the grid are assumed
    # to be smooth enough such that 4*ecutwfc is enough.
    #
    # This is to be contrasted with "dense" grid.
    #
    # Probably it is better to name is a "coarse" grid.
    #
    # The usual gvec will be treated as "dense" grid.
    if using_dual_grid
        # Use 4*ecutwfc as the cutoff
        Ns1 = 2*round( Int64, sqrt(4*ecutwfc/2)*LatVecsLen[1]/pi ) + 1
        Ns2 = 2*round( Int64, sqrt(4*ecutwfc/2)*LatVecsLen[2]/pi ) + 1
        Ns3 = 2*round( Int64, sqrt(4*ecutwfc/2)*LatVecsLen[3]/pi ) + 1
        Ns1 = good_fft_order(Ns1)
        Ns2 = good_fft_order(Ns2)
        Ns3 = good_fft_order(Ns3)
        Nss = (Ns1,Ns2,Ns3)
        # XXX Setting Nss from the constructor is not yet supported
        #
        # TODO: simply copy from gvec instead of calling this function again
        #       Need to calculate Ngs (no. of G-vectors for smooth grid)
        #       gvecs will be a subset of gvec
        gvecs = init_gvec( Nss, RecVecs, 4*ecutwfc )
        planfws = plan_fft!( zeros(ComplexF64,Nss), flags=FFTW.MEASURE )
        planbws = plan_ifft!( zeros(ComplexF64,Nss), flags=FFTW.MEASURE )
    else
        Nss = nothing
        gvecs = nothing
        planfws = nothing
        planbws = nothing
    end

    if kpoints == nothing
        kpoints = KPoints( 1, (1,1,1), zeros(3,1), [1.0], RecVecs )
    end

    if using_dual_grid
        gvecw = init_gvecw( ecutwfc, gvecs, kpoints )
    else
        gvecw = init_gvecw( ecutwfc, gvec, kpoints )
    end

    return PWGrid(
        using_dual_grid,
        ecutwfc, ecutrho, Ns, Nss, LatVecs, RecVecs, CellVolume,
        gvec, gvecs, gvecw, planfw, planbw, planfws, planbws
    )
end

"""
Flip real space indices for FFT.
"""
function mm_to_nn(mm::Int64,S::Int64)
    if mm > S/2
        return mm - S
    else
        return mm
    end
end

"""
Calculates number of G-vectors satisfying |G|^2 <= 2*ecutrho.
This function is used by function `init_gvec`.
"""
function calc_Ng( Ns, RecVecs, ecutrho )
    ig = 0
    Ng = 0
    #
    G = zeros(Float64,3)
    #
    for k in 0:Ns[3]-1
    for j in 0:Ns[2]-1
    for i in 0:Ns[1]-1
        ig = ig + 1
        gi = mm_to_nn( i, Ns[1] )
        gj = mm_to_nn( j, Ns[2] )
        gk = mm_to_nn( k, Ns[3] )
        G[1] = RecVecs[1,1]*gi + RecVecs[1,2]*gj + RecVecs[1,3]*gk
        G[2] = RecVecs[2,1]*gi + RecVecs[2,2]*gj + RecVecs[2,3]*gk
        G[3] = RecVecs[3,1]*gi + RecVecs[3,2]*gj + RecVecs[3,3]*gk
        G2 = G[1]^2 + G[2]^2 + G[3]^2
        if 0.5*G2 <= ecutrho
            Ng = Ng + 1
        end
    end
    end
    end
    return Ng
end


"""
Creates an instance of `GVectors`, given the following inputs:

- `Ns`: sampling points
- `RecVecs`: reciprocal lattice vectors
- `ecutrho`: cutoff energy (in hartree)
"""
function init_gvec( Ns, RecVecs, ecutrho )

    Ng = calc_Ng( Ns, RecVecs, ecutrho )

    G_temp = zeros(Float64,3)

    G  = Array{Float64}(undef,3,Ng)
    G2 = Array{Float64}(undef,Ng)
    idx_g2r = Array{Int64}(undef,Ng)
    idx_g2miller = Vector{Tuple{Int64,Int64,Int64}}(undef,Ng) # Array{Int64,2}(undef,3,Ng)

    ig = 0
    ip = 0
    for k in 0:Ns[3]-1
    for j in 0:Ns[2]-1
    for i in 0:Ns[1]-1
        ip = ip + 1
        gi = mm_to_nn( i, Ns[1] )
        gj = mm_to_nn( j, Ns[2] )
        gk = mm_to_nn( k, Ns[3] )
        # XXX RecVecs are stored by rows
        G_temp[1] = RecVecs[1,1]*gi + RecVecs[1,2]*gj + RecVecs[1,3]*gk
        G_temp[2] = RecVecs[2,1]*gi + RecVecs[2,2]*gj + RecVecs[2,3]*gk
        G_temp[3] = RecVecs[3,1]*gi + RecVecs[3,2]*gj + RecVecs[3,3]*gk
        G2_temp = G_temp[1]^2 + G_temp[2]^2 + G_temp[3]^2
        if 0.5*G2_temp <= ecutrho
            ig = ig + 1
            G[:,ig] = G_temp[:]
            G2[ig] = G2_temp
            idx_g2r[ig] = ip
            idx_g2miller[ig] = (gi, gj, gk)
        end
    end
    end
    end

    # if sorted
    idx_sorted = sortperm(G2)
    G = G[:,idx_sorted]
    G2 = G2[idx_sorted]
    idx_g2r = idx_g2r[idx_sorted]
    idx_g2miller = idx_g2miller[idx_sorted]

    G2_shells, idx_g2shells = init_Gshells( G2 )

    return GVectors( Ng, G, G2, idx_g2r, G2_shells, idx_g2shells, idx_g2miller )

end


function init_Gshells( G2_sorted::Array{Float64,1} )

    eps8 = 1e-8

    Ng = length(G2_sorted)

    ngl = 1
    for ig = 2:Ng
        if G2_sorted[ig] > (G2_sorted[ig-1] + eps8)
            ngl = ngl + 1
        end
    end

    G2_shells = zeros(ngl)
    idx_g2shells = zeros(Int64,Ng)
    
    G2_shells[1] = G2_sorted[1]
    idx_g2shells[1] = 1

    igl = 1
    for ig = 2:Ng
        if G2_sorted[ig] > (G2_sorted[ig-1] + eps8)
            igl = igl + 1
            G2_shells[igl] = G2_sorted[ig]
        end
        idx_g2shells[ig] = igl
    end

    return G2_shells, idx_g2shells

end




"""
Creates an instance of `GVectorsW`, given the following inputs

- `ecutwfc`: cutoff for wave function (in hartree)

- `gvec`: an instance of `GVectors`

- `kpoints`: an instance of `KPoints`

This function will loop over all kpoints and determine a set of G+k vectors
which has magnitude less than 2*ecutwfc.
"""
function init_gvecw( ecutwfc::Float64, gvec::GVectors, kpoints::KPoints )
    G = gvec.G
    Ng = gvec.Ng
    idx_g2r = gvec.idx_g2r
    #
    kpts = kpoints.k
    Nkpt = kpoints.Nkpt
    #
    Gk2 = zeros(Float64,Ng)
    Gk = zeros(Float64,3)
    idx_gw2g = Array{Array{Int64,1},1}(undef,Nkpt)
    idx_gw2r = Array{Array{Int64,1},1}(undef,Nkpt)
    Ngw = Array{Int64,1}(undef,Nkpt)
    #
    for ik = 1:Nkpt
        for ig = 1:Ng
            Gk[1] = G[1,ig] + kpts[1,ik]
            Gk[2] = G[2,ig] + kpts[2,ik]
            Gk[3] = G[3,ig] + kpts[3,ik]
            Gk2[ig] = Gk[1]^2 + Gk[2]^2 + Gk[3]^2
        end
        idx_gw2g[ik] = findall( 0.5*Gk2 .<= ecutwfc )
        idx_gw2r[ik] = idx_g2r[idx_gw2g[ik]]
        Ngw[ik] = length(idx_gw2g[ik])
    end
    
    Ngwx = maximum(Ngw)

    return GVectorsW( Ngwx, Ngw, idx_gw2g, idx_gw2r, kpoints )

end

"""
Creates uniform real-space grid points for a given sampling points `Ns`
and `LatVecs`

NOTE: This is no longer used in the PWGrid constructor. It is kept here for
in case we need an explicit real space grid.
"""
function init_grid_R( Ns, LatVecs )
    #
    Npoints = prod(Ns)
    #
    R = Array{Float64}(undef,3,Npoints)
    ip = 0
    for k in 0:Ns[3]-1
    for j in 0:Ns[2]-1
    for i in 0:Ns[1]-1
        ip = ip + 1
        R[1,ip] = LatVecs[1,1]*i/Ns[1] + LatVecs[1,2]*j/Ns[2] + LatVecs[1,3]*k/Ns[3]
        R[2,ip] = LatVecs[2,1]*i/Ns[1] + LatVecs[2,2]*j/Ns[2] + LatVecs[2,3]*k/Ns[3]
        R[3,ip] = LatVecs[3,1]*i/Ns[1] + LatVecs[3,2]*j/Ns[2] + LatVecs[3,3]*k/Ns[3]
    end
    end
    end
    #
    return R
end

# ------------------------------------------------------------
# Some differential operators
#-------------------------------------------------------------

function op_nabla( pw::PWGrid, Rhoe::AbstractVector{Float64} )
    G = pw.gvec.G
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    Npoints = prod(pw.Ns)

    RhoeG = R_to_G(pw,Rhoe)[idx_g2r]

    ∇RhoeG_full = zeros(ComplexF64,3,Npoints)
    ∇Rhoe = zeros(Float64,3,Npoints)
    
    for ig = 1:Ng
        ip = idx_g2r[ig]
        ∇RhoeG_full[1,ip] = im*G[1,ig]*RhoeG[ig]
        ∇RhoeG_full[2,ip] = im*G[2,ig]*RhoeG[ig]
        ∇RhoeG_full[3,ip] = im*G[3,ig]*RhoeG[ig]
    end

    ∇Rhoe[1,:] = real(G_to_R(pw,∇RhoeG_full[1,:]))
    ∇Rhoe[2,:] = real(G_to_R(pw,∇RhoeG_full[2,:]))
    ∇Rhoe[3,:] = real(G_to_R(pw,∇RhoeG_full[3,:]))
    return ∇Rhoe

end


function op_nabla_dot( pw::PWGrid, h::Array{Float64,2} )
    G = pw.gvec.G
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    Npoints = prod(pw.Ns)

    hG = zeros(ComplexF64,3,Ng)
    hG[1,:] = R_to_G( pw, h[1,:] )[idx_g2r]
    hG[2,:] = R_to_G( pw, h[2,:] )[idx_g2r]
    hG[3,:] = R_to_G( pw, h[3,:] )[idx_g2r]

    divhG_full = zeros(ComplexF64,Npoints)
    
    for ig = 1:Ng
        ip = idx_g2r[ig]
        divhG_full[ip] = im*( G[1,ig]*hG[1,ig] + G[2,ig]*hG[2,ig] + G[3,ig]*hG[3,ig] )
    end

    divh = real( G_to_R( pw, divhG_full ) )
    return divh

end



function op_nabla_dot(
    pw::PWGrid,
    hGx::Array{ComplexF64,3},
    hGy::Array{ComplexF64,3},
    hGz::Array{ComplexF64,3}
)
    G = pw.gvec.G
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    Npoints = prod(pw.Ns)

    R_to_G!(pw, hGx)
    R_to_G!(pw, hGy)
    R_to_G!(pw, hGz)    

    divhG_full = zeros(ComplexF64, pw.Ns)
    
    for ig = 1:Ng
        ip = idx_g2r[ig]
        divhG_full[ip] = im*( G[1,ig]*hGx[ip] + G[2,ig]*hGy[ip] + G[3,ig]*hGz[ip] )
    end
    G_to_R!(pw, divhG_full)
    return real(divhG_full)
end

include("PWGrid_io.jl")
