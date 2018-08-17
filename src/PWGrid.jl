"""
The type for set of G-vectors for describing density and potentials
"""
struct GVectors
    Ng::Int64
    G::Array{Float64,2}
    G2::Array{Float64,1}
    idx_g2r::Array{Int64,1}
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

"""
The type for describing plane wave basis set for a given periodic unit cell.
"""
struct PWGrid
    ecutwfc::Float64
    ecutrho::Float64
    Ns::Tuple{Int64,Int64,Int64}
    LatVecs::Array{Float64,2}
    RecVecs::Array{Float64,2}
    CellVolume::Float64
    r::Array{Float64,2}
    gvec::GVectors
    gvecw::GVectorsW
    planfw::FFTW.cFFTWPlan{Complex{Float64},-1,false,3}
    planbw::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,false,3},Float64}
end


"""
Creates an instance of `PWGrid` given the following inputs:

- `ecutwfc`: cutoff energy for wave function

- `LatVecs`: unit cell lattice vectors.
  Convention used: The three lattice vectors v1, v2, and v3 arranged are
  arranged **by column**.

- `kpoints`: optional, an instance of `KPoints`.
"""
function PWGrid( ecutwfc::Float64, LatVecs::Array{Float64,2}; kpoints=nothing )

    ecutrho = 4.0*ecutwfc
    #
    #RecVecs = 2*pi*inv(LatVecs')
    RecVecs = 2*pi*invTrans_m3x3(LatVecs)

    CellVolume = det(LatVecs)
    #
    LatVecsLen = Array{Float64}(undef,3)
    LatVecsLen[1] = norm(LatVecs[:,1])
    LatVecsLen[2] = norm(LatVecs[:,2])
    LatVecsLen[3] = norm(LatVecs[:,3])

    Ns1 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[1]/pi ) + 1
    Ns2 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[2]/pi ) + 1
    Ns3 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[3]/pi ) + 1

    Ns1 = good_fft_order(Ns1)
    Ns2 = good_fft_order(Ns2)
    Ns3 = good_fft_order(Ns3)

    Ns = (Ns1,Ns2,Ns3)

    Npoints = prod(Ns)
    r = init_grid_R( Ns, LatVecs )
    
    gvec = init_gvec( Ns, RecVecs, ecutrho )

    if kpoints==nothing
        kpoints = KPoints( 1, zeros(3,1), [1.0], RecVecs )
    end

    gvecw = init_gvecw( ecutwfc, gvec, kpoints )

    planfw = plan_fft( zeros(Ns) )
    planbw = plan_ifft( zeros(Ns) )

    return PWGrid( ecutwfc, ecutrho, Ns, LatVecs, RecVecs, CellVolume, r, gvec, gvecw,
                   planfw, planbw )
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
        if 0.5*G2 < ecutrho
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

    ig = 0
    ip = 0
    for k in 0:Ns[3]-1
    for j in 0:Ns[2]-1
    for i in 0:Ns[1]-1
        ip = ip + 1
        gi = mm_to_nn( i, Ns[1] )
        gj = mm_to_nn( j, Ns[2] )
        gk = mm_to_nn( k, Ns[3] )
        G_temp[1] = RecVecs[1,1]*gi + RecVecs[1,2]*gj + RecVecs[1,3]*gk
        G_temp[2] = RecVecs[2,1]*gi + RecVecs[2,2]*gj + RecVecs[2,3]*gk
        G_temp[3] = RecVecs[3,1]*gi + RecVecs[3,2]*gj + RecVecs[3,3]*gk
        G2_temp = G_temp[1]^2 + G_temp[2]^2 + G_temp[3]^2
        if 0.5*G2_temp <= ecutrho
            ig = ig + 1
            G[:,ig] = G_temp[:]
            G2[ig] = G2_temp
            idx_g2r[ig] = ip
        end
    end
    end
    end

    return GVectors( Ng, G, G2, idx_g2r )
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
            Gk[:] = G[:,ig] .+ kpts[:,ik]
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


# Overloaded println

import Base: println

"""
Display some information about `pw::PWGrid`. This function calls
`println(gvec::GVectors)` and `println(gvecw::GVectorsW)`.
"""
function println( pw::PWGrid; header=true )
    if header
        @printf("\n")
        @printf("                                     ------\n")
        @printf("                                     PWGrid\n")
        @printf("                                     ------\n")
        @printf("\n")
    end
    LatVecs = pw.LatVecs
    RecVecs = pw.RecVecs
    @printf("Direct lattice vectors:\n")
    @printf("\n")
    for i = 1:3
        @printf("%18.10f %18.10f %18.10f\n", LatVecs[i,1], LatVecs[i,2], LatVecs[i,3])
    end
    @printf("\n")
    @printf("Reciprocal lattice vectors:\n")
    @printf("\n")
    for i = 1:3
        @printf("%18.10f %18.10f %18.10f\n", RecVecs[i,1], RecVecs[i,2], RecVecs[i,3])
    end
    @printf("\n")
    @printf("Direct lattive volume = %18.10f bohr^3\n", pw.CellVolume )
    @printf("ecutwfc               = %18.10f Ha\n", pw.ecutwfc)
    @printf("ecutrho               = %18.10f Ha\n", pw.ecutrho)    
    @printf("Sampling points       = (%5d,%5d,%5d)\n", pw.Ns[1], pw.Ns[2], pw.Ns[3])
    #
    println( pw.gvec )
    println( pw.gvec, pw.gvecw )
end

"""
Display some information about `gvec::GVectors`.
"""
function println( gvec::GVectors )
    Ng = gvec.Ng
    G = gvec.G
    G2 = gvec.G2
    
    @printf("\n")
    @printf("                                    --------\n")
    @printf("                                    GVectors\n")
    @printf("                                    --------\n")
    @printf("\n")
    @printf("Ng = %12d\n", Ng)
    @printf("\n")
    for ig = 1:3
        @printf("%8d [%18.10f,%18.10f,%18.10f] : %18.10f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])        
    end
    @printf(" ....... \n")
    for ig = Ng-3:Ng
        @printf("%8d [%18.10f.%18.10f,%18.10f] : %18.10f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])
    end
    @printf("\n")
    @printf("Max G2 = %18.10f\n", maximum(G2))
end

"""
Display some information about `gvecw::GVectorsW`.
"""
function println( gvec::GVectors, gvecw::GVectorsW )
    G = gvec.G
    G2 = gvec.G2

    Ngwx = gvecw.Ngwx
    Ngw = gvecw.Ngw
    k = gvecw.kpoints.k
    Nkpt = gvecw.kpoints.Nkpt

    @printf("\n")
    @printf("                                    ---------\n")
    @printf("                                    GVectorsW\n")
    @printf("                                    ---------\n")
    @printf("\n")
    @printf("Ngwx = %12d\n", Ngwx)
        
    for ik = 1:Nkpt
        idx_gw2g = gvecw.idx_gw2g[ik]
        Gw = zeros(3,Ngw[ik])
        Gw2 = zeros(Ngw[ik])
        for igk = 1:Ngw[ik]
            ig = idx_gw2g[igk]
            Gw[:,igk] = G[:,ig] + k[:,ik]
            Gw2[igk] = Gw[1,igk]^2 + Gw[2,igk]^2 + Gw[3,igk]^2
        end
        @printf("Max Gw2[%3d] = %18.10f\n", ik, maximum(Gw2))    
    end
end

# ------------------------------------------------------------
# Probably should be moved to PWGrid
#-------------------------------------------------------------

function op_nabla( pw::PWGrid, Rhoe::Array{Float64,1} )
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
