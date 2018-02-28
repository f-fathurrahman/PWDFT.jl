struct GVectorsW
    Ngwx::Int
    idx_gw2r::Array{Int64}
end

struct GVectors
    Ng::Int64
    G::Array{Float64,2}
    G2::Array{Float64}
end

struct PWGrid
    ecutwfc::Float64
    ecutrho::Float64
    Ns::Tuple{Int64,Int64,Int64}
    LatVecs::Array{Float64,2}
    RecVecs::Array{Float64,2}
    Ω::Float64
    r::Array{Float64,2}
    gvec::GVectors
    gvecw::GVectorsW
    planfw::Base.DFT.FFTW.cFFTWPlan{Complex{Float64},-1,false,3}
    planbw::Base.DFT.ScaledPlan{Complex{Float64},Base.DFT.FFTW.cFFTWPlan{Complex{Float64},1,false,3},Float64}
end


function PWGrid( ecutwfc::Float64, LatVecs::Array{Float64,2} )

    ecutrho = 4.0*ecutwfc
    #
    RecVecs = 2*pi*inv(LatVecs')
    Ω = det(LatVecs)
    #
    LatVecsLen = Array{Float64}(3)
    LatVecsLen[1] = norm(LatVecs[1,:])
    LatVecsLen[2] = norm(LatVecs[2,:])
    LatVecsLen[3] = norm(LatVecs[3,:])

    Ns1 = 2*round( Int, sqrt(ecutrho/2)*LatVecsLen[1]/pi ) + 1
    Ns2 = 2*round( Int, sqrt(ecutrho/2)*LatVecsLen[2]/pi ) + 1
    Ns3 = 2*round( Int, sqrt(ecutrho/2)*LatVecsLen[3]/pi ) + 1

    # Use even sampling numbers
    Ns1 = Ns1 % 2 == 1 ? Ns1 + 1 : Ns1
    Ns2 = Ns2 % 2 == 1 ? Ns2 + 1 : Ns2
    Ns3 = Ns3 % 2 == 1 ? Ns3 + 1 : Ns3

    Ns = (Ns1,Ns2,Ns3)

    Npoints = prod(Ns)
    r = init_grid_R( Ns, LatVecs )

    gvec = init_grid_G( Ns, RecVecs )
    gvecw = init_gvecw( ecutwfc, gvec.G2 )

    planfw = plan_fft( zeros(Ns) )
    planbw = plan_ifft( zeros(Ns) )

    return PWGrid( ecutwfc, ecutrho, Ns, LatVecs, RecVecs, Ω, r, gvec, gvecw,
                   planfw, planbw )
end

function mm_to_nn(mm::Int64,S::Int64)
    if mm > S/2
        return mm - S
    else
        return mm
    end
end


function init_grid_G( Ns, RecVecs )

    Ng = prod(Ns)

    G  = Array{Float64}(3,Ng)
    G2 = Array{Float64}(Ng)

    ig = 0
    for k in 0:Ns[3]-1
    for j in 0:Ns[2]-1
    for i in 0:Ns[1]-1
        ig = ig + 1
        gi = mm_to_nn( i, Ns[1] )
        gj = mm_to_nn( j, Ns[2] )
        gk = mm_to_nn( k, Ns[3] )
        G[1,ig] = RecVecs[1,1]*gi + RecVecs[2,1]*gj + RecVecs[3,1]*gk
        G[2,ig] = RecVecs[1,2]*gi + RecVecs[2,2]*gj + RecVecs[3,2]*gk
        G[3,ig] = RecVecs[1,3]*gi + RecVecs[2,3]*gj + RecVecs[3,3]*gk
        G2[ig] = G[1,ig]^2 + G[2,ig]^2 + G[3,ig]^2
    end
    end
    end

    return GVectors( Ng, G, G2 )
end

function init_gvecw( ecutwfc, G2 )
    idx_gw2r = findn( 0.5*G2 .< ecutwfc )
    Ngwx = length(idx_gw2r)
    return GVectorsW( Ngwx, idx_gw2r )
end


function init_grid_R( Ns, LatVecs )
    #
    Npoints = prod(Ns)
    #
    R = Array{Float64}(3,Npoints)
    ip = 0
    for k in 0:Ns[3]-1
    for j in 0:Ns[2]-1
    for i in 0:Ns[1]-1
        ip = ip + 1
        R[1,ip] = LatVecs[1,1]*i/Ns[1] + LatVecs[2,1]*j/Ns[2] + LatVecs[3,1]*k/Ns[3]
        R[2,ip] = LatVecs[1,2]*i/Ns[1] + LatVecs[2,2]*j/Ns[2] + LatVecs[3,2]*k/Ns[3]
        R[3,ip] = LatVecs[1,3]*i/Ns[1] + LatVecs[2,3]*j/Ns[2] + LatVecs[3,3]*k/Ns[3]
    end
    end
    end
    #
    return R
end

import Base.println
function println( pw::PWGrid )
    @printf("\nPlane wave grid\n\n")
    @printf("ecutwfc = %10.3f Ha\n", pw.ecutwfc)
    @printf("ecutrho = %10.3f Ha\n", pw.ecutrho)
    LatVecs = pw.LatVecs
    RecVecs = pw.RecVecs
    @printf("Direct lattice vectors:\n")
    for i = 1:3
        @printf("%18.10f %18.10f %18.10f\n", LatVecs[i,1], LatVecs[i,2], LatVecs[i,3])
    end
    @printf("\nReciprocal lattice vectors:\n")
    for i = 1:3
        @printf("%18.10f %18.10f %18.10f\n", RecVecs[i,1], RecVecs[i,2], RecVecs[i,3])
    end
    @printf("\n")
    @printf("Direct lattive volume = %18.10f\n", pw.Ω )
    @printf("Sampling points: [%5d,%5d,%5d]\n", pw.Ns[1], pw.Ns[2], pw.Ns[3])
    @printf("Ngwx = %10d\n", pw.gvecw.Ngwx)
end
