using Printf
using LinearAlgebra
using FFTW

using PWDFT

include("PWGridGamma.jl")
include("wrappers_fft_gamma.jl")
include("Poisson_solve_gamma.jl")

function gen_dr( r, center )
    Npoints = size(r)[2]
    dr = Array{Float64}(undef,Npoints)
    #
    for ip=1:Npoints
        dx2 = ( r[1,ip] - center[1] )^2
        dy2 = ( r[2,ip] - center[2] )^2
        dz2 = ( r[3,ip] - center[3] )^2
        dr[ip] = sqrt( dx2 + dy2 + dz2 )
    end
    return dr
end

function gen_rho( dr, σ1, σ2 )
    Npoints = size(dr)[1]
    rho = zeros(Float64,Npoints,1) # Add second dim to rho
    c1 = 2*σ1^2
    c2 = 2*σ2^2
    cc1 = sqrt(2*pi*σ1^2)^3
    cc2 = sqrt(2*pi*σ2^2)^3
    for ip=1:Npoints
        g1 = exp(-dr[ip]^2/c1)/cc1
        g2 = exp(-dr[ip]^2/c2)/cc2
        rho[ip] = g2 - g1
    end
    return rho
end


function test_main( ecutwfc_Ry::Float64 )
    #
    LatVecs = gen_lattice_sc(16.0)
    #
    pw = PWGridGamma( ecutwfc_Ry*0.5, LatVecs )
    #
    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    r = PWDFT.init_grid_R( pw.Ns, LatVecs )
    Ns = pw.Ns
    #
    # Generate array of distances
    #
    center = sum(LatVecs,dims=2)/2
    dr = gen_dr( r, center )
    #
    # Generate charge density
    #
    σ1 = 0.75
    σ2 = 0.50
    rho = gen_rho( dr, σ1, σ2 )
    #
    # Solve Poisson equation and calculate Hartree energy
    #
    phiG = Poisson_solve( pw, rho )
    phi = real( G_to_R(pw, phiG) )
    Ehartree = 0.5*dot( phi, rho ) * CellVolume/Npoints
    #
    Uanal = ( (1/σ1 + 1/σ2)/2 - sqrt(2) / sqrt( σ1^2 + σ2^2 ) ) / sqrt(pi)
    @printf("Num, ana, diff = %18.10f %18.10f %18.10e\n", Ehartree, Uanal, abs(Ehartree-Uanal))
    #
    # Calculate Ehartree using reciprocal space
    #
    rhoG = R_to_G( pw, rho[:] ) / Npoints  # XXX normalize by Npoints
    EhartreeG = 0.0
    Ng = pw.gvec.Ng
    G = pw.gvec.G
    G2 = pw.gvec.G2
    idx_g2r = pw.gvec.idx_g2r
    idx_g2rm = pw.gvec.idx_g2rm
    for ig = 2:Ng
        ip = idx_g2r[ig]
        EhartreeG = EhartreeG + 0.5*real(phiG[ip]*conj(rhoG[ip]))*CellVolume/Npoints
        ipm = idx_g2rm[ig]
        EhartreeG = EhartreeG + 0.5*real(phiG[ipm]*conj(rhoG[ipm]))*CellVolume/Npoints
    end
    println("Using G-space summation")
    @printf("Num, ana, diff = %18.10f %18.10f %18.10e\n", EhartreeG, Uanal, abs(EhartreeG-Uanal))
end

#test_main(10.0)
#test_main(15.0)
test_main(30.0)
