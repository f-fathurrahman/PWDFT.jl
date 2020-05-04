using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

include("PWGridGammaOnly.jl")
include("BlochWavefuncGammaOnly.jl")

function test_01()

    Random.seed!(1234)

    LatVecs = gen_lattice_sc(6.0)

    pw = PWGrid(5.0, LatVecs)

    pw_gamma = PWGridGammaOnly(5.0, LatVecs)

    Npoints = prod(pw.Ns)
    dVol = pw.CellVolume/Npoints

    psiR = randn(Float64,Npoints)
    integPsiR = dot(psiR, psiR)*dVol
    psiR = psiR/sqrt(integPsiR)
    
    integPsiR = sum( psiR .* psiR )*dVol
    println("integPsiR = ", integPsiR)
    
    ctmp = R_to_G(pw, psiR)

    for ig in 2:4
        ip = pw_gamma.gvec.idx_g2r[ig]
        ipm = pw_gamma.gvec.idx_g2rm[ig]
        @printf("+G ig=%3d: ", ig); print(ctmp[ip]); println()
        @printf("-G ig=%3d: ", ig); print(ctmp[ipm]); println()
    end

    psiG = zeros(ComplexF64,pw.gvec.Ng)
    for ig = 1:pw.gvec.Ng
        ip = pw.gvec.idx_g2r[ig]
        psiG[ig] = ctmp[ip]
    end
    println()
    integPsi = dot(psiG,psiG)
    psiG = psiG/sqrt(integPsi)
    println("psiG[1] = ", psiG[1])
    println("dot psiG = ", dot(psiG,psiG))
    println("dot psiG = ", dot(psiG[2:end],psiG[2:end]))

    psiG_gamma = zeros(ComplexF64, pw_gamma.gvec.Ng)
    psiG_gamma[1] = ctmp[1]
    for ig = 2:pw_gamma.gvec.Ng
        ip = pw_gamma.gvec.idx_g2r[ig]
        psiG_gamma[ig] = ctmp[ip]
    end
    println()
    # This is how the normalization is done for Gamma-only trick
    integPsi = 2*dot(psiG_gamma,psiG_gamma) - psiG_gamma[1]*psiG_gamma[1]
    psiG_gamma = psiG_gamma/sqrt(integPsi)
    println("psiG_gamma[1] = ", psiG_gamma[1])
    println("dot psiG_gamma = ", dot(psiG_gamma,psiG_gamma))
    println("dot psiG_gamma*2 = ", 2*dot(psiG_gamma[2:end],psiG_gamma[2:end]))

    Δ = 2*dot(psiG_gamma,psiG_gamma) - dot(psiG,psiG)
    println("Δ = ", Δ)
    println("psiG_gamma[1]^2 = ", psiG_gamma[1]^2) # should be the same as Δ

end

#test_01()

function gen_two_wavefunc( pw::PWGrid, pw_gamma::PWGridGammaOnly )
    
    Npoints = prod(pw.Ns)
    dVol = pw.CellVolume/Npoints

    psiR = randn(Float64,Npoints)
    integPsiR = dot(psiR, psiR)*dVol
    psiR = psiR/sqrt(integPsi)
    
    integPsi = sum( psiR .* psiR )*dVol
    println("integPsi = ", integPsi)
    
    ctmp = R_to_G(pw, psiR)

    # psi for pw
    psi1 = zeros(ComplexF64, pw.gvec.Ng)
    for ig = 1:pw.gvec.Ng
        ip = pw.gvec.idx_g2r[ig]
        psi1[ig] = ctmp[ip]
    end

    psi2 = zeros(ComplexF64, pw_gamma.gvec.Ng)
    for ig = 1:pw_gamma.gvec.Ng
        ip = pw_gamma.gvec.idx_g2r[ig]
        psi2[ig] = ctmp[ip]
    end

    return psi1, psi2

end

function test_02()

end