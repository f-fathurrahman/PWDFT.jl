using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

include("PWGridGammaOnly.jl")

function main()
    
    Random.seed!(1234)

    LatVecs = gen_lattice_sc(6.0)

    pw = PWGrid(5.0, LatVecs)
    Ns = pw.Ns
    RecVecs = pw.RecVecs
    ecutrho = pw.ecutrho

    gvec_gamma = GVectorsGammaOnly(Ns, RecVecs, ecutrho)

    G = gvec_gamma.G
    G2 = gvec_gamma.G2
    for ig = 1:4
        @printf("%4d [%10.5f,%10.5f,%10.5f] %10.5f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])
    end

    # Ordinary G
    println()
    for ig in 1:4
        @printf("%4d [%10.5f,%10.5f,%10.5f] %10.5f\n", ig,
            pw.gvec.G[1,ig], pw.gvec.G[2,ig], pw.gvec.G[3,ig], pw.gvec.G2[ig])
    end

    Npoints = prod(Ns)
    dVol = pw.CellVolume/Npoints
    psi = randn(Float64,Npoints) # TODO: normalize in R-space
    integPsi = sum(psi.*psi)*dVol
    psi = psi/sqrt(integPsi)
    integPsi = sum(psi.*psi)*dVol
    println("integPsi = ", integPsi)
    ctmp = R_to_G(pw, psi)

    #for ig in 2:4
    #    ip = gvec_gamma.idx_g2r[ig]
    #    ipm = gvec_gamma.idx_g2rm[ig]
    #    @printf("+G ig=%3d: ", ig); print(ctmp[ip]); println()
    #    @printf("-G ig=%3d: ", ig); print(ctmp[ipm]); println()
    #end

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

    psiG_gamma = zeros(ComplexF64,gvec_gamma.Ng)
    psiG_gamma[1] = ctmp[1]
    for ig = 2:gvec_gamma.Ng
        ip = gvec_gamma.idx_g2r[ig]
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
    println("psiG_gamma[1]^2 = ", psiG_gamma[1]^2)
end

main()