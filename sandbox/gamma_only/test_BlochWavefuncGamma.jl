using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

include("PWGridGamma.jl")
include("ortho_GS_gamma.jl")
include("BlochWavefuncGamma.jl")

function test_01()

    Random.seed!(1234)

    LatVecs = gen_lattice_sc(6.0)

    pw = PWGrid(5.0, LatVecs)

    pw_gamma = PWGridGamma(5.0, LatVecs)

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
    println("diff = ", Δ - psiG_gamma[1]^2) # should be the same as Δ
end

#test_01()

# Generate normalize psi
# No orthogonalization is imposed between psi
function gen_two_wavefunc( pw::PWGrid, pw_gamma::PWGridGamma; Nstates=1 )
    
    Npoints = prod(pw.Ns)
    dVol = pw.CellVolume/Npoints

    psiR = randn(Float64,Npoints,Nstates)
    integPsiR = dot(psiR, psiR)*dVol
    psiR = psiR/sqrt(integPsiR)
    
    integPsiR = sum( psiR .* psiR )*dVol
    
    ctmp = zeros(ComplexF64,Npoints)
    psi1 = zeros(ComplexF64, pw.gvec.Ng, Nstates)
    psi2 = zeros(ComplexF64, pw_gamma.gvec.Ng, Nstates)

    for ist in 1:Nstates
        
        ctmp = R_to_G(pw, psiR[:,ist])
        
        # psi for pw
        for ig = 1:pw.gvec.Ng
            ip = pw.gvec.idx_g2r[ig]
            psi1[ig,ist] = ctmp[ip]
        end
        integPsi = dot( psi1[:,ist], psi1[:,ist] )
        psi1[:,ist] = psi1[:,ist]/sqrt(integPsi)

        # psi for pw_gamma
        for ig = 1:pw_gamma.gvec.Ng
            ip = pw_gamma.gvec.idx_g2r[ig]
            psi2[ig,ist] = ctmp[ip]
        end
        integPsi = 2*dot(psi2[:,ist], psi2[:,ist]) - psi2[1,ist]*psi2[1,ist]
        psi2[:,ist] = psi2[:,ist]/sqrt(integPsi)

    end

    return psi1, BlochWavefuncGamma(psi2)

end

function test_02()
    Random.seed!(1234)

    ecutwfc = 5.0
    LatVecs = gen_lattice_sc(6.0)
    
    pw = PWGrid(ecutwfc, LatVecs)
    pw_gamma = PWGridGamma(ecutwfc, LatVecs)

    psi1, psi2 = gen_two_wavefunc(pw, pw_gamma)

    println()

    println("size psi1 = ", Base.summarysize(psi1))
    println("size psi2 = ", Base.summarysize(psi2))

    res1 = dot(psi1,psi1)
    res2 = dot(psi2,psi2)
    println("dot(psi1,psi1) = ", res1)
    println("dot(psi2,psi2) = ", res2)
    println("diff = ", res1 - res2)
end

#test_02()

function test_03()
    Random.seed!(1234)

    ecutwfc = 5.0
    LatVecs = gen_lattice_sc(6.0)
    
    pw = PWGrid(ecutwfc, LatVecs)
    pw_gamma = PWGridGamma(ecutwfc, LatVecs)

    psi1, psi2 = gen_two_wavefunc(pw, pw_gamma, Nstates=4)

    println()

    println("size psi1 = ", Base.summarysize(psi1))
    println("size psi2 = ", Base.summarysize(psi2))

    res1 = dot(psi1,psi1)
    res2 = dot(psi2,psi2)
    println("dot(psi1,psi1) = ", res1)
    println("dot(psi2,psi2) = ", res2)
    println("diff = ", res1 - res2)
    
    ss = 0.0
    for ist in 1:4
        integPsi = 2*dot(psi2.data[:,ist], psi2.data[:,ist]) - psi2.data[1,ist]*psi2.data[1,ist]
        #println("psi2.data[1] = ", psi2.data[1,ist])
        ss = ss + integPsi
        println("integPsi = ", integPsi)
    end
    println("ss = ", ss)

end
#test_03()

function test_04()
    psi = randn_BlochWavefuncGamma(6,4)
    ortho_check(psi)
    #println(psi); println()
    println("dot 1 2:", 2*dot(psi.data[:,1],psi.data[:,2]) )
    ortho_GS_gamma!(psi.data)
    println("dot 1 2:", 2*dot(psi.data[:,1],psi.data[:,2]) )
end
test_04()