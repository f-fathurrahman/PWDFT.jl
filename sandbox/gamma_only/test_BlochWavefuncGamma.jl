using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

include("PWGridGamma.jl")
include("wrappers_fft_gamma.jl")
include("ortho_GS_gamma.jl")
include("PsPotNLGamma.jl")
include("HamiltonianGamma.jl")
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

    println("Before normalization:")
    println("psiG[1] = ", psiG[1])
    println("dot psiG = ", dot(psiG,psiG))

    integPsi = dot(psiG,psiG)
    psiG = psiG/sqrt(integPsi)

    println("After normalization:")
    println("psiG[1] = ", psiG[1])
    println("dot psiG = ", dot(psiG,psiG), " (should be close to one)")

    psiG_gamma = zeros(ComplexF64, pw_gamma.gvec.Ng)
    psiG_gamma[1] = ctmp[1]
    for ig = 2:pw_gamma.gvec.Ng
        ip = pw_gamma.gvec.idx_g2r[ig]
        psiG_gamma[ig] = ctmp[ip]
    end
    println()


    integPsi = 2*dot(psiG_gamma,psiG_gamma) - conj(psiG_gamma[1])*psiG_gamma[1]

    println("Before normalization:")
    println("psiG_gamma[1] = ", psiG_gamma[1])
    println("dot psiG_gamma = ", integPsi)
    
    # This is how the normalization is done for Gamma-only trick
    psiG_gamma = psiG_gamma/sqrt(integPsi)

    integPsi = 2*dot(psiG_gamma,psiG_gamma) - conj(psiG_gamma[1])*psiG_gamma[1]

    println("After normalization:")
    println("psiG_gamma[1] = ", psiG_gamma[1])
    println("dot psiG_gamma = ", integPsi, " (should be close to one)")

    Δ = integPsi - dot(psiG,psiG)
    println()
    println("Δ = ", Δ, " (should be very small)")
end
#test_01()



# Generate normalized psi
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

    # Should set DC component of psi1 and psi2 to zero

    for ist in 1:Nstates
        
        ctmp = R_to_G(pw, psiR[:,ist])

        ctmp[1] = 0.0 + im*0.0
        
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
        integPsi = 2*dot(psi2[:,ist], psi2[:,ist]) - conj(psi2[1,ist])*psi2[1,ist]
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
    println("dot(psi1,psi1) = ", res1, " (should be close to one)")
    println("dot(psi2,psi2) = ", res2, " (should be close to one)")
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
        integPsi = 2*dot(psi2.data[1][:,ist], psi2.data[1][:,ist]) - psi2.data[1][1,ist]*psi2.data[1][1,ist]
        #println("psi2.data[1] = ", psi2.data[1,ist])
        ss = ss + integPsi
        println("integPsi = ", integPsi)
    end
    println("ss = ", ss)

end
#test_03()



function test_04()
    psis = randn_BlochWavefuncGamma(10,4)

    println()
    println("dot(psis,psis) = ", dot(psis,psis), " (should be close to Nstates)")

    psi1 = psis.data[1]

    println()
    println(2*dot(psi1[:,1],psi1[:,1]), " (should be close to 1)")
    c = dot(psi1[:,1],psi1[:,2])
    println(2*c)
    println(c + conj(c), " (should be close to 0)")

    ortho_check(psis)
end
#test_04()



function test_05()
    psi1 = randn_BlochWavefuncGamma(6,5)
    psi2 = randn_BlochWavefuncGamma(6,5)
    
    #ortho_check(psi1)
    #ortho_check(psi2)

    psi3 = psi1 + psi2
    ortho_gram_schmidt!(psi3)
    ortho_check(psi3)

    ortho_check(psi1 - psi2)

    println()
    println(dot(psi1,psi1), " (should be close to Nstates)")
    println(dot(psi2,psi2), " (should be close to Nstates)")
    println(dot(psi1,psi2))
    println(dot(psi1,psi3))

end
#test_05()