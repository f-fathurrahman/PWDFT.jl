function calc_stress_hartree!(pw, Rhoe, E_hartree, stress_hartree)

    G = pw.gvec.G
    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    Npoints = size(Rhoe, 1)

    fill!(stress_hartree, 0.0)
    ctmp = zeros(ComplexF64, Npoints)
    # Use Nspin=1
    for ip in 1:Npoints
        ctmp[ip] = Rhoe[ip,1]
    end
    # We create again plan_fft
    planfw = plan_fft!(zeros(ComplexF64, pw.Ns)) # using default plan
    ff = reshape(ctmp, pw.Ns)
    planfw*ff
    @views ctmp[:] /= Npoints # rescale

    for ig in 2:Ng
        ip = idx_g2r[ig]
        shart = ctmp[ip] * conj(ctmp[ip]) / G2[ig]
        for l in 1:3, m in 1:l
            stress_hartree[l,m] += 2 * shart * G[l,ig] * G[m,ig] / G2[ig]
        end
    end
    #
    # stress_hartree[:,:] *= 4π # gamma only
    stress_hartree[:,:] *= 2π
    #
    println("E_hartree = ", E_hartree)
    for l in 1:3
        stress_hartree[l,l] -= E_hartree / pw.CellVolume
    end
    #
    for l in 1:3, m in 1:(l-1)
        stress_hartree[m,l] = stress_hartree[l,m]
    end
    #
    println("Before multiplied by -1")
    display(stress_hartree); println()
    stress_hartree[:,:] *= -1
    println("After multiplied by -1")
    display(stress_hartree); println()
    #
    return
end