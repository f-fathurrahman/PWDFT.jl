function calc_E_NN_simple( LatVecs::Array{Float64,2},
                           atoms::Atoms, Zv::Array{Float64,1};
                           TOL=1e-7, verbose=false )

    TOL = 1e-7

    E_NN_old = 0.0
    E_NN = 0.0
    for i = 1:20
        ecutwfc = 15.0 + (i-1)*5.0
        pw = PWGrid( ecutwfc, LatVecs )
        E_NN = calc_E_NN_simple( pw, atoms, Zv, verbose=verbose )
        diff_E_NN = abs( E_NN - E_NN_old )
        if verbose
            @printf("%18.10f %18.10f %18.10e\n", ecutwfc, E_NN, diff_E_NN)
        end
        if diff_E_NN < TOL
            @printf("Convergence achived in Ewald energy calculation\n")
            break
        end
        E_NN_old = E_NN
    end

    return E_NN

end

#
# Implementation of a simple method to calculate Ewald energy based on Arias' notes
#
function calc_E_NN_simple( pw::PWGrid, atoms::Atoms,
                           Zv::Array{Float64,1}; sigma=nothing, verbose=false )
    #
    Ω  = pw.Ω
    r  = pw.r
    Ns = pw.Ns
    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    Npoints = prod(Ns)
    LatVecs = pw.LatVecs
    Xpos = atoms.positions
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

    Sf = calc_strfact( atoms, pw )

    #
    # Generate array of distances
    dr = gen_dr_center( r, LatVecs )
    #
    # Generate charge density
    #
    if sigma==nothing
      sigma = 0.25*ones(Nspecies)
    end
    #
    g1  = zeros( Float64, Npoints )
    rho_is = zeros( Float64, Npoints, Nspecies )
    Rho = zeros(Float64, Npoints)
    #

    if verbose
        println("\nCalculating E_NN")
    end

    for isp = 1:Nspecies
        c1 = 2*sigma[isp]^2
        cc1 = sqrt(2*pi*sigma[isp]^2)^3
        for ip=1:Npoints
            g1[ip] = exp(-dr[ip]^2/c1)/cc1
        end
        #
        g1 = Zv[isp] * g1[:]
        #
        ctmp = R_to_G( pw, g1 )
        for ig = 1:Ng
            ip = idx_g2r[ig]
            ctmp[ip] = ctmp[ip]*Sf[ig,isp]
        end
        rho_is[:,isp] = real( G_to_R(pw, ctmp) )
        intrho = sum(rho_is[:,isp])*Ω/Npoints
        if verbose
            @printf("Species %d, intrho: %18.10f\n", isp, intrho)
        end
        Rho[:] = Rho[:] + rho_is[:,isp]
    end
    #
    if verbose
        intrho = sum(Rho)*Ω/Npoints
        @printf("Integrated total Gaussian rho: %18.10f\n", intrho)
    end
    #
    # Solve Poisson equation and calculate Hartree energy
    ctmp = 4.0*pi*R_to_G( pw, Rho )
    ctmp[1] = 0.0
    for ig = 2:Ng
        ip = idx_g2r[ig]
        ctmp[ip] = ctmp[ip]/G2[ig]
    end
    phi = real( G_to_R( pw, ctmp ) )
    Ehartree = 0.5*dot( phi, Rho ) * Ω/Npoints
    #
    Eself = 0.0
    Natoms = size(Xpos)[2]
    for ia = 1:Natoms
        isp = atm2species[ia]
        Eself = Eself + Zv[isp]^2/(2*sqrt(pi))*(1.0/sigma[isp])
    end
    E_nn = Ehartree - Eself

    if verbose
        @printf("Ehartree, Eself, E_nn = %18.10f %18.10f %18.10f\n", Ehartree, Eself, E_nn)
    end

    return E_nn
end

# distance to center of the simulation cell
function gen_dr_center( r, LatVecs )
    #
    cx = 0.5*sum(LatVecs[:,1])
    cy = 0.5*sum(LatVecs[:,2])
    cz = 0.5*sum(LatVecs[:,3])
    center = [cx,cy,cz]
    #
    Npoints = size(r)[2]
    dr = Array{Float64}(Npoints)
    for ip=1:Npoints
        dr[ip] = norm( r[:,ip] - center )
    end
    return dr
end
