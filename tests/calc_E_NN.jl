#
# Implementation of a simple method to calculate Ewald energy based on Arias' notes
#
function calc_E_NN( pw::PWGrid, Sf, Xpos, Nspecies::Int, atm2species,
                     Zv::Array{Float64}; sigma=nothing, verbose=false )
    #
    Ω  = pw.Ω
    r  = pw.r
    Ns = pw.Ns
    G2 = pw.gvec.G2
    Npoints = prod(Ns)
    LatVecs = pw.LatVecs
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
        ctmp = R_to_G( Ns, g1 )
        for ip = 1:Npoints
            ctmp[ip] = ctmp[ip]*Sf[ip,isp]
        end
        rho_is[:,isp] = real( G_to_R(Ns, ctmp) )
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
    ctmp = 4.0*pi*R_to_G( Ns, Rho )
    ctmp[1] = 0.0
    for ip = 2:Npoints
        ctmp[ip] = ctmp[ip] / G2[ip]
    end
    phi = real( G_to_R( Ns, ctmp ) )
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
