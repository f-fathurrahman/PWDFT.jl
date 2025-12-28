function atomic_rho_g(
    Ham::Hamiltonian{PsPot_UPF};
    starting_magn::Union{Nothing,Vector{Float64}} = nothing,
    angle1::Union{Nothing,Vector{Float64}} = nothing,
    angle2::Union{Nothing,Vector{Float64}} = nothing
)

    atoms = Ham.atoms
    pw = Ham.pw
    pspots = Ham.pspots
    Nspecies = atoms.Nspecies
    Ng = pw.gvec.Ng
    idx_g2shells = pw.gvec.idx_g2shells
    G2_shells = pw.gvec.G2_shells
    Npoints = prod(pw.Ns)
    idx_g2r = pw.gvec.idx_g2r
    CellVolume = pw.CellVolume
    Nelectrons = Ham.electrons.Nelectrons
    Nspin = Ham.electrons.Nspin_comp
    Ngl = length(Ham.pw.gvec.G2_shells)

    eps8 = 1e-8
    strf = calc_strfact( atoms, pw )

    # Determine maximum radial points for every atomic species
    NpointsMax = 0
    for isp in 1:Nspecies
        if NpointsMax < pspots[isp].Nr
            NpointsMax = pspots[isp].Nr
        end
    end
    #println("NpointsMax = ", NpointsMax)

    aux = zeros(Float64, NpointsMax)
    rhocgnt = zeros(Float64, Ngl)
    rhocg = zeros(ComplexF64, Npoints, Nspin)

    # No starting magnetization is give, set them to a default value
    if (Nspin == 2) && isnothing(starting_magn)
        starting_magn = 0.1*ones(Nspecies)
        @info "Using default starting magnetization = $(starting_magn)"
    end

    if (Nspin == 4) && isnothing(angle1) && isnothing(angle2)
        starting_magn = 0.1*ones(Nspecies)
        angle1 = 0.4*pi*ones(Nspecies)
        angle2 = 0.5*pi*ones(Nspecies)
        @info "Using default angle1 = $(angle1)"
        @info "Using default angle2 = $(angle2)"
    end

    angular = zeros(Float64, 3)

    for isp in 1:Nspecies

        psp = pspots[isp]

        #println("sum psp.rhoatom = ", sum(psp.rhoatom))

        # G == 0 term
        for ir in 1:psp.Nr
            aux[ir] = psp.rhoatom[ir]
        end
        rhocgnt[1] = PWDFT.integ_simpson( psp.Nr, aux, psp.rab )

        # G != 0 terms
        for igl in 2:Ngl
            gx = sqrt(G2_shells[igl])
            for ir in 1:psp.Nr
                if psp.r[ir] < eps8
                   aux[ir] = psp.rhoatom[ir]
                else
                   aux[ir] = psp.rhoatom[ir]*sin(gx*psp.r[ir])/(psp.r[ir]*gx)
                end
            end
            rhocgnt[igl] = PWDFT.integ_simpson( psp.Nr, aux, psp.rab )
        end

        #println("sum rhocgnt = ", sum(rhocgnt))

        for ig in 1:Ng
            ip = idx_g2r[ig]
            igl = idx_g2shells[ig]
            rhocg[ip,1] += strf[ig,isp]*rhocgnt[igl]/CellVolume
        end

        # Collinear magnetism
        if Nspin == 2
            for ig in 1:Ng
                ip = idx_g2r[ig]
                igl = idx_g2shells[ig]
                rhocg[ip,2] += starting_magn[isp]*strf[ig,isp]*rhocgnt[igl]/CellVolume
            end
        end
        
        # ffr: spherical coord?
        if Nspin == 4
            angular[1] = sin(angle1[isp])*cos(angle2[isp])
            angular[2] = sin(angle1[isp])*sin(angle2[isp])
            angular[3] = cos(angle1[isp])
            # For spin-polarized case
            for ispin in 2:Nspin
                for ig in 1:Ng
                    ip = idx_g2r[ig]
                    igl = idx_g2shells[ig]
                    # rhoscale is set to 1
                    rho1 = angular[ispin-1] * strf[ig,isp] * rhocgnt[igl] / CellVolume
                    rhocg[ip,ispin] += starting_magn[isp] * rho1
                end
            end
        end # if

    end # loop over species

    println("Some rhocg: ")
    for ip in 1:5
        @printf("%3d %18.10f %18.10f\n", ip, real(rhocg[ip,1]), imag(rhocg[ip,1]))
    end
    if Nspin == 4
        println("Some magn cg: ")
        for ip in 1:5
            println("$ip $(rhocg[ip,2]) $(rhocg[ip,3]) $(rhocg[ip,4])" )
        end
    end

    println("sum abs rhocg = ", sum(abs.(rhocg)))

    charge = rhocg[1,1]*CellVolume
    SMALL_CHARGE = 1e-10
    if abs(charge) < SMALL_CHARGE
        @info "Small total charge is detected in atomic_rho_g: probably no rhoatom"
        @info "Initializing to uniform charge density"
        fill!(rhocg, 0.0 + im*0.0)
        rhocg[1,1] = Nelectrons/CellVolume
        # Is this OK?
        if Nspin == 2
            rhocg[1,2] = sum(starting_magn)
        end
    else
        print("atomic_rho_g: Initial charge = ", charge)
        # Renormalize
        println(" Renormalized to ", Nelectrons)
        rhocg .*= Nelectrons/charge
    end

    Rhoe_tot = real(G_to_R(pw,rhocg[:,1]))*Npoints

    # Convert to Rhoe_up and Rhoe_dn    
    # Rhoe_tot = Rhoe_up + Rhoe_dn
    # magn = Rhoe_up - Rhoe_dn
    # 2*Rhoe_up = Rhoe_tot + magn
    # 2*Rhoe_dn = Rhoe_tot - magn
    Rhoe = zeros(Float64, Npoints, Nspin) # perspin
    if Nspin == 2
        magn = real(G_to_R(pw,rhocg[:,2]))*Npoints
        Rhoe[:,1] = 0.5*(Rhoe_tot + magn)
        Rhoe[:,2] = 0.5*(Rhoe_tot - magn)
    elseif Nspin == 4
        for i in 1:4
            Rhoe[:,i] = real(G_to_R(pw,rhocg[:,i]))*Npoints
        end
    else
        Rhoe[:,1] = Rhoe_tot
    end

    println("atomic_rho_g: integ rhoe = ", sum(Rhoe)*CellVolume/Npoints)

    return Rhoe, rhocg

end

function _check_negative_rhoe(Rhoe, CellVolume)
    Nspin = size(Rhoe, 2)
    Npoints = size(Rhoe, 1)
    for ispin in 1:Nspin
        rhoneg = 0.0
        for ip in 1:Npoints
            rhoneg = rhoneg + min(0.0, Rhoe[ip,ispin])
        end
        rhoneg = rhoneg * CellVolume/Npoints
        @printf("Negative rho for ispin %d: %15.10e\n", ispin, rhoneg)
    end
    return
end