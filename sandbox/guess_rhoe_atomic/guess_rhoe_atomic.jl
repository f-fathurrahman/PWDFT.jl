function guess_rhoe_atomic( Ham::Hamiltonian; starting_magnetization=nothing )

    atoms = Ham.atoms
    pw = Ham.pw
    pspots = Ham.pspots

    Nspecies = atoms.Nspecies

    xmin = -7.0
    amesh = 0.0125
    rmax = 100.0

    Zatoms = get_Zatoms( atoms )

    r = Array{Array{Float64,1},1}(undef,Nspecies)
    rab = Array{Array{Float64,1},1}(undef,Nspecies)
    rho_at = Array{Array{Float64,1},1}(undef,Nspecies)

    decay_lengths = PWDFT.get_atmlength( atoms, pspots )

    NptsRadial = zeros(Int64,Nspecies)

    eps8 = 1e-8

    strf = calc_strfact( atoms, pw )

    for isp = 1:Nspecies

        two_l2 = 2.0*decay_lengths[isp]^2
        znorml = 4.0*pi*pspots[isp].zval/(pi*two_l2)^1.5

        Npoints = round(Int64, 1.0 + (log(Zatoms[isp]*rmax)-xmin)/amesh)
        Npoints = round(Int64, Npoints/2)*2 + 1

        NptsRadial[isp] = Npoints

        r[isp] = zeros(Npoints)
        rab[isp] = zeros(Npoints)
        rho_at[isp] = zeros(Npoints)

        for ir = 1:Npoints
            rr = exp(xmin + (ir-1)*amesh)/Zatoms[isp]
            rab[isp][ir] = rr*amesh
            r[isp][ir] = rr
            rho_at[isp][ir] = znorml*exp(-rr^2/two_l2)*rr^2
        end
    end

    NpointsMax = maximum(NptsRadial)

    aux = zeros(NpointsMax)
    
    Ngl = length(Ham.pw.gvec.G2_shells)
    rhocgnt = zeros(Ngl)

    Ng = pw.gvec.Ng
    idx_g2shells = pw.gvec.idx_g2shells
    G2_shells = pw.gvec.G2_shells
    Npoints = prod(pw.Ns)
    idx_g2r = pw.gvec.idx_g2r
    CellVolume = pw.CellVolume
    Nspin = Ham.electrons.Nspin
    rhocg = zeros(ComplexF64,Npoints,Nspin)

    for isp = 1:Nspecies

        for ir = 1:NptsRadial[isp]
            aux[ir] = rho_at[isp][ir]
        end
        rhocgnt[1] = integ_simpson( NptsRadial[isp], aux, rab[isp] )

        for igl = 2:Ngl
            gx = sqrt(G2_shells[igl])
            for ir = 1:NptsRadial[isp]
                if r[isp][ir] < eps8
                   aux[ir] = rho_at[isp][ir]
                else
                   aux[ir] = rho_at[isp][ir]*sin(gx*r[isp][ir])/(r[isp][ir]*gx)
                end
            end
            rhocgnt[igl] = integ_simpson( NptsRadial[isp], aux, rab[isp] )
        end

        for ig = 1:Ng
            ip = idx_g2r[ig]
            rhocg[ip,1] = rhocg[ip,1] + strf[ig,isp]*rhocgnt[idx_g2shells[ig]]/CellVolume
        end

        if Nspin == 2
            for ig = 1:Ng
                ip = idx_g2r[ig]
                rhocg[ip,2] = rhocg[ip,2] + starting_magnetization[isp]*strf[ig,isp]*rhocgnt[idx_g2shells[ig]]/CellVolume
            end
        end
    end

    Rhoe = zeros(Npoints,Nspin)
    Rhoe_tot = real(G_to_R(pw,rhocg[:,1]))*Npoints
    
    # Rhoe_tot = Rhoe_up + Rhoe_dn
    # magn = Rhoe_up - Rhoe_dn
    # 2*Rhoe_up = Rhoe_tot + magn
    # 2*Rhoe_dn = Rhoe_tot - magn

    if Nspin == 2
        magn = real(G_to_R(pw,rhocg[:,2]))*Npoints
        Rhoe[:,1] = 0.5*(Rhoe_tot + magn)
        Rhoe[:,2] = 0.5*(Rhoe_tot - magn)
    else
        Rhoe[:,1] = Rhoe_tot
    end

    for ispin = 1:Nspin
        println("guess rhoe atomic: integ rhoe = ", sum(Rhoe[:,ispin])*CellVolume/Npoints)
    end

    return Rhoe
end

function integ_simpson(Npoints, f, rab)
#
# simpson's rule integration. On input:
#   mesh = the number of grid points (should be odd)
#   func(i)= function to be integrated
#   rab(i) = r(i) * dr(i)/di * di
#
# For the logarithmic grid not including r=0 :
#   r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
#
# For the logarithmic grid including r=0 :
#   r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
#
# Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr
# where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
#
    asum = 0.0
    r12 = 1.0/3.0
    f3 = f[1]*rab[1]*r12

    for i = 2:2:Npoints-1
        f1 = f3
        f2 = f[i]*rab[i]*r12
        f3 = f[i+1]*rab[i+1]*r12
        asum = asum + f1 + 4.0*f2 + f3
    end

    return asum
end