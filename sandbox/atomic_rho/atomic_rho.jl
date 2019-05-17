function guess_rhoe_atomic( Ham::Hamiltonian )

    atoms = Ham.atoms
    Nspecies = atoms.Nspecies

    pspots = Ham.pspots

    xmin = -7.0
    amesh = 0.0125
    rmax = 100.0

    Zatoms = get_Zatoms( atoms )

    r = Array{Array{Float64,1},1}(undef,Nspecies)
    rab = Array{Array{Float64,1},1}(undef,Nspecies)
    rho_at = Array{Array{Float64,1},1}(undef,Nspecies)

    decay_lengths = PWDFT.get_atmlength( atoms, pspots )

    NptsRadial = zeros(Int64,Nspecies)

    for isp = 1:Nspecies

        two_l2 = 2.0*decay_lengths[isp]^2
        znorml = 4.0*pi*pspots[isp].zval/(pi*two_l2)^1.5

        Npoints = round( Int64, 1.0 + (log(Zatoms[isp]*rmax)-xmin)/amesh )
        Npoints = round( Int64, (Npoints/2)*2 + 1 )

        NptsRadial[isp] = Npoints

        println(isp, " ", Npoints)

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
    println("NpointsMax = ", NpointsMax)

    aux = zeros(NpointsMax)
    
    Ng = Ham.pw.gvec.Ng
    rhocgnt = zeros(Ng)

    for isp = 1:Nspecies
        for ir = 1:NptsRadial[isp]
            aux[ir] = rho_at[isp][ir]
        end
        rhocgnt[1] = integ_simpson(NptsRadial[isp], aux, rab[isp])
    end

    return
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