function guess_rhoe( Ham::Hamiltonian )
    atoms = Ham.atoms
    pw = Ham.pw
    pspots = Ham.pspots

    atm2species = atoms.atm2species
    idx_g2r = pw.gvec.idx_g2r
    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng

    Natoms = atoms.Natoms
    Npoints = prod(pw.Ns)
    Rhoe = zeros(Float64,Npoints)
    RhoeG = zeros(ComplexF64,Npoints)

    decay_lengths = get_atmlength( atoms, pspots )

    strf = calc_strfact( atoms, pw )

    println("Generating Gaussian Rhoe")

    for ia = 1:Natoms
        isp = atm2species[ia]
        Zval = pspots[isp].zval
        l = decay_lengths[isp]
        for ig = 1:Ng
            ip = idx_g2r[ig]
            RhoeG[ip] = RhoeG[ip] + Zval*strf[ig,isp]*exp(-G2[ig]*l^2)
        end
    end
    Rhoe = real( G_to_R(pw, RhoeG) )

    Nelectrons = Ham.electrons.Nelectrons
    CellVolume = pw.CellVolume

    intRhoe = sum(Rhoe)*CellVolume/Npoints
    Rhoe = Rhoe*Nelectrons/intRhoe
    intRhoe = sum(Rhoe)*CellVolume/Npoints
    @printf("Integrated Gaussian Rhoe: %18.10f\n", intRhoe)

    return Rhoe
end

function get_atmlength( atoms::Atoms, pspots::Array{PsPot_GTH} )
    Nspecies = atoms.Nspecies
    Lengths = zeros(Nspecies)
    Zatoms = get_Zatoms(atoms)
    for isp = 1:Nspecies
        Zval = pspots[isp].zval
        Lengths[isp] = atmlength( 0.0, Float64(Zval), Zatoms[isp] )
    end
    return Lengths
end


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

    if (Nspin == 2) && starting_magnetization==nothing
        starting_magnetization = 0.1*ones(Nspecies)
    end

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


#=*****************************************************************************************
!{\src2tex{textfont=tt}}
!!****f* ABINIT/atmlength
!! NAME
!! atmlength
!!
!! FUNCTION
!! Return atomic decay Length for one given type of atom.
!! This Length is used to generate an approximate atomic gaussian density
!! in reciprocal space:   n^AT(G)=exp[-(2pi.length.G)^2]
!!
!! COPYRIGHT
!! Copyright (C) 2000-2010 ABINIT group (DCA,XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! densty=parameter for initialisation of the density of this atom type
!!        if densty>0, returned decay Length if densty !
!! zion=charge on current type of atom (real number)
!! znucl=atomic number, for current type of atom
!!
!! OUTPUT
!! length=decay length
!!
!! PARENTS
!!      initro
!!
!! CHILDREN
!!
!! SOURCE

modified by Fadjar Fathurrahman
=#

function atmlength( densty::Float64, zion::Float64, znucl::Float64 )
  
    tol10  = 1.e-8

    data_Length = zeros(16)

    # Either use the input value, or the default value, tabulated now.
    if abs(densty) > tol10 
        return densty
    end

    # Count the number of core electrons.
    coreel = znucl - zion
    
    # Round the number of valence electrons
    nval = round(Int64, zion) #nint(zion)

    Length = 0.0

    # For each set of core electron numbers, there are different decay lengths,
    # they start from nval=1, and proceed by group of 5, until a default is used    
    if nval == 0
        Length = 0.0
    
    # Bare ions : adjusted on 1H and 2He only
    elseif coreel < 0.5
        data_Length[1:4] = [0.6, 0.4, 0.3, 0.25]
        Length = 0.2
        if nval <= 4
            Length = data_Length[nval]
        end
    
    # 1s2 core : adjusted on 3Li, 6C, 7N, and 8O
    elseif coreel < 2.5
        data_Length[1:8] = [1.8, 1.4, 1.0, 0.7, 0.6, 0.5, 0.4, 0.35]
        Length = 0.3
        if nval <= 8
            Length = data_Length[nval]
        end

    # Ne core (1s2 2s2 2p6) : adjusted on 11na, 13al, 14si and 17cl
    elseif coreel < 10.5
        data_Length[1:10] = [2.0, 1.6, 1.25, 1.1, 1.0,
                             0.9, 0.8, 0.7 , 0.7, 0.7]
        Length = 0.6
        if nval <= 10
            Length = data_Length[nval]
        end

    # Mg core (1s2 2s2 2p6 3s2) : adjusted on 19k, and on coreel==10
    elseif coreel < 12.5
        data_Length[1:10] = [1.9, 1.5, 1.15, 1.0, 0.9,
                             0.8, 0.7, 0.6 , 0.6, 0.6]
        Length = 0.5
        if nval <= 10
            Length = data_Length[nval]
        end

    # Ar core (Ne + 3s2 3p6) : adjusted on 20ca, 25mn and 30zn
    elseif coreel < 18.5
        data_Length[1:12] = [2.0,  1.8,  1.5, 1.2,  1.0,
                             0.9,  0.85, 0.8, 0.75, 0.7,
                             0.65, 0.65]
        Length = 0.6
        if nval <= 12
            Length = data_Length[nval]
        end

    # Full 3rd shell core (Ar + 3d10) : adjusted on 31ga, 34se and 38sr
    elseif coreel < 28.5
        data_Length[1:14] = [1.5 , 1.25, 1.15, 1.05, 1.00,
                             0.95, 0.95, 0.9 , 0.9 , 0.85,
                             0.85, 0.80, 0.8 , 0.75]
        Length = 0.7
        if nval <= 14
            Length = data_Length[nval]
        end

    # Krypton core (Ar + 3d10 4s2 4p6) : adjusted on 39y, 42mo and 48cd
    elseif coreel < 36.5
        data_Length[1:12] = [2.0 , 2.00, 1.60, 1.40, 1.25,
                             1.10, 1.00, 0.95, 0.90, 0.85,
                             0.80, 0.75]
        Length = 0.7
        if nval <= 12
            Length = data_Length[nval]
        end

    # For the remaining elements, consider a function of nval only
    else
        data_Length[1:12] = [2.0 , 2.00, 1.55, 1.25, 1.15,
                             1.10, 1.05, 1.0 , 0.95, 0.9,
                             0.85, 0.85]
        Length = 0.8
        if nval <= 12
            Length = data_Length[nval]
        end
    end

    return Length
end
