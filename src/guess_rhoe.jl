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
