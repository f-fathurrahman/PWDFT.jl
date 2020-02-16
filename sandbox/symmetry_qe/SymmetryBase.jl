# Based on the file PW/src/symm_base.f90 of Quantum ESPRESS
#
# Copyright (C) 2010-2017 Quantum ESPRESSO group
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#

using Printf
using LinearAlgebra

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

const GBL_eps1 = 1.0e-6
const GBL_eps2 = 1.0e-5

const GBL_sin3  =  0.866025403784438597
const GBL_cos3  =  0.5
const GBL_msin3 = -0.866025403784438597
const GBL_mcos3 = -0.5

const GBL_s0 = reshape([
    1.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  1.0,
   -1.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,  1.0,
   -1.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0, -1.0,
    1.0,  0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0, -1.0,
    0.0,  1.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0, -1.0,
    0.0, -1.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0, -1.0,
    0.0, -1.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  1.0,
    0.0,  1.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0,
    0.0,  0.0,  1.0,  0.0, -1.0,  0.0,  1.0,  0.0,  0.0,
    0.0,  0.0, -1.0,  0.0, -1.0,  0.0, -1.0,  0.0,  0.0,
    0.0,  0.0, -1.0,  0.0,  1.0,  0.0,  1.0,  0.0,  0.0,
    0.0,  0.0,  1.0,  0.0,  1.0,  0.0, -1.0,  0.0,  0.0,
   -1.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  1.0,  0.0,
   -1.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0, -1.0,  0.0,
    1.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  1.0,  0.0,
    1.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0, -1.0,  0.0,
    0.0,  0.0,  1.0,  1.0,  0.0,  0.0,  0.0,  1.0,  0.0,
    0.0,  0.0, -1.0, -1.0,  0.0,  0.0,  0.0,  1.0,  0.0,
    0.0,  0.0, -1.0,  1.0,  0.0,  0.0,  0.0, -1.0,  0.0,
    0.0,  0.0,  1.0, -1.0,  0.0,  0.0,  0.0, -1.0,  0.0,
    0.0,  1.0,  0.0,  0.0,  0.0,  1.0,  1.0,  0.0,  0.0,
    0.0, -1.0,  0.0,  0.0,  0.0, -1.0,  1.0,  0.0,  0.0,
    0.0, -1.0,  0.0,  0.0,  0.0,  1.0, -1.0,  0.0,  0.0,
    0.0,  1.0,  0.0,  0.0,  0.0, -1.0, -1.0,  0.0,  0.0,
   GBL_cos3,  GBL_sin3, 0.0, GBL_msin3,  GBL_cos3, 0.0, 0.0, 0.0,  1.0,
   GBL_cos3, GBL_msin3, 0.0,  GBL_sin3,  GBL_cos3, 0.0, 0.0, 0.0,  1.0,
   GBL_mcos3,  GBL_sin3, 0.0, GBL_msin3, GBL_mcos3, 0.0, 0.0, 0.0,  1.0,
   GBL_mcos3, GBL_msin3, 0.0,  GBL_sin3, GBL_mcos3, 0.0, 0.0, 0.0,  1.0,
   GBL_cos3, GBL_msin3, 0.0, GBL_msin3, GBL_mcos3, 0.0, 0.0, 0.0, -1.0,
   GBL_cos3,  GBL_sin3, 0.0,  GBL_sin3, GBL_mcos3, 0.0, 0.0, 0.0, -1.0,
   GBL_mcos3, GBL_msin3, 0.0, GBL_msin3,  GBL_cos3, 0.0, 0.0, 0.0, -1.0,
   GBL_mcos3,  GBL_sin3, 0.0,  GBL_sin3,  GBL_cos3, 0.0, 0.0, 0.0, -1.0
], (3,3,32))

function find_symm_bravais_latt( LatVecs_ )

    alat = norm(LatVecs_[:,1])
    LatVecs = LatVecs_[:,:]/alat

    rot = zeros(Float64,3,3)
    overlap = zeros(Float64,3,3)

    # calculate LatVecs' * LatVecs
    for j in 1:3
        for k in 1:3
            rot[k,j] = LatVecs[1,k]*LatVecs[1,j] +
                       LatVecs[2,k]*LatVecs[2,j] +
                       LatVecs[3,k]*LatVecs[3,j]
        end
    end

    overlap = inv(rot)

    #display(LatVecs); println()
    #display(LatVecs'*LatVecs); println()
    #display(rot); println()
    #println("\noverlap matrix: ")
    #display(overlap); println()

    #display(GBL_s0[:,:,1]); println()
    #display(GBL_s0[:,:,2]); println()


    imat = zeros(Int64,32)

    rat = zeros(Float64,3)
    rot = zeros(Float64,3,3)

    s = zeros(Int64,3,3,48)

    nrot = 1
    for irot in 1:32

        # for each possible symmetry
        for jpol in 1:3
            for mpol in 1:3
                # Cartesian coordinates the rotated vector
                rat[mpol] = GBL_s0[mpol,1,irot]*LatVecs[1,jpol] +
                            GBL_s0[mpol,2,irot]*LatVecs[2,jpol] +
                            GBL_s0[mpol,3,irot]*LatVecs[3,jpol]
            end

            for kpol in 1:3
                # The rotated vector is projected on the direct lattice
                rot[kpol,jpol] = LatVecs[1,kpol]*rat[1] +
                                 LatVecs[2,kpol]*rat[2] +
                                 LatVecs[3,kpol]*rat[3]
            end
        end

        status = _check_set_s!( overlap, rot, s, nrot )

        if status == -1
            continue
        end

        imat[nrot] = irot
        nrot = nrot + 1

    end

    println("nrot after loop = ", nrot)

    nrot = nrot - 1

    if !(nrot in [1, 2, 4, 6, 8, 12, 24])
        println("WARNING: Bravais lattice has wrong number of symmetries: ", nrot)
        println("symmetries are disabled")
        nrot = 1
    end

    #display(s[:,:,2]); println()

    # Set the inversion symmetry (Bravais lattices have always inversion symmetry)
    for irot = 1:nrot
        #sname[irot+nrot] = s0name(imat(irot)+32)
        for kpol in 1:3
            for jpol in 1:3
                s[kpol,jpol,irot+nrot] = -s[kpol,jpol,irot]
            end
        end
    end


    nrot = 2*nrot
    
    # Reset fractional translations to zero before checking the group
    ft = zeros(Float64,3,48)    

    if !is_group(nrot, s, ft)
        # ... This happens for instance for an hexagonal lattice with one axis 
        # oriented at 15 degrees from the x axis, the other along (-1,1,0)
        println("NOTICE: Symmetry group for Bravais lattic is not a group")
        nrot = 1
    end

    println("nrot = ", nrot)

end

function _check_set_s!( overlap, rot, s, nrot )

    # The inverse of the overlap matrix is applied
    for jpol in 1:3
        for kpol in 1:3
            value = overlap[jpol,1]*rot[1,kpol] +
                    overlap[jpol,2]*rot[2,kpol] +
                    overlap[jpol,3]*rot[3,kpol]
            if( abs(round(value) - value) > GBL_eps1 )
                # If a noninteger is obtained, this implies that this operation
                # is not a symmetry operation for the given lattice
                return -1
            end
            #
            s[kpol,jpol,nrot] = round(Int64, value)
        end
    end

    return 0
end


function is_group( nsym_::Int64, s, ft )
    # Checks that {S} is a group.
    #
    #!
    #INTEGER, INTENT(IN) :: nsym_
    #INTEGER :: isym, jsym, ksym, ss(3,3)
    #REAL(DP) :: st(3), dt(3)
    #LOGICAL :: found
    #!

    ss = zeros(Int64,3,3)
    st = zeros(Float64,3)
    dt = zeros(Float64,3)   

    for isym in 1:nsym_
        for jsym in 1:nsym_
            
            ss = s[:,:,isym] * s[:,:,jsym]

            st[:] = ft[:,jsym] + s[1,:,jsym]*ft[1,isym] +
                                 s[2,:,jsym]*ft[2,isym] +
                                 s[3,:,jsym]*ft[3,isym]
            #
            # ... here we check that the input matrices really form a group:
            # S(k) = S(i)*S(j)
            # ftau_k = S(j)*ftau_i+ftau_j (modulo a lattice vector)
            #
            found = false
            
            for ksym in 1:nsym_
                dt[1] = ft[1,ksym] - st[1] - round( ft[1,ksym] - st[1] )
                dt[2] = ft[2,ksym] - st[2] - round( ft[2,ksym] - st[2] )
                dt[3] = ft[3,ksym] - st[3] - round( ft[3,ksym] - st[3] )

                cond1 = all( s[:,:,ksym] .== ss[:,:] )
                cond2 = all( abs.(dt .< GBL_eps2) )

                if ( cond1 && cond2 )
                    if found
                        return false
                    end
                    found = true
                end
            end
            
            if !found
               return false
            end

        end
    end
    
    return true

end


function init_Ham_Si_fcc( meshk::Array{Int64,1} )
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.1  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.431*ANG2BOHR))
    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk )
end

function main()
    Ham = init_Ham_Si_fcc([3,3,3])

    find_symm_bravais_latt( Ham.atoms.LatVecs )

    println("Pass here")
end

main()