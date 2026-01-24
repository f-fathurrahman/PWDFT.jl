# Based on the file PW/src/symm_base.f90 of Quantum ESPRESSO
#
# Copyright (C) 2010-2017 Quantum ESPRESSO group
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#

const _GBL_eps1 = 1.0e-6
const _GBL_eps2 = 1.0e-5

const _GBL_sin3  =  0.866025403784438597
const _GBL_cos3  =  0.5
const _GBL_msin3 = -0.866025403784438597
const _GBL_mcos3 = -0.5

const _GBL_s0 = reshape([
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
   _GBL_cos3,  _GBL_sin3,  0.0, _GBL_msin3, _GBL_cos3,  0.0, 0.0, 0.0,  1.0,
   _GBL_cos3,  _GBL_msin3, 0.0, _GBL_sin3,  _GBL_cos3,  0.0, 0.0, 0.0,  1.0,
   _GBL_mcos3, _GBL_sin3,  0.0, _GBL_msin3, _GBL_mcos3, 0.0, 0.0, 0.0,  1.0,
   _GBL_mcos3, _GBL_msin3, 0.0, _GBL_sin3,  _GBL_mcos3, 0.0, 0.0, 0.0,  1.0,
   _GBL_cos3,  _GBL_msin3, 0.0, _GBL_msin3, _GBL_mcos3, 0.0, 0.0, 0.0, -1.0,
   _GBL_cos3,  _GBL_sin3,  0.0, _GBL_sin3,  _GBL_mcos3, 0.0, 0.0, 0.0, -1.0,
   _GBL_mcos3, _GBL_msin3, 0.0, _GBL_msin3, _GBL_cos3,  0.0, 0.0, 0.0, -1.0,
   _GBL_mcos3, _GBL_sin3,  0.0, _GBL_sin3,  _GBL_cos3,  0.0, 0.0, 0.0, -1.0
], (3,3,32))

const _GBL_s0name = [
    "identity                                     ",
    "180 deg rotation - cart. axis [0,0,1]        ",
    "180 deg rotation - cart. axis [0,1,0]        ",
    "180 deg rotation - cart. axis [1,0,0]        ",
    "180 deg rotation - cart. axis [1,1,0]        ",
    "180 deg rotation - cart. axis [1,-1,0]       ",
    " 90 deg rotation - cart. axis [0,0,-1]       ",
    " 90 deg rotation - cart. axis [0,0,1]        ",
    "180 deg rotation - cart. axis [1,0,1]        ",
    "180 deg rotation - cart. axis [-1,0,1]       ",
    " 90 deg rotation - cart. axis [0,1,0]        ",
    " 90 deg rotation - cart. axis [0,-1,0]       ",
    "180 deg rotation - cart. axis [0,1,1]        ",
    "180 deg rotation - cart. axis [0,1,-1]       ",
    " 90 deg rotation - cart. axis [-1,0,0]       ",
    " 90 deg rotation - cart. axis [1,0,0]        ",
    "120 deg rotation - cart. axis [-1,-1,-1]     ",
    "120 deg rotation - cart. axis [-1,1,1]       ",
    "120 deg rotation - cart. axis [1,1,-1]       ",
    "120 deg rotation - cart. axis [1,-1,1]       ",
    "120 deg rotation - cart. axis [1,1,1]        ",
    "120 deg rotation - cart. axis [-1,1,-1]      ",
    "120 deg rotation - cart. axis [1,-1,-1]      ",
    "120 deg rotation - cart. axis [-1,-1,1]      ",
    " 60 deg rotation - cryst. axis [0,0,1]       ",
    " 60 deg rotation - cryst. axis [0,0,-1]      ",
    "120 deg rotation - cryst. axis [0,0,1]       ",
    "120 deg rotation - cryst. axis [0,0,-1]      ",
    "180 deg rotation - cryst. axis [1,-1,0]      ",
    "180 deg rotation - cryst. axis [2,1,0]       ",
    "180 deg rotation - cryst. axis [0,1,0]       ",
    "180 deg rotation - cryst. axis [1,1,0]       ",
    "inversion                                    ",
    "inv. 180 deg rotation - cart. axis [0,0,1]   ",
    "inv. 180 deg rotation - cart. axis [0,1,0]   ",
    "inv. 180 deg rotation - cart. axis [1,0,0]   ",
    "inv. 180 deg rotation - cart. axis [1,1,0]   ",
    "inv. 180 deg rotation - cart. axis [1,-1,0]  ",
    "inv.  90 deg rotation - cart. axis [0,0,-1]  ",
    "inv.  90 deg rotation - cart. axis [0,0,1]   ",
    "inv. 180 deg rotation - cart. axis [1,0,1]   ",
    "inv. 180 deg rotation - cart. axis [-1,0,1]  ",
    "inv.  90 deg rotation - cart. axis [0,1,0]   ",
    "inv.  90 deg rotation - cart. axis [0,-1,0]  ",
    "inv. 180 deg rotation - cart. axis [0,1,1]   ",
    "inv. 180 deg rotation - cart. axis [0,1,-1]  ",
    "inv.  90 deg rotation - cart. axis [-1,0,0]  ",
    "inv.  90 deg rotation - cart. axis [1,0,0]   ",
    "inv. 120 deg rotation - cart. axis [-1,-1,-1]",
    "inv. 120 deg rotation - cart. axis [-1,1,1]  ",
    "inv. 120 deg rotation - cart. axis [1,1,-1]  ",
    "inv. 120 deg rotation - cart. axis [1,-1,1]  ",
    "inv. 120 deg rotation - cart. axis [1,1,1]   ",
    "inv. 120 deg rotation - cart. axis [-1,1,-1] ",
    "inv. 120 deg rotation - cart. axis [1,-1,-1] ",
    "inv. 120 deg rotation - cart. axis [-1,-1,1] ",
    "inv.  60 deg rotation - cryst. axis [0,0,1]  ",
    "inv.  60 deg rotation - cryst. axis [0,0,-1] ",
    "inv. 120 deg rotation - cryst. axis [0,0,1]  ",
    "inv. 120 deg rotation - cryst. axis [0,0,-1] ",
    "inv. 180 deg rotation - cryst. axis [1,-1,0] ",
    "inv. 180 deg rotation - cryst. axis [2,1,0]  ",
    "inv. 180 deg rotation - cryst. axis [0,1,0]  ",
    "inv. 180 deg rotation - cryst. axis [1,1,0]  "
]

"""
    find_symm_bravais_latt( LatVecs_ )

Find symmetry operations of a given Bravais lattice defined by `LatVecs_`.
"""
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

    imat = zeros(Int64,32)

    rat = zeros(Float64,3)
    rot = zeros(Float64,3,3)

    s = zeros(Int64,3,3,48)
    sname = [repeat(" ", 45) for i in 1:48]

    nrot = 1
    for irot in 1:32

        # for each possible symmetry
        for jpol in 1:3
            for mpol in 1:3
                # Cartesian coordinates the rotated vector
                rat[mpol] = _GBL_s0[mpol,1,irot]*LatVecs[1,jpol] +
                            _GBL_s0[mpol,2,irot]*LatVecs[2,jpol] +
                            _GBL_s0[mpol,3,irot]*LatVecs[3,jpol]
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

        sname[nrot] = _GBL_s0name[irot]
        imat[nrot] = irot
        nrot = nrot + 1

    end

    nrot = nrot - 1

    if !(nrot in [1, 2, 4, 6, 8, 12, 24])
        println("WARNING: Bravais lattice has wrong number of symmetries: ", nrot)
        println("symmetries are disabled")
        nrot = 1
    end

    # Set the inversion symmetry (Bravais lattices have always inversion symmetry)
    for irot = 1:nrot
        sname[irot+nrot] = _GBL_s0name[imat[irot]+32]
        for kpol in 1:3
            for jpol in 1:3
                s[kpol,jpol,irot+nrot] = -s[kpol,jpol,irot]
            end
        end
    end

    nrot = 2*nrot
    
    # Reset fractional translations to zero before checking the group
    ft_zero = zeros(Float64, 3, 48)    

    if !is_group(nrot, s, ft_zero)
        # This happens for instance for an hexagonal lattice with one axis 
        # oriented at 15 degrees from the x axis, the other along (-1,1,0)
        println("NOTICE: Symmetry group for Bravais lattice is not a group")
        nrot = 1
    end

    return nrot, s, sname

end

function _check_set_s!( overlap, rot, s, nrot )

    # The inverse of the overlap matrix is applied
    for jpol in 1:3
        for kpol in 1:3
            value = overlap[jpol,1]*rot[1,kpol] +
                    overlap[jpol,2]*rot[2,kpol] +
                    overlap[jpol,3]*rot[3,kpol]
            if( abs(round(value) - value) > _GBL_eps1 )
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


# Given the point group of the Bravais lattice, this routine finds
# the subgroup which is the point group of the considered crystal.  
# Non symmorphic groups are allowed, provided that fractional
# translations are allowed (nofrac=.false) and that the unit cell
# is not a supercell.
# On output, the array sym is set to .TRUE.. for each operation
# of the original point group that is also a symmetry operation
# of the crystal symmetry point group.
function sgam_at!(
    atoms::Atoms, sym, nrot, s, ft, irt;
    no_z_inv = false,
    fractional_translations = true
)
    
    nat = atoms.Natoms
    tau = atoms.positions
    atm2species = atoms.atm2species

    bg = inv(Matrix(atoms.LatVecs'))  # follows QE convention (no 2*pi factor)

    xau = zeros(Float64,3,nat)
    rau = zeros(Float64,3,nat)
    ft_ = zeros(Float64,3)

    # Compute the coordinates of each atom in the basis of the direct lattice vectors
    for ia = 1:nat
        for i = 1:3
            xau[i,ia] = bg[1,i]*tau[1,ia] + bg[2,i]*tau[2,ia] + bg[3,i]*tau[3,ia]
        end
    end

    #
    # Check if the identity has fractional translations (this means
    # that the cell is actually a supercell). When this happens, fractional
    # translations are disabled, because there is no guarantee that the 
    # generated sym.ops. form a group.
    nb = 1
    irot = 1

    if fractional_translations
        for na in 2:nat
            if atm2species[nb] != atm2species[na]
                continue
            end 
            for i in 1:3
                ft_[i] = xau[i,na] - xau[i,nb] - round( xau[i,na] - xau[i,nb] )
            end
            sym[irot] = checksym!( irot, nat, atm2species, xau, xau, ft_, irt )    
            if sym[irot]
                fractional_translations = false
                #println("Found symmetry operation: I + ", ft_[:])
                #println("This is a supercell")
                #println("Fractional translations are disabled")
                break # go outside the loop over na
            end        
        end # for na
    end # if

    nsym_ns = 0
    fft_fact = [1,1,1]
    ftaux = zeros(Float64,3) 

    for irot in 1:nrot

        # rau is rotated atom coordinates
        for na in 1:nat
            for i in 1:3
               rau[i,na] = s[1,i,irot] * xau[1,na] + s[2,i,irot] * xau[2,na] + s[3,i,irot] * xau[3,na]
            end
        end
        
        # first attempt: no fractional translation
        for i in 1:3
            ft[i,irot] = 0.0
            ft_[i] = 0.0
        end
        
        sym[irot] = checksym!( irot, nat, atm2species, xau, rau, ft_, irt )
        
        if ( !sym[irot] && fractional_translations )
            nb = 1
            for na in 1:nat
                if atm2species[nb] == atm2species[na]
                 
                    # Second attempt: check all possible fractional translations
                    for i in 1:3
                        ft_[i] = rau[i,na] - xau[i,nb] - round( rau[i,na] - xau[i,nb] )
                    end
                    
                    #  ft_ is in crystal axis and is a valid fractional translation
                    # only if ft_(i)=0 or ft_(i)=1/n, with n=2,3,4,
                    for i in 1:3
                        if abs(ft_[i]) > _GBL_eps2
                            ftaux[i] = abs( 1.0/ft_[i] - round(Int64, 1.0/ft_[i]) )
                            nfrac = round( Int64, 1.0/abs(ft_[i]) )
                            if( (ftaux[i] < _GBL_eps2) && !(nfrac in [2, 3, 4, 6]) )
                                ftaux[i] = 2*_GBL_eps2
                            end # if
                        else
                            ftaux[i] = 0.0
                        end # if
                    end # for
                    
                    if( any( ftaux .> _GBL_eps2 ) )
                        continue
                    end
                    
                    sym[irot] = checksym!( irot, nat, atm2species, xau, rau, ft_, irt )
                 
                    if sym[irot]
                        nsym_ns = nsym_ns + 1
                        ft[:,irot] = ft_[:]
                    
                        # Find factors that must be present in FFT grid dimensions
                        # in order to ensure that fractional translations are
                        # commensurate with FFT grids.
                        for i in 1:3
                            if abs(ft_[i]) > _GBL_eps2
                                nfrac = round( Int64, 1.0/abs(ft_[i]) )
                            else
                                nfrac = 0
                            end # if
                            fft_fact[i] = _mcm( fft_fact[i], nfrac)
                        end # do

                        break # break from loop over na
                    end # if
                end # if

            end # do
           
        end # if

    end

    # disable z -> -z symmetry
    if no_z_inv
        for irot in 1:nrot
            if s[3,3,irot] == -1
                sym[irot] = false
            end
        end
    end

    return

end


function sgam_at_mag!( atoms::Atoms, m_loc, sym, nrot, s, ft, sname, irt, t_rev, nsym_ns )
#=
    Find magnetic symmetries, i.e. point-group symmetries that are
    also symmetries of the local magnetization - including
    rotation + time reversal operations.
=#

    EPS2 = 1e-5

    Natoms = atoms.Natoms
    mxau = zeros(Float64, 3, Natoms) # magnetization
    mrau = zeros(Float64, 3, Natoms) # rotated magnetization

    no_t_rev = false # XXX
    
    bg = inv(Matrix(atoms.LatVecs'))  # follows QE convention (no 2*pi factor)
    # The exact value of bg does not need to coincide with QE

    println("in sgam_at_mag!: m_loc = ", m_loc)

    # Compute the local magnetization of each atom in the basis of
    # the direct lattice vectors
    for ia in 1:Natoms
        mxau[:,ia] = bg[1,:]*m_loc[1,ia] + bg[2,:]*m_loc[2,ia] + bg[3,:]*m_loc[3,ia]
        println("ia = $ia mxau = $mxau")
    end
    
    for irot in 1:nrot
        #
        println("\nirot = ", irot)
        #
        t_rev[irot] = false
        if sym[irot]
            # rotated local magnetization
            for ia in 1:Natoms
                mrau[:,ia]= s[1,:,irot]*mxau[1,ia] +
                            s[2,:,irot]*mxau[2,ia] +
                            s[3,:,irot]*mxau[3,ia]
                println("ia=$ia mrau = $mrau")
            end
            if sname[irot][1:3] == "inv"
                @info "inversion sym, set mrau to -mrau"
                mrau[:,:] = -mrau[:,:]
            end
            # check if this a magnetic symmetry
            t1 = true
            t2 = true  
            for ia in 1:Natoms
                ja = irt[irot,ia]
                if (ja < 1) || (ja > Natoms)
                    error("out-of-bound atomic index: $ja")
                end
                ss1 = sum(abs.(mrau[:,ia] - mxau[:,ja]))
                ss2 = sum(abs.(mrau[:,ia] + mxau[:,ja]))
                t1 = ( ss1 < EPS2 ) && t1
                t2 = ( ss2 < EPS2 ) && t2
                println("ss1=$ss1 t1=$t1 , ss2=$ss2 t2=$t2")
            end
           
            if (!t1 && !t2)
                # not a magnetic symmetry
                sym[irot] = false
                println("sym irot is set to false: ", irot)
            elseif (t2 && !t1)
                # magnetic symmetry with time reversal, if allowed
                if no_t_rev
                    sym[irot] = false
                else
                    t_rev[irot] = true
                end
            end # if 
            #
            if ( !sym[irot] && any( abs.(ft[:,irot]) .> EPS2 ) )
                nsym_ns -= 1
            end
        end # if
    end # nrot

    @info "at the end of sgam_at_mag!: $(count(sym))"

    return nsym_ns
end



# Returns minimum common multiple of two integers
# if i=0, returns j, and vice versa; if i<0 or j<0, returns -1.
function _mcm( i, j )
    res = -1
    if (i < 0) || (j < 0)
        res = -1
    elseif (i == 0) && (j == 0)
        res = 0
    else
        n1 = min(i,j)
        n2 = max(i,j)
        for k in 1:n1
            res = k*n2 
            if res%n1 == 0
                return res
            end
        end
        res = n2
    end
    return res
end

# This function receives as input all the atomic positions xau,
# and the rotated rau by the symmetry operation ir. It returns
# .TRUE. if, for each atom na, it is possible to find an atom nb
# which is of the same type of na, and coincides with it after the
# symmetry operation. Fractional translations are allowed.
function checksym!( irot, Natoms, ityp, xau, rau, ft_, irt )

    ACCEPT = 1e-5

    continue_ia = zeros(Bool, Natoms)
    
    for ia in 1:Natoms

        for ib in 1:Natoms
            
            if ityp[ib] == ityp[ia]
                is_equal =  eqvect( rau[:,ia], xau[:,ib], ft_ , ACCEPT )
                if is_equal
                   # the rotated atom does coincide with one of the like atoms
                   # keep track of which atom the rotated atom coincides with
                   irt[irot,ia] = ib
                   continue_ia[ia] = true
                   break # move outside loop over nb
                end # if
            end # if

        end # for
        
        # The rotated atom does not coincide with any of the like atoms
        # s(ir) + ft is not a symmetry operation
        if !continue_ia[ia]
            return false
        else
            continue
        end
    end

    # s(ir) + ft is a symmetry operation
    return true
end


# This function test if the difference x-y-f is an integer.
# x, y = 3d vectors in crystal axis, f = fractional translation
# adapted from QE-6.4
function eqvect(x, y, f, accep)

    cond1 = ( abs(x[1]-y[1]-f[1] - round(x[1]-y[1]-f[1]) ) < accep )
    cond2 = ( abs(x[2]-y[2]-f[2] - round(x[2]-y[2]-f[2]) ) < accep )
    cond3 = ( abs(x[3]-y[3]-f[3] - round(x[3]-y[3]-f[3]) ) < accep )

    res = cond1 && cond2 && cond3
  
    return res
end



# Copy symmetry operations in sequential order so that:
#
# * \(s(i,j,\text{irot})\), with \(\text{irot} \leq \text{nsym}\) are the symmetry
#   operations of the crystal;
# * \(s(i,j,\text{irot})\), with \(\text{nsym}+1<\text{irot}\leq \text{nrot}\) are 
#   the symmetry operations of the lattice.
#
# On exit \(\textrm{copy_sym}\) returns nsym.
function copy_sym!( nrot_, sym, s, ft, sname, irt )
    
    stemp = zeros(Int64,3,3)
    ft_ = zeros(Float64,3)

    Natoms = size(irt,2)
    irtemp = zeros(Int64, Natoms )
    
    jrot = 0
    for irot in 1:nrot_
        if sym[irot]
            jrot = jrot + 1
            if irot > jrot
                for j in 1:3, i in 1:3
                    stemp[i,j] = s[i,j,jrot]
                    s[i,j,jrot] = s[i,j,irot]
                    s[i,j,irot] = stemp[i,j]
                end
                
                for i in 1:3
                    ft_[i] = ft[i,jrot]
                    ft[i,jrot] = ft[i,irot]
                    ft[i,irot] = ft_[i]
                end
                
                for ia in 1:Natoms
                    irtemp[ia] = irt[jrot,ia]
                    irt[jrot,ia] = irt[irot,ia]
                    irt[irot,ia] = irtemp[ia]
                end
                
                nametemp = sname[jrot]
                sname[jrot] = sname[irot]
                sname[irot] = nametemp

            end # if
        end # if
    end # for
    
    for i in 1:jrot
        sym[i] = true
    end
    for i in jrot+1:nrot_
        sym[i] = false
    end

    return jrot
end


# Checks that {S} is a group.
function is_group( nsym_::Int64, s, ft )

    ss = zeros(Int64,3,3)
    st = zeros(Float64,3)
    dt = zeros(Float64,3)   

    for isym in 1:nsym_
        for jsym in 1:nsym_
            
            ss = s[:,:,isym] * s[:,:,jsym]

            for i in 1:3
                st[i] = ft[i,jsym] + s[1,i,jsym]*ft[1,isym] +
                                     s[2,i,jsym]*ft[2,isym] +
                                     s[3,i,jsym]*ft[3,isym]
            end
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
                cond2 = all( abs.(dt .< _GBL_eps2) )

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


struct SymmetryInfo
    Nrots::Int64
    Nsyms::Int64
    s::Array{Int64,3}
    inv_s::Array{Int64,3}
    idx_inv_s::Vector{Int64}
    sname::Array{String}
    sr::Array{Float64,3}
    ft::Array{Float64,2}
    non_symmorphic::Vector{Bool}
    t_rev::Vector{Bool}
    irt::Array{Int64,2}
    D1::Array{Float64,3}
    D2::Array{Float64,3}
    D3::Array{Float64,3}
end

function SymmetryInfo()
    Nrots = 1
    Nsyms = 1
    s = zeros(Int64,3,3,1)
    s[1,1,1] = 1
    s[2,2,1] = 1
    s[3,3,1] = 1
    #
    sr = zeros(Float64,3,3,1) # XXX: This should not be used
    # We need Atoms.LatVecs to properly set up this
    #
    inv_s = copy(s)
    idx_inv_s = [1]
    #
    sname = ["identity"]
    #
    ft = zeros(3,1)
    #
    non_symmorphic = [false]
    t_rev = [false]
    #
    irt = zeros(Int64,1,1)
    #
    D1 = zeros(1,1,1)
    D2 = zeros(1,1,1)
    D3 = zeros(1,1,1)
    #
    return SymmetryInfo( Nrots, Nsyms, s, inv_s, idx_inv_s, sname, sr, ft, non_symmorphic, t_rev, irt, D1, D2, D3 )
end

function SymmetryInfo( atoms::Atoms; magnetic_sym = false, m_loc = nothing )

    Nrots, s, sname = find_symm_bravais_latt( atoms.LatVecs )

    # Allocate memory
    sym = zeros(Bool, 48)
    irt = zeros(Int64, 48, atoms.Natoms)
    ft = zeros(Float64, 3, 48)

    sgam_at!( atoms, sym, Nrots, s, ft, irt,
        fractional_translations = true, no_z_inv = false
    )

    nsym_ns = 0 # nonsymmorphic (fractional translation) symms
    t_rev = zeros(Bool, 48)
    if magnetic_sym
        @assert !isnothing(m_loc)
        nsym_ns = sgam_at_mag!( atoms, m_loc, sym, Nrots, s, ft, sname, irt, t_rev, nsym_ns )
    end

    Nsyms = copy_sym!( Nrots, sym, s, ft, sname, irt )
    
    inv_s = zeros(Int64,3,3,Nsyms)
    for isym = 1:Nsyms
        inv_s[:,:,isym] = Base.convert(Array{Int64,2}, inv(s[:,:,isym]))
    end

    idx_inv_s = zeros(Int64, Nsyms)
    ss = zeros(Int64, 3, 3)
    for isym in 1:Nsyms
        found = false
        for jsym in 1:Nsyms
            ss[:,:] = s[:,:,jsym] * s[:,:,isym]
            # s(:,:,1) is the identity
            if all( s[:,:,1] .== ss )
                idx_inv_s[isym] = jsym
                found = true
            end
        end
        if !found
            @error("idx_inv_s is not found")
        end
    end

    non_symmorphic = zeros(Bool,Nsyms)
    SMALL = 1e-10
    for isym = 1:Nsyms
        non_symmorphic[isym] = ( (abs(ft[1,isym]) >= SMALL) ||
                                 (abs(ft[2,isym]) >= SMALL) ||
                                 (abs(ft[3,isym]) >= SMALL) )
    end

    sr = zeros(Float64,3,3,Nsyms) # s in Cartesian
    RecVecs = inv(Matrix(atoms.LatVecs')) # exclude 2π factor
    sb = zeros(Float64, 3, 3)
    for isym in 1:Nsyms
        @views sb[:,:] = RecVecs * s[:,:,isym]
        @views sr[:,:,isym] = atoms.LatVecs * sb'
    end

    D1, D2, D3 = _init_D_matrices(sr)

    return SymmetryInfo(
        Nrots, Nsyms,
        s[:,:,1:Nsyms], inv_s, idx_inv_s, sname[1:Nsyms],
        sr,
        ft[:,1:Nsyms], non_symmorphic, t_rev,
        irt[1:Nsyms,:],
        D1, D2, D3
    )

end


# From d_matrix.f90 in QE
# sr: rotation matrices in Cartesian
function _init_D_matrices( sr::Array{Float64,3} )
    # !! Calculates the d-matrices.
    # !
    # USE kinds,            ONLY: DP
    # USE symm_base,        ONLY: nsym, sr
    # USE random_numbers,   ONLY: randy
    # USE matrix_inversion

    Nsyms = size(sr, 3)

    # Hardcoded parameters 
    # MAXL = max value of l allowed
    # MAXL = number of m components for l=MAXL
    # MAXLM = number of l,m spherical harmonics for l <= MAXL
    MAXL = 3
    MAXM = 2*MAXL + 1
    MAXLM = (MAXL + 1)^2

    SMALL = 1.0e-9

    # These arrays will be returned
    dy1 = zeros(Float64, 3, 3, Nsyms)
    dy2 = zeros(Float64, 5, 5, Nsyms)
    dy3 = zeros(Float64, 7, 7, Nsyms)

    #
    # randomly distributed points on a sphere
    #
    rl = zeros(Float64, 3, MAXM)
    for m in 1:MAXM
        rl[1,m] = rand() - 0.5
        rl[2,m] = rand() - 0.5
        rl[3,m] = rand() - 0.5
    end
  
    # CALL ylmr2( maxlm, 2*maxl+1, rl, rrl, ylm )
    ylm = zeros(Float64, MAXM, MAXLM)
    Ylm_real_qe!(MAXL, rl, ylm) # Ylm_real_qe accept l value starting from 0


    # invert Yl for each block of definite l (note the transpose operation)
    
    # l = 1 block
    yl1 = zeros(Float64, 3, 3)
    for m in 1:3, n in 1:3
        yl1[m,n] = ylm[n,1+m]
    end
    yl1_inv = inv(yl1)  
    
    #  l = 2 block
    yl2 = zeros(Float64, 5, 5)
    for m in 1:5, n in 1:5
        yl2[m,n] = ylm[n,4+m]
    end
    yl2_inv = inv(yl2)


    #  l = 3 block
    yl3 = zeros(Float64, 7, 7)
    for m in 1:7, n in 1:7
        yl3[m,n] = ylm[n,9+m]
    end
    yl3_inv = inv(yl3)



    # now for each symmetry operation of the point-group

    srl = zeros(Float64, 3, MAXM)
    delt = zeros(Float64, 7, 7)
    ylms = zeros(Float64, MAXM, MAXLM)

    for isym in 1:Nsyms
        #
        # srl[:,m] = rotated rl[:,m] vectors
        #
        @views srl[:,:] = sr[:,:,isym] * rl
        # srl = MATMUL( sr(:,:,isym), rl )
        
        # CALL ylmr2( maxlm, maxm, srl, rrl, ylms )
        Ylm_real_qe!(MAXL, srl, ylms)

        #  find  D_S = Yl_S * Yl_inv (again, beware the transpose)
        
        # l = 1
        for m in 1:3, n in 1:3
            yl1[m,n] = ylms[n,1+m]
        end
        @views dy1[:,:,isym] = yl1[:,:] * yl1_inv[:,:]

        # l = 2 block
        for m in 1:5, n in 1:5
            yl2[m,n] = ylms[n,4+m]
        end
        @views dy2[:,:,isym] = yl2[:,:] * yl2_inv[:,:]


        # l = 3 block
        for m in 1:7, n in 1:7
            yl3[m,n] = ylms[n,9+m]
        end
        @views dy3[:,:,isym] = yl3[:,:] * yl3_inv[:,:]
    end


    # check that D_S matrices are orthogonal as they should if Ylm are correctly defined.
    fill!( delt, 0.0 )
    for m in 1:7
        delt[m,m] = 1.0
    end
    
    for isym in 1:Nsyms
        _check_is_d_orthogonal(isym, dy1, 1, 3, delt, SMALL)
        _check_is_d_orthogonal(isym, dy2, 2, 5, delt, SMALL)
        _check_is_d_orthogonal(isym, dy3, 3, 7, delt, SMALL)
    end
    return dy1, dy2, dy3
end


function _check_is_d_orthogonal(isym, D, l, Nsize, delt, SMALL)
    # l = 3 block
    cc = 0.0
    for m in 1:Nsize, n in 1:Nsize
        dd = dot(D[:,m,isym], D[:,n,isym])
        cc += ( dd - delt[m,n] )^2
    end
    if cc > SMALL
        @printf("D_S (l=%d) for this symmetry operation (isym=%d) is not orthogonal\n", l, isym)
        error("Error in init_d_matrices")
    end
    return
end



function symmetrize_vector!( LatVecs_, sym_info::SymmetryInfo, v::Array{Float64,2} )

    Nsyms = sym_info.Nsyms
    if Nsyms == 1
        return
    end

    alat = norm(LatVecs_[:,1])
    LatVecs = LatVecs_[:,:]/alat

    RecVecs = inv(Matrix(LatVecs'))
    
    s = convert(Array{Float64,3}, sym_info.s)
    
    irt = sym_info.irt

    Nvecs = size(v)[2]

    tmp = zeros(3,Nvecs)
    
    # bring vector to crystal axis
    for i = 1:Nvecs
        tmp[:,i] = v[1,i]*LatVecs[1,:] + v[2,i]*LatVecs[2,:] + v[3,i]*LatVecs[3,:]
    end

    # symmetrize in crystal axis
    v[:,:] .= 0.0
    dv = zeros(3)
    for i = 1:Nvecs
        for isym = 1:Nsyms
            iar = irt[isym,i]
            v[:,i] = v[:,i] + s[:,1,isym]*tmp[1,iar] + s[:,2,isym]*tmp[2,iar] + s[:,3,isym]*tmp[3,iar]
        end
    end
    
    tmp[:,:] = v[:,:]/Nsyms

    # bring vector back to cartesian axis
    for i = 1:Nvecs
        v[:,i] = tmp[1,i]*RecVecs[:,1] + tmp[2,i]*RecVecs[:,2] + tmp[3,i]*RecVecs[:,3]
    end

    return
end


function symmetrize_matrix!( LatVecs_, sym_info::SymmetryInfo, M )
    # only for 3x3 matrix
    # Symmetrize a function \(f(i,j)\) (e.g. stress, dielectric tensor in
    # cartesian axis), where \(i,j\) are the cartesian components.
    
    Nsyms = sym_info.Nsyms
    if Nsyms == 1
        return
    end

    work = zeros(Float64, 3, 3)

    alat = norm(LatVecs_[:,1])
    LatVecs = LatVecs_[:,:]/alat
    RecVecs = inv(Matrix(LatVecs'))

    # bring matrix to crystal axis
    # CALL cart_to_crys( matr )
    #@info "Original M"
    #display(M)

    fill!(work, 0.0)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        work[i,j] += M[k,l] * LatVecs[k,i] * LatVecs[l,j]
    end
    M[:,:] .= work[:,:]
    # XXX: Check this!

    #@info "M"
    #display(M)

    s = convert(Array{Float64,3}, sym_info.s)
    # symmetrize in crystal axis
    fill!(work, 0.0)
    for isym in 1:Nsyms
        for i in 1:3, j in 1:3, k in 1:3, l in 1:3
            work[i,j] += s[i,k,isym] * s[j,l,isym] * M[k,l]
        end
    end
    M[:,:] = work[:,:] / Nsyms
    
    # bring matrix back to cartesian axis
    fill!(work, 0.0)
    for i in 1:3, j in 1:3, k in 1:3, l in 1:3
        work[i,j] += M[k,l] * RecVecs[i,k] * RecVecs[j,l]
    end
    M[:,:] .= work[:,:]

    #@info "After M = "
    #display(M)

    return
end
