# Based on the file PW/src/symm_base.f90 of Quantum ESPRESS
#
# Copyright (C) 2010-2017 Quantum ESPRESSO group
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#

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

const GBL_s0name = [
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
    sname = [repeat(" ", 45) for i in 1:48]

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

        sname[nrot] = GBL_s0name[irot]
        imat[nrot] = irot
        nrot = nrot + 1

    end

    #println("nrot after loop = ", nrot)

    nrot = nrot - 1

    if !(nrot in [1, 2, 4, 6, 8, 12, 24])
        println("WARNING: Bravais lattice has wrong number of symmetries: ", nrot)
        println("symmetries are disabled")
        nrot = 1
    end

    #display(s[:,:,2]); println()

    # Set the inversion symmetry (Bravais lattices have always inversion symmetry)
    for irot = 1:nrot
        sname[irot+nrot] = GBL_s0name[imat[irot]+32]
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
        println("NOTICE: Symmetry group for Bravais lattice is not a group")
        nrot = 1
    end

    println("nrot = ", nrot)

    return nrot, s, ft, sname

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


# Given the point group of the Bravais lattice, this routine finds
# the subgroup which is the point group of the considered crystal.  
# Non symmorphic groups are allowed, provided that fractional
# translations are allowed (nofrac=.false) and that the unit cell
# is not a supercell.
# On output, the array sym is set to .TRUE.. for each operation
# of the original point group that is also a symmetry operation
# of the crystal symmetry point group.
function sgam_at!( atoms::Atoms, sym, no_z_inv, nrot, s, ft, sname, irt )
    
    nat = atoms.Natoms
    tau = atoms.positions
    atm2species = atoms.atm2species

    bg = inv(Matrix(atoms.LatVecs'))  # follows QE convention (no 2*pi factor)

    xau = zeros(Float64,3,nat)
    rau = zeros(Float64,3,nat)
    ft_ = zeros(Float64,3)

    # ... Compute the coordinates of each atom in the basis of
    # the direct lattice vectors
    for ia = 1:nat
       xau[:,ia] = bg[1,:]*tau[1,ia] + bg[2,:]*tau[2,ia] + bg[3,:]*tau[3,ia]
    end
    
    #println("xau matrix")
    #display(xau'); println()

    #
    # ... check if the identity has fractional translations (this means
    # that the cell is actually a supercell). When this happens, fractional
    # translations are disabled, because there is no guarantee that the 
    # generated sym.ops. form a group.
    #
    nb = 1
    irot = 1

    fractional_translations = true # FIXME true by default

    if fractional_translations
        for na in 2:nat
            if atm2species[nb] == atm2species[na]
                #
                ft_[:] = xau[:,na] - xau[:,nb] - round.( xau[:,na] - xau[:,nb] )

                sym[irot] = checksym!( irot, nat, atm2species, xau, xau, ft_, irt )
                
                if sym[irot]
                    fractional_translations = false
                    println("Found symmetry operation: I + ", ft)
                    println("This is a supercell")
                    println("Fractional translations are disabled")
                    
                    #println("Should GOTO 10")
                    break # go outside the loop over na
                end # if
            
            end # if
        
        end # for
    
    end # if

    #display(irt); println();
    # continue 10
    #println("This is GOTO 10")


    nsym_ns = 0
    fft_fact = [1,1,1]
    ftaux = zeros(Float64,3) 

    for irot in 1:nrot
        
        #println("\nirot = ", irot)

        #println("rau matrix:")
        for na in 1:nat
           # rau = rotated atom coordinates
           rau[:,na] = s[1,:,irot] * xau[1,na] + s[2,:,irot] * xau[2,na] + s[3,:,irot] * xau[3,na]
           #@printf("%18.10f %18.10f %18.10f\n", rau[1,na], rau[2,na], rau[3,na])
        end
        
        # first attempt: no fractional translation
        ft[:,irot] .= 0.0 #
        ft_[:] .= 0.0
        
        sym[irot] = checksym!( irot, nat, atm2species, xau, rau, ft_, irt )
        #if sym[irot]
        #    println("true at the first attempt")
        #end
        
        if ( !sym[irot] && fractional_translations )
            nb = 1
            for na in 1:nat
                if atm2species[nb] == atm2species[na]
                 
                    # ... second attempt: check all possible fractional translations
                    ft_[:] = rau[:,na] - xau[:,nb] - round.( rau[:,na] - xau[:,nb] )
                    #@printf("ft_ = %18.10f %18.10f %18.10f\n", ft_[1], ft_[2], ft_[3])
                    
                    # ... ft_ is in crystal axis and is a valid fractional translation
                    # only if ft_(i)=0 or ft_(i)=1/n, with n=2,3,4,
                    
                    for i in 1:3
                        if abs(ft_[i]) > GBL_eps2
                            ftaux[i] = abs( 1.0/ft_[i] - round(Int64, 1.0/ft_[i]) )
                            nfrac = round( Int64, 1.0/abs(ft_[i]) )
                            if( (ftaux[i] < GBL_eps2) && (nfrac != 2) && (nfrac != 3) && (nfrac != 4) && (nfrac != 6) )
                                ftaux[i] = 2*GBL_eps2
                            end # if
                        else
                            ftaux[i] = 0.0
                        end # if
                    end # for
                    
                    if( any( ftaux .> GBL_eps2 ) )
                        #println("Should be doing cycle here irot ", irot)
                        continue #CYCLE
                    end
                    #println("should not after cycle")
                    
                    sym[irot] = checksym!( irot, nat, atm2species, xau, rau, ft_, irt )
                 
                    if sym[irot]
                        nsym_ns = nsym_ns + 1
                        ft[:,irot] = ft_[:]
                    
                        # ... Find factors that must be present in FFT grid dimensions
                        # in order to ensure that fractional translations are
                        # commensurate with FFT grids.
                        for i in 1:3
                            if abs(ft_[i]) > GBL_eps2
                                nfrac = round( Int64, 1.0/abs(ft_[i]) )
                            else
                                nfrac = 0
                            end # if
                            fft_fact[i] = _mcm( fft_fact[i], nfrac)
                        end # do
                        #println("Should GOTO 20")
                        break # break from loop over na
                    end # if
                end # if

            end # do
           
        end # if
        
        #20   CONTINUE next irot
    end

    #for i in 1:48
    #    if sym[i]
    #        println(i, " ", sym[i], " sname = ", sname[i])
    #    end
    #end
    #println("count(Nsym) = ", count(sym))
    
    println("fft_fact = ", fft_fact)

    # ... disable all symmetries z -> -z
    # skipped

    return

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
        mcm = n2
    end
    return res
end

# This function receives as input all the atomic positions xau,
# and the rotated rau by the symmetry operation ir. It returns
# .TRUE. if, for each atom na, it is possible to find an atom nb
# which is of the same type of na, and coincides with it after the
# symmetry operation. Fractional translations are allowed.
function checksym!( irot, nat, ityp, xau, rau, ft_, irt )

    ACCEPT = 1e-5

    continue_na = zeros(Bool, nat)
    
    for na in 1:nat

        for nb in 1:nat
            if ityp[nb] == ityp[na]

                #@printf("na = %3d nb = %3d\n", na, nb)
                #@printf("rau = %18.10f %18.10f %18.10f\n", rau[1,na], rau[2,na], rau[3,na])
                #@printf("xau = %18.10f %18.10f %18.10f\n", xau[1,na], xau[2,na], xau[3,na])
                #@printf("ft_ = %18.10f %18.10f %18.10f\n", ft_[1], ft_[2], ft_[3])

                is_equal =  eqvect( rau[:,na], xau[:,nb], ft_ , ACCEPT )
                if is_equal
                   #!
                   #! ... the rotated atom does coincide with one of the like atoms
                   #! keep track of which atom the rotated atom coincides with
                   irt[irot, na] = nb
                   
                   #println("Will GOTO 10")
                   continue_na[na] = true
                   break # move outside loop over nb
                end # if
            end # if
        end # for
        
        # ... the rotated atom does not coincide with any of the like atoms
        # s(ir) + ft is not a symmetry operation
        if !continue_na[na]
            #println("eqvect immediately returning false")
            return false
        else
            #10   CONTINUE
            #println("This should be after GOTO 10")
            continue
        end
    end

    # ... s(ir) + ft is a symmetry operation
    return true
end


function eqvect(x, y, f, accep)
  
    # This function test if the difference x-y-f is an integer.
    # x, y = 3d vectors in crystal axis, f = fractional translation
    # adapted from QE-6.4

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
    # INTEGER :: stemp(3,3), ttemp, irot, jrot
    # REAL(DP) :: ft_(3)
    # INTEGER, ALLOCATABLE :: irtemp(:)
    # CHARACTER(LEN=45) :: nametemp
    
    stemp = zeros(Int64,3,3)
    ft_ = zeros(Float64,3)

    irtemp = zeros(Int64, size(irt,2) )
    #ALLOCATE ( irtemp(SIZE(irt,2)) )
    
    jrot = 0
    for irot in 1:nrot_
        if sym[irot]
            jrot = jrot + 1
            if irot > jrot
                stemp[:,:] = s[:,:,jrot]
                s[:,:,jrot] = s[:,:,irot]
                s[:,:,irot] = stemp[:,:]
                
                ft_[:] = ft[:,jrot]
                ft[:,jrot] = ft[:,irot]
                ft[:,irot] = ft_[:]
                
                irtemp[:] = irt[jrot,:]
                irt[jrot,:] = irt[irot,:]
                irt[irot,:] = irtemp[:]
                
                nametemp = sname[jrot]
                sname[jrot] = sname[irot]
                sname[irot] = nametemp
                
                # for collinear magnetism, SKIPPED
                #ttemp = t_rev[jrot]
                #t_rev[jrot] = t_rev[irot]
                #t_rev[irot] = ttemp
            end # if
        end # if
    end # for
    
    sym[1:jrot] .= true
    sym[jrot+1:nrot_] .= false

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

