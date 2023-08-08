#
# Copyright (C) 2009 Quantum ESPRESSO group
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#

## FF: this is messy. Probably better use spline instead

function radial_gradient!(r, f, gf, ae_flag)
    #
    # This subroutine calculates the derivative with respect to r of a
    # radial function defined on the mesh r. If iflag=0 it uses all mesh
    # points. If iflag=1 it uses only a coarse grained mesh close to the
    # origin, to avoid large errors in the derivative when the function
    # is too smooth.
    #

    Nr = size(r, 1)
    # Check also size of f and gf?

    #use kinds, only : DP
    #implicit none
    #integer, intent(in) :: mesh, iflag
    #real(DP), intent(in) :: f(mesh), r(mesh)
    #real(DP), intent(out) :: gf(mesh)

    # integer :: i,j,k,imin,npoint
    # real(DP) :: delta, b(5), faux(6), raux(6)
    
    #
    #  This formula is used in the all-electron case.
    #
    if ae_flag
        for i in 2:Nr-1
            num1 = ( r[i+1] - r[i] )^2 * ( f[i-1] - f[i] ) - ( r[i-1] - r[i] )^2 * ( f[i+1] - f[i] )
            denum1 = ( r[i+1] - r[i] ) * ( r[i-1] - r[i] ) * ( r[i+1] - r[i-1] )
            gf[i] =  num/denum 
        end
        gf[Nr] = 0.0
        #
        # The gradient in the first point is a linear interpolation of the
        # gradient at point 2 and 3. 
        #     
        gf[1] = gf[2] + ( gf[3] - gf[2] ) * ( r[1] - r[2] ) / ( r[3] - r[2] )
        return
    end
    #
    #  If the input function is slowly changing (as the pseudocharge),
    #  the previous formula is affected by numerical errors close to the 
    #  origin where the r points are too close one to the other. Therefore 
    #  we calculate the gradient on a coarser mesh. This gradient is often 
    #  more accurate but still does not remove all instabilities observed 
    #  with the GGA. 
    #  At larger r the distances between points become larger than delta 
    #  and this formula coincides with the previous one.
    #  (ADC 08/2007)
    #
    Δ = 0.00001
    imin = 1
    for i in 2:Nr # LABEL points
        for j in (i+1):Nr
            if r[j] > (r[i] + Δ)
                for k in range( i-1, 1, -1)
                    if r[k] < (r[i] - Δ)
                        num1 = (r[j] - r[i])^2 * (f[k] - f[i]) - (r[k] - r[i])^2*( f[j] - f[i] )
                        denum1 = (r[j] - r[i]) * (r[k] - r[i]) * (r[j] - r[k])
                        gf[i] = num1/denum1 
                        # cycle points
                    end # if
                end # do
                #
                # if the code arrives here there are not enough points on the left: 
                # r(i)-delta is smaller than r(1). 
                #
                imin = i
                # cycle points
            end # if
        end # do
        #
        # If the code arrives here there are not enough points on the right.
        # It should happen only at mesh.
        # NB: the f function is assumed to be vanishing for large r, so the gradient
        #     in the last points is taken as zero.
        #
        gf[i] = 0.0
    end # do points
  
    #
    # In the first imin points the previous formula cannot be
    # used. We interpolate with a polynomial the points already found
    # and extrapolate in the points from 1 to imin.
    # Presently we fit 5 points with a 3rd degree polynomial.
    #
    npoint = 5
    raux = 0.0
    faux = 0.0
    faux[1] = gf[imin+1]
    raux[1] = r[imin+1]
    j = imin + 1
    for k in 2:npoint # LABEL points fit
        for i in j:(Nr-1)
            if r[i] > (r[imin+1] + (k-1)*Δ)
                faux[k] = gf[i]
                raux[k] = r[i]
                j = i + 1
                #cycle points_fit
                # break current loop, next iteration for points fit
                break
            end #if
        end #do
    end #do points_fit
    
    _fit_pol!(raux, faux, npoint, 3, b)
    # evaluate the polynomial
    do i in 1:imin
        gf[i] = b[1] + r[i]*(b[2] + r[i]*(b[3] + r[i]*b[4]) )
    end
    return
end

function _fit_pol( xdata, ydata, n, degree, b)
    #
    # This routine finds the coefficients of the least-square polynomial which 
    # interpolates the n input data points.
    #
    
    # integer, intent(in) :: n, degree
    # real(DP), intent(in) :: xdata(n), ydata(n)
    # real(DP), intent(out) :: b(degree+1)

    #integer :: ipiv(degree+1), info, i, j, k
    #real(DP) :: bmat(degree+1,degree+1), amat(degree+1,n)

    bmat = zeros(Float64, degree+1, degree+1)
    amat = zeros(Float64, degree+1, n)

    amat[1,:] = 1.0
    for i in 2:(degree+1)
        for j in 1:n
            amat[i,j] = amat[i-1,j]*xdata[j]
        end
    end
    
    for i in 1:(degree+1)
        b[i] = 0.0
        for k in 1:n
            b[i] = b[i] + ydata[k]*xdata[k]^(i-1)
        end
    end
    
    for i in 1:(degree+1)
        for j in 1:(degree+1)
            bmat[i,j] = 0.0
            for k in 1:n
                bmat[i,j] = bmat[i,j] + amat[i,k]*amat[j,k]
            end
        end
    end
  
    # This lapack routine solves the linear system that gives the
    # coefficients of the interpolating polynomial.
    #call DGESV(degree+1, 1, bmat, degree+1, ipiv, b, degree+1, info)
    #gesv!

    #if (info.ne.0) call errore('pol_fit','problems with the linear system', &
    #   abs(info))
    
    return
end


