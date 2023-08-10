#
# Copyright (C) 2009 Quantum ESPRESSO group
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#

# Essentially radial gradient when ae_flag is true
function radial_gradient_AE!(r, f, gf)

    Nr = size(r, 1)

    @assert size(f, 1) == Nr
    @assert size(gf, 1) == Nr

    #
    # This formula is used in the all-electron case.
    #
    for i in 2:(Nr-1)
        num1 = ( r[i+1] - r[i] )^2 * ( f[i-1] - f[i] ) - ( r[i-1] - r[i] )^2 * ( f[i+1] - f[i] )
        denum1 = ( r[i+1] - r[i] ) * ( r[i-1] - r[i] ) * ( r[i+1] - r[i-1] )
        gf[i] = num/denum 
    end
    gf[Nr] = 0.0
    #
    # The gradient in the first point is a linear interpolation of the
    # gradient at point 2 and 3. 
    #     
    gf[1] = gf[2] + ( gf[3] - gf[2] ) * ( r[1] - r[2] ) / ( r[3] - r[2] )
    return

end