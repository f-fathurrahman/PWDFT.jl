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