# Based on PWSCF documentation (version 6.2)

function gen_lattice_cubic( a::Float64 )
    v1 = a*[1,0,0]
    v2 = a*[0,1,0]
    v3 = a*[0,0,1]
    #
    LL = zeros(3,3)
    LL[1,:] = v1
    LL[2,:] = v2
    LL[3,:] = v3
    return LL
end

function gen_lattice_fcc( a::Float64 )
    v1 = 0.5*a*[-1,0,1]
    v2 = 0.5*a*[0,1,1]
    v3 = 0.5*a*[-1,1,0]
    #
    LL = zeros(3,3)
    LL[1,:] = v1
    LL[2,:] = v2
    LL[3,:] = v3
    return LL
end

function gen_lattice_bcc( a::Float64 )
    v1 = 0.5*a*[1,1,1]
    v2 = 0.5*a*[-1,1,1]
    v3 = 0.5*a*[-1,-1,1]
    #
    LL = zeros(3,3)
    LL[1,:] = v1
    LL[2,:] = v2
    LL[3,:] = v3
    return LL    
end

# more symmetric axis:
function gen_lattice_bcc_v2( a::Float64 )
    v1 = 0.5*a*[-1,1,1]
    v2 = 0.5*a*[1,-1,1]
    v3 = 0.5*a*[1,1,-1]
    #
    LL = zeros(3,3)
    LL[1,:] = v1
    LL[2,:] = v2
    LL[3,:] = v3
    return LL
end

# also for trigonal P
function gen_lattice_hexagonal( a::Float64, c::Float64 )
    v1 = a*[1,0,0]
    v2 = a*[-1/2,sqrt(3)/2,0]
    v3 = [0,0,c]
    #
    LL = zeros(3,3)
    LL[1,:] = v1
    LL[2,:] = v2
    LL[3,:] = v3
    return LL    
end


# 5          Trigonal R, 3fold axis c        celldm(4)=cos(gamma)
# The crystallographic vectors form a three-fold star around
# the z-axis, the primitive cell is a simple rhombohedron:
# v1 = a(tx,-ty,tz),   v2 = a(0,2ty,tz),   v3 = a(-tx,-ty,tz)
# where c=cos(gamma) is the cosine of the angle gamma between
# any pair of crystallographic vectors, tx, ty, tz are:
#   tx=sqrt((1-c)/2), ty=sqrt((1-c)/6), tz=sqrt((1+2c)/3)
function gen_lattice_trigonal( a::Float64, gamma_degree::Float64 )
    c = cos( gamma_degree*pi/180 )
    tx = sqrt((1-c)/2)
    ty = sqrt((1-c)/6)
    tz = sqrt((1+2c)/3)
    #
    v1 = a*[tx,-ty,tz]
    v2 = a*[0,2*ty,tz]
    v3 = a*[-tx,-ty,tz]
    #
    LL = zeros(3,3)
    LL[1,:] = v1
    LL[2,:] = v2
    LL[3,:] = v3
    return LL
end

#   -5          Trigonal R, 3fold axis <111>    celldm(4)=cos(gamma)
# The crystallographic vectors form a three-fold star around
# <111>. Defining a' = a/sqrt(3) :
# v1 = a' (u,v,v),   v2 = a' (v,u,v),   v3 = a' (v,v,u)
# where u and v are defined as
#   u = tz - 2*sqrt(2)*ty,  v = tz + sqrt(2)*ty
# and tx, ty, tz as for case ibrav=5
# Note: if you prefer x,y,z as axis in the cubic limit,
#       set  u = tz + 2*sqrt(2)*ty,  v = tz - sqrt(2)*ty
#       See also the note in Modules/latgen.f90
function gen_lattice_trigonal_v2( a::Float64, gamma_degree::Float64 )
    c = cos( gamma_degree*pi/180 )    
    tx = sqrt((1-c)/2)
    ty = sqrt((1-c)/6)
    tz = sqrt((1+2c)/3)
    u = tz - 2*sqrt(2)*ty
    v = tz + sqrt(2)*ty
    ap = a/sqrt(3)
    #
    v1 = ap*[u,v,v]
    v2 = ap*[v,u,v]
    v3 = ap*[v,v,u]
    #
    LL = zeros(3,3)
    LL[1,:] = v1
    LL[2,:] = v2
    LL[3,:] = v3
    return LL    
end
