# Based on PWSCF documentation (version 6.2)

function gen_lattice_cubic( a::Float64 )
    v1 = a*[1,0,0]
    v2 = a*[0,1,0]
    v3 = a*[0,0,1]
    #
    LL = zeros(3,3)
    LL[:,1] = v1
    LL[:,2] = v2
    LL[:,3] = v3
    return LL
end

const gen_lattice_sc = gen_lattice_cubic

function gen_lattice_fcc( a::Float64 )
    v1 = 0.5*a*[-1,0,1]
    v2 = 0.5*a*[0,1,1]
    v3 = 0.5*a*[-1,1,0]
    #
    LL = zeros(3,3)
    LL[:,1] = v1
    LL[:,2] = v2
    LL[:,3] = v3
    return LL
end

function gen_lattice_bcc( a::Float64 )
    v1 = 0.5*a*[1,1,1]
    v2 = 0.5*a*[-1,1,1]
    v3 = 0.5*a*[-1,-1,1]
    #
    LL = zeros(3,3)
    LL[:,1] = v1
    LL[:,2] = v2
    LL[:,3] = v3
    return LL    
end

# more symmetric axis:
function gen_lattice_bcc_v2( a::Float64 )
    v1 = 0.5*a*[-1,1,1]
    v2 = 0.5*a*[1,-1,1]
    v3 = 0.5*a*[1,1,-1]
    #
    LL = zeros(3,3)
    LL[:,1] = v1
    LL[:,2] = v2
    LL[:,3] = v3
    return LL
end

# also for trigonal P
function gen_lattice_hexagonal( a::Float64, c::Float64 )
    v1 = a*[1,0,0]
    v2 = a*[-0.5,sqrt(3)/2,0]
    v3 = [0,0,c]
    #
    LL = zeros(3,3)
    LL[:,1] = v1
    LL[:,2] = v2
    LL[:,3] = v3
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
    LL[:,1] = v1
    LL[:,2] = v2
    LL[:,3] = v3
    return LL
end

gen_lattice_rhombohedral = gen_lattice_trigonal

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
    LL[:,1] = v1
    LL[:,2] = v2
    LL[:,3] = v3
    return LL
end

#Tetragonal P (st)               celldm(3)=c/a
#v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,c/a)
function gen_lattice_tetragonal_P( a::Float64, c::Float64)
    v1 = a*[1,0,0]
    v2 = a*[0,1,0]
    v3 = [0,0,c]
    #
    LL = zeros(3,3)
    LL[:,1] = v1
    LL[:,2] = v2
    LL[:,3] = v3
    return LL
end


#7          Tetragonal I (bct)              celldm(3)=c/a
#v1=(a/2)(1,-1,c/a),  v2=(a/2)(1,1,c/a),  v3=(a/2)(-1,-1,c/a)
function gen_lattice_tetragonal_I( a::Float64, c::Float64 )
    v1 = 0.5*[a,-a,c]
    v2 = 0.5*[a,a,c]
    v3 = 0.5*[-a,-a,c]
    #
    LL = zeros(3,3)
    LL[:,1] = v1
    LL[:,2] = v2
    LL[:,3] = v3
    return LL
end

# 8          Orthorhombic P                  celldm(2)=b/a
#                                              celldm(3)=c/a
#       v1 = (a,0,0),  v2 = (0,b,0), v3 = (0,0,c)
function gen_lattice_orthorhombic( a::Float64, b::Float64, c::Float64 )
    v1 = [a,0,0]
    v2 = [0,b,0]
    v3 = [0,0,c]
    #
    LL = zeros(3,3)
    LL[:,1] = v1
    LL[:,2] = v2
    LL[:,3] = v3
    return LL
end

# 12  Monoclinic P, unique axis c     celldm(2)=b/a
#                                     celldm(3)=c/a,
#                                     celldm(4)=cos(ab)
#       v1=(a,0,0), v2=(b*cos(gamma),b*sin(gamma),0),  v3 = (0,0,c)
#       where gamma is the angle between axis a and b.
function gen_lattice_monoclinic( a::Float64, b::Float64, c::Float64, gamma_degree::Float64 )
    gamma = gamma_degree*pi/180
    v1 = [a,0,0]
    v2 = [b*cos(gamma), b*sin(gamma), 0]
    v3 = [0,0,c]
    #
    LL = zeros(3,3)
    LL[:,1] = v1
    LL[:,2] = v2
    LL[:,3] = v3
    return LL
end

# 14 Triclinic celldm(2)= b/a,
#              celldm(3)= c/a,
#              celldm(4)= cos(bc),
#              celldm(5)= cos(ac),
#              celldm(6)= cos(ab)
#       v1 = (a, 0, 0),
#       v2 = (b*cos(gamma), b*sin(gamma), 0)
#       v3 = (c*cos(beta),  c*(cos(alpha)-cos(beta)cos(gamma))/sin(gamma),
#            c*sqrt( 1 + 2*cos(alpha)cos(beta)cos(gamma)
#                      - cos(alpha)^2-cos(beta)^2-cos(gamma)^2 )/sin(gamma) )
#       where alpha is the angle between axis b and c
#              beta is the angle between axis a and c
#             gamma is the angle between axis a and b

function gen_lattice_triclinic( a::Float64, b::Float64, c::Float64,
                                alpha_degree::Float64, beta_degree::Float64, gamma_degree::Float64)

    if alpha_degree + beta_degree + gamma_degree <= 180.0
        error("sum of angles must be larger than 180°")
    end

    alpha = alpha_degree*pi/180
    beta = beta_degree*pi/180
    gamma = gamma_degree*pi/180

    #
    v1 = [a, 0, 0]
    v2 = [b*cos(gamma), b*sin(gamma), 0]
    t1 = c*(cos(alpha)-cos(beta)cos(gamma))/sin(gamma)
    t2 = c*sqrt( 1 + 2*cos(alpha)cos(beta)cos(gamma) -
                 cos(alpha)^2-cos(beta)^2-cos(gamma)^2 )/sin(gamma)
    v3 = [c*cos(beta), t1, t2]
    #
    LL = zeros(3,3)
    LL[:,1] = v1
    LL[:,2] = v2
    LL[:,3] = v3
    return LL
end
