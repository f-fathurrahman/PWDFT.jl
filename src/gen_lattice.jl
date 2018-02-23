function gen_lattice_hexagonal(a; coa=8.0/3.0)
    LL = zeros(3,3)
    LL[1,:] = [1.0, 0.0, 0.0]
    LL[2,:] = [cos(pi/3.0), sin(pi/3.0), 0.0]
    LL[3,:] = [0.0, 0.0, coa]
    return a*LL
end

function gen_lattice_fcc(a)
    LL = zeros(3,3)
    LL[1,:] = [-1.0, 0.0, 1.0]
    LL[2,:] = [ 0.0, 1.0, 1.0]
    LL[3,:] = [-1.0, 1.0, 0.0]
    return 0.5*a*LL
end

function gen_lattice_bcc(a)
    LL = zeros(3,3)
    LL[1,:] = [ 1.0,  1.0, 1.0]
    LL[2,:] = [-1.0,  1.0, 1.0]
    LL[3,:] = [-1.0, -1.0, 1.0]
    return 0.5*a*LL
end
