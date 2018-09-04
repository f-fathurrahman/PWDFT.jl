using Printf
using LinearAlgebra

using PWDFT

function gen_Rhoe_aux(
    eta::Float64, atoms::Atoms, Zvals::Array{Float64,1}, pw::PWGrid
)
    #
    Npoints = prod(pw.Ns)
    Rhoe_aux = zeros(Float64,Npoints)
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    #
    for ip = 1:Npoints
        r = pw.r[:,ip]
        for ia = 1:Natoms
            isp = atm2species[ia]
            R = atoms.positions[:,ia]
            Z = Zvals[isp]
            dr2 = dot(r-R,r-R)
            #Rhoe_aux[ip] = Rhoe_aux[ip] + Z*exp(-0.5*eta^2*dr2)
            Rhoe_aux[ip] = Rhoe_aux[ip] + Z*exp(-2*eta^2*dr2)
        end
    end
    #return -eta^3/((2*pi)^1.5)*Rhoe_aux
    #return -(2*eta)^3/((2*pi)^1.5)*Rhoe_aux
    return -2*sqrt(2)*eta^3/(pi^1.5)*Rhoe_aux
    #return -eta^3/(pi^1.5)*Rhoe_aux*4/sqrt(2)
    #return Rhoe_aux*(-eta^3/((2*pi)^1.5))*(2^3)
    #return Rhoe_aux
end


function gen_Rhoe_aux_G(
    eta::Float64, atoms::Atoms, Zvals::Array{Float64,1}, pw::PWGrid
)
    Sf = calc_strfact( atoms, pw )
    
    Ng = pw.gvec.Ng
    G2 = pw.gvec.G2
    Nspecies = atoms.Nspecies
    Npoints = prod(pw.Ns)
    idx_g2r = pw.gvec.idx_g2r
    CellVolume = pw.CellVolume
    
    Rhoe_aux_G = zeros(ComplexF64,Npoints)
    Rhoe_aux = zeros(Float64,Npoints)

    for ig = 1:Ng
        ip = idx_g2r[ig]
        for isp = 1:Nspecies
            Rhoe_aux_G[ip] = Rhoe_aux_G[ip] + Zvals[isp]*exp(-0.125*G2[ig]/eta^2)*Sf[ig,isp]/CellVolume
        end
    end
    Rhoe_aux = real( G_to_R(pw,Rhoe_aux_G) )
    #return -eta^3/((2*pi)^1.5)*Rhoe_aux*Npoints/CellVolume
    #return -eta^3/((2*pi)^1.5)*Rhoe_aux*CellVolume
    #return -eta^3/((2*pi)^1.5)*Rhoe_aux*CellVolume
    #return -Rhoe_aux*(pi/2/eta^2)^(3/2)/CellVolume*Npoints*(2*eta)^3/((2*pi)^(3/2))
    return -Rhoe_aux*Npoints
end


function test_H()
    
    atoms = Atoms(xyz_string=
        """
        1

        H  8.0  8.0  8.0
        """, in_bohr=true,
        LatVecs=gen_lattice_sc(16.0))
    
    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )
    println(Ham)

    pw = Ham.pw
    Npoints = prod(pw.Ns)
    dVol = pw.CellVolume/Npoints

    println("dVol = ", dVol)
    
    Zvals = get_Zvals(Ham.pspots) # need to call this ?
    println("Zvals = ", Zvals)
    println("Q = ", sum(Zvals))

    Gcut = 2*ecutwfc/(2*pi)
    TOL = 1e-8
    eta = Gcut^2/-log(TOL)/2
    println("eta = ", eta)

    Rhoe_aux = gen_Rhoe_aux(eta, atoms, Zvals, pw)
    println("R space: integ Rhoe_aux = ", sum(Rhoe_aux)*dVol)

    Rhoe_aux = gen_Rhoe_aux_G(eta, atoms, Zvals, pw)
    println("G space: integ Rhoe_aux = ", sum(Rhoe_aux)*dVol)
end

function test_Si()
    atoms = init_atoms_xyz_string(
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(10.2631)
    atoms.positions = atoms.LatVecs*atoms.positions
    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Si-q4.gth"]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
    #
    pw = Ham.pw
    Npoints = prod(pw.Ns)
    dVol = pw.CellVolume/Npoints

    println("dVol = ", dVol)
    
    Zvals = get_Zvals(Ham.pspots) # need to call this ?
    println("Zvals = ", Zvals)
    println("Q = ", sum(Zvals))

    Gcut = 2*ecutwfc/(2*pi)
    TOL = 1e-8
    eta = Gcut^2/-log(TOL)/2
    println("eta = ", eta)

    Rhoe_aux = gen_Rhoe_aux_G(eta, atoms, Zvals, pw)
    println("G space: integ Rhoe_aux = ", sum(Rhoe_aux)*dVol)
end

test_H()
#test_Si()
