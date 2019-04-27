using Printf
using LinearAlgebra
using PWDFT

function main()

    atoms = Atoms( xyz_string="""
    2

    H  0.0   0.0  0.0
    H  1.75  0.0  0.0
    """, in_bohr=true )

    pw = PWGrid( 30.0, gen_lattice_sc(16.0) )
    #println(pw)

    Sf = calc_strfact(atoms, pw)

    Npoints = prod(pw.Ns)
    dVol = pw.CellVolume/Npoints

    Rhoe = zeros(Npoints)

    sigma1 = 0.25
    center = [8.0, 8.0, 8.0]
    Zval = 1.0

    for ip = 1:Npoints
        dr2 = norm( pw.r[:,ip] - center )^2
        Rhoe[ip] = Zval*exp( -dr2/(2*sigma1^2) ) / sqrt(2*pi*sigma1^2)^3
    end

    @printf("integ Rhoe = %18.10f\n", sum(Rhoe)*dVol)

    # This is not strictly needed for this particular case
    RhoeG = R_to_G(pw, Rhoe)
    for ig = 1:pw.gvec.Ng
        ip = pw.gvec.idx_g2r[ig]
        RhoeG[ip] = RhoeG[ip]*Sf[ig,1]
    end

    RhoeR = real( G_to_R(pw, RhoeG) )
    @printf("integ RhoeR = %18.10f\n", sum(RhoeR)*dVol)
    #


    V_Hartree = real( G_to_R(pw, Poisson_solve(pw, RhoeR) ) )
    E_Hartree = 0.5*sum( V_Hartree .* RhoeR)*dVol
    println("E_Hartree = ", E_Hartree)

    E_self = Zval^2/(2*sqrt(pi)*sigma1)*atoms.Natoms
    println("E_self    = ", E_self)

    E_ewald = E_Hartree - E_self
    println("E_ewald   = ", E_ewald)

end

@time main()
@time main()