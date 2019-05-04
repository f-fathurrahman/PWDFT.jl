using Printf
using LinearAlgebra
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

function main()
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], Ns_=(32,32,32) )
    println(Ham)

    psiks = rand_BlochWavefunc(Ham)
    Rhoe = calc_rhoe(Ham, psiks)

    pw = Ham.pw
    Npoints = prod(Ham.pw.Ns)
    CellVolume = Ham.pw.CellVolume
    dVol = CellVolume/Npoints

    println("Npoints = ", Npoints)
    println("dVol    = ", dVol)
    println("dVol^-1 = ", 1/dVol)
    println(sum(Rhoe)*dVol)

    RhoeG = R_to_G(Ham.pw, Rhoe[:,1])
    println(RhoeG[1]*dVol)

    Rhoe[:,1] .= 1.0
    Rhoe[4,1] = 2.0
    Rhoe[3,1] = 4.0

    phiG = Poisson_solve( pw, Rhoe[:,1] )
    phi = real( G_to_R(pw, phiG) )
    Ehartree = 0.5*dot( phi, Rhoe[:,1] )*dVol

    println(R_to_G(pw, Rhoe[:,1])[1]*dVol)
    println("Ehartree (in Ry) = ", Ehartree/2) # convert to Rydberg

    Rhoe[:,1] .= 1.0
    Rhoe[4,1] = 2.0
    Rhoe[3,1] = 4.0
    cRhoeG = R_to_G(pw, Rhoe[:,1])
    Ehartree = 0.0
    for ig = 2:Ham.pw.gvec.Ng
        ip = Ham.pw.gvec.idx_g2r[ig]
        Ehartree = Ehartree + abs(cRhoeG[ip])^2/Ham.pw.gvec.G2[ig]        
    end
    Ehartree = 0.5*Ehartree*dVol
    println("Ehartree (in Ry) v2 = ", Ehartree/2)
end

main()
