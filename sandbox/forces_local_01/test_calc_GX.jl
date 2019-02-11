using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

function create_GaAs_v1()

    LatVecs = zeros(3,3)
    LatVecs[:,1] = [0.5, 0.5, 0.0]
    LatVecs[:,2] = [0.5, 0.0, 0.5]
    LatVecs[:,3] = [0.0, 0.5, 0.5]
    LatVecs = LatVecs*5.6537*ANG2BOHR

    atoms = Atoms(xyz_string_frac=
        """
        2

        Ga  0.0  0.0  0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=LatVecs)

    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]

    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )

    return Ham.atoms, Ham.pw
end


function create_GaAs_v2()
    
    atoms = Atoms(xyz_string_frac=
        """
        2

        Ga  0.0  0.0  0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.6537*ANG2BOHR))

    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]

    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )

    return Ham.atoms, Ham.pw
end


function test_main()
    atoms1, pw1 = create_GaAs_v1()
    atoms2, pw2 = create_GaAs_v2()

    Natoms = atoms1.Natoms
    Ng = pw1.gvec.Ng

    atpos1 = atoms1.positions
    atpos2 = atoms2.positions

    for ia = 1:Natoms
        isp = atoms1.atm2species[ia]        
        for ig = 1:Ng

            GX1 = atpos1[1,ia]*pw1.gvec.G[1,ig] +
                  atpos1[2,ia]*pw1.gvec.G[2,ig] +
                  atpos1[3,ia]*pw1.gvec.G[3,ig]
            Sf1 = cos(GX1) - im*sin(GX1)
            
            GX2 = atpos2[1,ia]*pw2.gvec.G[1,ig] +
                  atpos2[2,ia]*pw2.gvec.G[2,ig] +
                  atpos2[3,ia]*pw2.gvec.G[3,ig]
            Sf2 = cos(GX2) - im*sin(GX2)
            @printf("%18.10f %18.10f diff = %18.10e\n", GX1, GX2, abs(GX1-GX2))
        end
    end

    println("Ng1 = ", pw1.gvec.Ng)
    println("Ng2 = ", pw2.gvec.Ng)

end

test_main()