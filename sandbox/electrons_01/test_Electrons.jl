using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

function test_CuSO4()
    atoms = Atoms(
        xyz_file=joinpath(DIR_STRUCTURES,"CuSO4.xyz"),
        LatVecs=gen_lattice_sc(20.0)
    )
    println(atoms)
    
    Nspecies = atoms.Nspecies
    
    Pspots = Array{PsPot_GTH}(undef,Nspecies)
    pspfiles = [joinpath(DIR_PSP, "Cu-q11.gth"),
                joinpath(DIR_PSP, "O-q6.gth"),
                joinpath(DIR_PSP, "S-q6.gth")]
    for isp = 1:Nspecies
        Pspots[isp] = PsPot_GTH(pspfiles[isp])
        println(Pspots[isp])
    end

    electrons = Electrons( atoms, Pspots, Nstates_empty=1 )
    println(electrons)

    electrons = Electrons( atoms, Pspots, Nstates_empty=1, Nspin=2 )
    println(electrons)

end


function test_Ni_fcc()
    atoms = Atoms(xyz_string=
        """
        1

        Ni  0.0  0.0  0.0
        """, LatVecs = gen_lattice_fcc(5.0))
    println(atoms)
    
    Nspecies = atoms.Nspecies
    
    Pspots = Array{PsPot_GTH}(undef,Nspecies)
    pspfiles = [joinpath(DIR_PSP, "Ni-q18.gth")]

    for isp = 1:Nspecies
        Pspots[isp] = PsPot_GTH(pspfiles[isp])
        println(Pspots[isp])
    end

    electrons = Electrons( atoms, Pspots, Nkpt=14, Nstates_empty=1 )
    println(electrons)

    electrons = Electrons( atoms, Pspots, Nkpt=14, Nstates_empty=1, Nspin=2 )
    println(electrons)

end

test_CuSO4()
test_Ni_fcc()
