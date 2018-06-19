using Printf
using PWDFT

function test_CuSO4()
    atoms = init_atoms_xyz("../structures/CuSO4.xyz")
    println(atoms)
    
    Nspecies = atoms.Nspecies
    
    Pspots = Array{PsPot_GTH}(undef,Nspecies)
    pspfiles = ["../pseudopotentials/pade_gth/Cu-q11.gth",
                "../pseudopotentials/pade_gth/O-q6.gth",
                "../pseudopotentials/pade_gth/S-q6.gth"]
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
    atoms = init_atoms_xyz_string(
        """
        1

        Ni  0.0  0.0  0.0
        """
    )
    atoms.LatVecs = gen_lattice_fcc(5.0)
    println(atoms)
    
    Nspecies = atoms.Nspecies
    
    Pspots = Array{PsPot_GTH}(undef,Nspecies)
    pspfiles = ["../pseudopotentials/pade_gth/Ni-q18.gth"]

    for isp = 1:Nspecies
        Pspots[isp] = PsPot_GTH(pspfiles[isp])
        println(Pspots[isp])
    end

    electrons = Electrons( atoms, Pspots, Nstates_empty=1 )
    println(electrons)

    electrons = Electrons( atoms, Pspots, Nstates_empty=1, Nspin=2 )
    println(electrons)

end

test_CuSO4()
test_Ni_fcc()
