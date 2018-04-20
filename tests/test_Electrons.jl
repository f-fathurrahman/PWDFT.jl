using PWDFT

function test_main()
    atoms = init_atoms_xyz("../structures/CuSO4.xyz")
    println(atoms)
    
    Nspecies = atoms.Nspecies
    
    Pspots = Array{PsPot_GTH}(Nspecies)
    pspfiles = ["../pseudopotentials/pade_gth/Cu-q11.gth",
                "../pseudopotentials/pade_gth/O-q6.gth",
                "../pseudopotentials/pade_gth/S-q6.gth"]
    for isp = 1:Nspecies
        Pspots[isp] = PsPot_GTH(pspfiles[isp])
        println(Pspots[isp])
    end

    electrons = Electrons( atoms, Pspots )
    println(electrons)

end

test_main()
