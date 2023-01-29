function test_O2()
    filename = joinpath(DIR_STRUCTURES, "DATA_G2_mols", "O2.xyz")
    atoms = Atoms(ext_xyz_file=filename)
    pspots = get_default_PsPot_GTH(atoms)

    electrons = Electrons( atoms, pspots )
    println(electrons)

    electrons = Electrons( atoms, pspots, (5,7) )
    println(electrons)

    electrons = Electrons( atoms, pspots, (7,5), Nstates_extra=1 )
    println(electrons)

    return nothing
end

function test_Ni_q18_fcc()
    atoms = Atoms(xyz_string=
        """
        1

        Ni  0.0  0.0  0.0
        """, LatVecs = gen_lattice_fcc(5.0))
    println(atoms)
    
    Nspecies = atoms.Nspecies
    
    pspots = Vector{PsPot_GTH}(undef,Nspecies)
    pspfiles = [joinpath(DIR_PSP_GTH_LDA, "Ni-q18.gth")]
    for isp in 1:Nspecies
        pspots[isp] = PsPot_GTH(pspfiles[isp])
        println(pspots[isp])
    end

    electrons = Electrons( atoms, pspots, Nkpt=14, Nstates_empty=1 )
    println(electrons)

    electrons = Electrons( atoms, pspots, Nkpt=14, Nstates_empty=1, Nspin=2 )
    println(electrons)

    electrons = Electrons( atoms, pspots, Nkpt=14, Nstates=13 )
    println(electrons)

    return nothing
end



@testset "Electrons" begin
    @test test_O2() == nothing
    @test test_Ni_q18_fcc() == nothing
end