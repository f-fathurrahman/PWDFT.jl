function init_Ham_Fe_bcc_PBE() 
    atoms = Atoms(xyz_string_frac=
        """
        1

        Fe  0.0  0.0  0.0
        """, in_bohr=true, LatVecs=gen_lattice_bcc(2.856*ANG2BOHR))

    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pbe_gth", "Fe-q16.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], xcfunc="PBE",
                        Nspin=2, extra_states=4 )
end

const Etot_Fe_bcc = -116.5421659437 # ABINIT

@testset "Fe bcc GGA spinpol rpulay" begin
    Ham = init_Ham_Fe_bcc_PBE()
    KS_solve_SCF!( Ham, mix_method="rpulay", betamix=0.1, use_smearing=true,
                   starting_magnetization=[0.5], verbose=true )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Fe_bcc atol=1e-5
end

@testset "Fe bcc GGA spinpol anderson" begin
    Ham = init_Ham_Fe_bcc_PBE()
    KS_solve_SCF!( Ham, mix_method="anderson", betamix=0.1, use_smearing=true,
                   starting_magnetization=[0.5], verbose=true )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Fe_bcc atol=1e-5
end


@testset "Fe bcc GGA spinpol potmix simple" begin
    Ham = init_Ham_Fe_bcc_PBE()
    KS_solve_SCF_potmix!( Ham, mix_method="simple", use_smearing=true,
                          starting_magnetization=[0.5], verbose=true )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Fe_bcc atol=1e-5
end

