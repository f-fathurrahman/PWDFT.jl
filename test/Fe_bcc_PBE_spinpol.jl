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

const Etot_Fe_bcc_PBE_spinpol = -116.5421659437 # ABINIT

@testset "Fe bcc PBE spinpol Rhoe mix simple" begin
    Random.seed!(1234)
    Ham = init_Ham_Fe_bcc_PBE()
    KS_solve_SCF!( Ham, mix_method="simple", betamix=0.5, use_smearing=true,
                   starting_magnetization=[0.5], verbose=true )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Fe_bcc_PBE_spinpol atol=1e-5
end

@testset "Fe bcc PBE spinpol Rhoe mix anderson" begin
    Random.seed!(1234)
    Ham = init_Ham_Fe_bcc_PBE()
    KS_solve_SCF!( Ham, mix_method="anderson", betamix=0.2, use_smearing=true,
                   starting_magnetization=[0.5], verbose=true )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Fe_bcc_PBE_spinpol atol=1e-5
end

@testset "Fe bcc PBE spinpol Rhoe mix pulay" begin
    Random.seed!(1234)
    Ham = init_Ham_Fe_bcc_PBE()
    KS_solve_SCF!( Ham, mix_method="pulay", betamix=0.2, use_smearing=true,
                   starting_magnetization=[0.5], verbose=true )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Fe_bcc_PBE_spinpol atol=1e-5
end

@testset "Fe bcc PBE spinpol Rhoe mix ppulay" begin
    Random.seed!(1234)
    Ham = init_Ham_Fe_bcc_PBE()
    KS_solve_SCF!( Ham, mix_method="ppulay", betamix=0.2, use_smearing=true,
                   starting_magnetization=[0.5], verbose=true )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Fe_bcc_PBE_spinpol atol=1e-5
end

@testset "Fe bcc PBE spinpol Rhoe mix rpulay" begin
    Random.seed!(1234)
    Ham = init_Ham_Fe_bcc_PBE()
    KS_solve_SCF!( Ham, mix_method="rpulay", betamix=0.2, use_smearing=true,
                   starting_magnetization=[0.5], verbose=true )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Fe_bcc_PBE_spinpol atol=1e-5
end

@testset "Fe bcc PBE spinpol Rhoe mix broyden" begin
    Random.seed!(1234)
    Ham = init_Ham_Fe_bcc_PBE()
    KS_solve_SCF!( Ham, mix_method="broyden", betamix=0.2, use_smearing=true,
                   starting_magnetization=[0.5], verbose=true )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Fe_bcc_PBE_spinpol atol=1e-5
end

@testset "Fe bcc PBE spinpol potential mix simple" begin
    Random.seed!(1234)
    Ham = init_Ham_Fe_bcc_PBE()
    KS_solve_SCF_potmix!( Ham, mix_method="simple", betamix=0.2, use_smearing=true,
                   starting_magnetization=[0.5], verbose=true )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Fe_bcc_PBE_spinpol atol=1e-5
end

@testset "Fe bcc PBE spinpol potential mix linear_adaptive" begin
    Random.seed!(1234)
    Ham = init_Ham_Fe_bcc_PBE()
    KS_solve_SCF_potmix!( Ham, mix_method="linear_adaptive", betamix=0.2, use_smearing=true,
                   starting_magnetization=[0.5], verbose=true )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Fe_bcc_PBE_spinpol atol=1e-5
end

@testset "Fe bcc PBE spinpol potential mix broyden" begin
    Random.seed!(1234)
    Ham = init_Ham_Fe_bcc_PBE()
    KS_solve_SCF_potmix!( Ham, mix_method="broyden", betamix=0.2, use_smearing=true,
                   starting_magnetization=[0.5], verbose=true )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Fe_bcc_PBE_spinpol atol=1e-5
end
