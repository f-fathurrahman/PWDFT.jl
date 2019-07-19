function init_Ham_Pt_fcc_PBE()
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))

    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pbe_gth", "Pt-q10.gth")]
    ecutwfc = 30.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], extra_states=4, xcfunc="PBE" )
end

const Etot_Pt_fcc_PBE = -26.2327669856 # ABINIT

@testset "Pt fcc PBE Rhoe mix simple" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc_PBE()
    KS_solve_SCF!( Ham, mix_method="simple", betamix=0.5, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc_PBE atol=1e-5
end

@testset "Pt fcc PBE Rhoe mix anderson" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc_PBE()
    KS_solve_SCF!( Ham, mix_method="anderson", betamix=0.2, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc_PBE atol=1e-5
end

@testset "Pt fcc PBE Rhoe mix pulay" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc_PBE()
    KS_solve_SCF!( Ham, mix_method="pulay", betamix=0.2, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc_PBE atol=1e-5
end

@testset "Pt fcc PBE Rhoe mix ppulay" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc_PBE()
    KS_solve_SCF!( Ham, mix_method="ppulay", betamix=0.2, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc_PBE atol=1e-5
end

@testset "Pt fcc PBE Rhoe mix rpulay" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc_PBE()
    KS_solve_SCF!( Ham, mix_method="rpulay", betamix=0.2, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc_PBE atol=1e-5
end

@testset "Pt fcc PBE Rhoe mix broyden" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc_PBE()
    KS_solve_SCF!( Ham, mix_method="broyden", betamix=0.2, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc_PBE atol=1e-5
end

@testset "Pt fcc PBE potential mix simple" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc_PBE()
    KS_solve_SCF_potmix!( Ham, mix_method="simple", betamix=0.2, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc_PBE atol=1e-5
end

@testset "Pt fcc PBE potential mix linear_adaptive" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc_PBE()
    KS_solve_SCF_potmix!( Ham, mix_method="linear_adaptive", betamix=0.2, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc_PBE atol=1e-5
end

@testset "Pt fcc PBE potential mix broyden" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc_PBE()
    KS_solve_SCF_potmix!( Ham, mix_method="broyden", betamix=0.2, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc_PBE atol=1e-5
end
