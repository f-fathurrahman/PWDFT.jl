function init_Ham_Pt_fcc()
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))

    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth", "Pt-q10.gth")]
    ecutwfc = 30.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], extra_states=4 )
end

const Etot_Pt_fcc = -26.2945054360 # ABINIT

@testset "Pt fcc Rhoe mix simple" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc()
    KS_solve_SCF!( Ham, mix_method="simple", betamix=0.5, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc atol=1e-5
end

@testset "Pt fcc Rhoe mix anderson" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc()
    KS_solve_SCF!( Ham, mix_method="anderson", betamix=0.2, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc atol=1e-5
end

@testset "Pt fcc Rhoe mix pulay" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc()
    KS_solve_SCF!( Ham, mix_method="pulay", betamix=0.2, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc atol=1e-5
end

@testset "Pt fcc Rhoe mix ppulay" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc()
    KS_solve_SCF!( Ham, mix_method="ppulay", betamix=0.2, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc atol=1e-5
end

@testset "Pt fcc Rhoe mix rpulay" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc()
    KS_solve_SCF!( Ham, mix_method="rpulay", betamix=0.2, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc atol=1e-5
end

@testset "Pt fcc Rhoe mix broyden" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc()
    KS_solve_SCF!( Ham, mix_method="broyden", betamix=0.2, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc atol=1e-5
end

@testset "Pt fcc potential mix simple" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc()
    KS_solve_SCF_potmix!( Ham, mix_method="simple", betamix=0.2, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc atol=1e-5
end

@testset "Pt fcc potential mix linear_adaptive" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc()
    KS_solve_SCF_potmix!( Ham, mix_method="linear_adaptive", betamix=0.2, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc atol=1e-5
end

@testset "Pt fcc potential mix broyden" begin
    Random.seed!(1234)
    Ham = init_Ham_Pt_fcc()
    KS_solve_SCF_potmix!( Ham, mix_method="broyden", betamix=0.2, use_smearing=true, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Pt_fcc atol=1e-5
end
