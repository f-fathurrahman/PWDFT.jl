function init_Ham_Si_fcc()
    atoms = Atoms(xyz_string_frac=
    """
    2

    Si  0.0   0.0   0.0
    Si  0.25  0.25  0.25
    """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth", "Si-q4.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
end

const Etot_Si_fcc = -7.9110087934 # ABINIT


@testset "Si fcc SCF Rhoe mix simple" begin
    Random.seed!(1234)
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF!( Ham, betamix=0.5, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

@testset "Si fcc SCF Rhoe mix linear_adaptive" begin
    Random.seed!(1234)
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF!( Ham, betamix=0.1, mix_method="linear_adaptive", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

@testset "Si fcc SCF Rhoe mix anderson" begin
    Random.seed!(1234)
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="anderson", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end


@testset "Si fcc SCF Rhoe mix pulay" begin
    Random.seed!(1234)
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="pulay", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end


@testset "Si fcc SCF Rhoe mix ppulay" begin
    Random.seed!(1234)
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="ppulay", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

@testset "Si fcc SCF Rhoe mix rpulay" begin
    Random.seed!(1234)
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="rpulay", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

@testset "Si fcc SCF Rhoe mix broyden" begin
    Random.seed!(1234)
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="broyden", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end


@testset "Si fcc SCF potential mix simple" begin
    Random.seed!(1234)
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF_potmix!( Ham, betamix=0.5, mix_method="simple", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

@testset "Si fcc SCF potential mix linear_adaptive" begin
    Random.seed!(1234)
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF_potmix!( Ham, betamix=0.1, mix_method="linear_adaptive", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

@testset "Si fcc SCF potential mix broyden" begin
    Random.seed!(1234)
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF_potmix!( Ham, betamix=0.5, mix_method="broyden", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end


@testset "H atom Emin PCG Fletcher-Reeves" begin
    Random.seed!(1234)
    Ham = init_Ham_Si_fcc()
    KS_solve_Emin_PCG!( Ham, i_cg_beta=1, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

@testset "Si fcc Emin PCG Polak-Ribiere" begin
    Random.seed!(1234)
    Ham = init_Ham_Si_fcc()
    KS_solve_Emin_PCG!( Ham, i_cg_beta=2, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

@testset "H atom Emin PCG Hestenes-Stiefeld" begin
    Random.seed!(1234)
    Ham = init_Ham_Si_fcc()
    KS_solve_Emin_PCG!( Ham, i_cg_beta=3, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

@testset "H atom Emin PCG Dai-Yuan" begin
    Random.seed!(1234)
    Ham = init_Ham_Si_fcc()
    KS_solve_Emin_PCG!( Ham, i_cg_beta=4, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end