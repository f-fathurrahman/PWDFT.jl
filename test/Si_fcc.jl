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

#=
@testset "Si fcc SCF Rhoe mix simple" begin
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF!( Ham, betamix=0.5, verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

@testset "Si fcc SCF Rhoe mix linear_adaptive" begin
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF!( Ham, betamix=0.1, mix_method="linear_adaptive", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

@testset "Si fcc SCF Rhoe mix anderson" begin
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="anderson", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end


@testset "Si fcc SCF Rhoe mix pulay" begin
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="pulay", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end


@testset "Si fcc SCF Rhoe mix ppulay" begin
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="ppulay", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

@testset "Si fcc SCF Rhoe mix rpulay" begin
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="rpulay", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

@testset "Si fcc SCF Rhoe mix broyden" begin
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="broyden", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end
=#

@testset "Si fcc SCF potential mix broyden" begin
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF_potmix!( Ham, betamix=0.5, mix_method="simple", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

@testset "Si fcc SCF potential mix broyden" begin
    Ham = init_Ham_Si_fcc()
    KS_solve_SCF_potmix!( Ham, betamix=0.5, mix_method="broyden", verbose=true )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

#=

# Test fail for i_cg_beta=1

@testset "Si fcc Emin PCG Polak-Ribiere" begin
    Ham = init_Ham_Si_fcc()
    KS_solve_Emin_PCG!( Ham, i_cg_beta=2 )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

@testset "H atom Emin PCG Hestenes-Stiefeld" begin
    Ham = init_Ham_Si_fcc()
    KS_solve_Emin_PCG!( Ham, i_cg_beta=3 )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_Si_fcc atol=1e-5
end

=#

# Probably similar problem with i_cg_beta=1
# i_cg_beta=4 is skipped
