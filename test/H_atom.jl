function init_Ham_H_atom()
    
    atoms = Atoms( xyz_string="""
    1

    H  0.0  0.0  0.0
    """, LatVecs=gen_lattice_sc(16.0) )
    
    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth", "H-q1.gth")]
    ecutwfc = 15.0

    return Hamiltonian(atoms, pspfiles, 15.0)
end

const Etot_H_atom = -0.4400340349 # ABINIT

#=
@testset "H atom SCF Rhoe mix simple" begin
    Ham = init_Ham_H_atom()
    KS_solve_SCF!( Ham, betamix=0.7, verbose=false )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom atol=1e-5
end

@testset "H atom SCF Rhoe mix linear adaptive" begin
    Ham = init_Ham_H_atom()
    KS_solve_SCF!( Ham, betamix=0.1, mix_method="linear_adaptive", verbose=false )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom atol=1e-5
end

@testset "H atom SCF Rhoe mix anderson" begin
    Ham = init_Ham_H_atom()
    KS_solve_SCF!( Ham, betamix=0.7, mix_method="anderson", verbose=false )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom atol=1e-5
end

@testset "H atom SCF Rhoe mix pulay" begin
    Ham = init_Ham_H_atom()
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="pulay", verbose=false )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom atol=1e-5
end

@testset "H atom SCF Rhoe mix ppulay" begin
    Ham = init_Ham_H_atom()
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="ppulay", verbose=false )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom atol=1e-5
end

@testset "H atom SCF Rhoe mix rpulay" begin
    Ham = init_Ham_H_atom()
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="rpulay", verbose=false )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom atol=1e-5
end

@testset "H atom SCF Rhoe mix broyden" begin
    Ham = init_Ham_H_atom()
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="broyden", verbose=false )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom atol=1e-5
end
=#


# H atom Emin PCG using Fletcher-Reeves does not pass (early convergence problem)
# might be related to line minimization


@testset "H atom Emin PCG Polak-Ribiere" begin
    Ham = init_Ham_H_atom()
    KS_solve_Emin_PCG!( Ham, i_cg_beta=2, verbose=false )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom atol=1e-5
end

@testset "H atom Emin PCG Hestenes-Stiefeld" begin
    Ham = init_Ham_H_atom()
    KS_solve_Emin_PCG!( Ham, i_cg_beta=3, verbose=false )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom atol=1e-5
end

@testset "H atom Emin PCG Dai-Yuan" begin
    Ham = init_Ham_H_atom()
    KS_solve_Emin_PCG!( Ham, i_cg_beta=4, verbose=false )
    println("")
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom atol=1e-5
end