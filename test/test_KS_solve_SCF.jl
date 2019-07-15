const DIR_PWDFT = joinpath( dirname(pathof(PWDFT)), "..")

Ham = Hamiltonian(    
    Atoms( xyz_string="""
    1

    H  0.0  0.0  0.0
    """, LatVecs=gen_lattice_sc(16.0) ),
    [joinpath(DIR_PWDFT, "pseudopotentials/pade_gth/H-q1.gth")],
    15.0
)

const Etot_H_atom = -0.4400340532

@testset "H atom SCF Rhoe mix simple" begin
    KS_solve_SCF!( Ham, betamix=0.7, verbose=false )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom atol=1e-5
end

@testset "H atom SCF Rhoe mix anderson" begin
    KS_solve_SCF!( Ham, betamix=0.7, mix_method="anderson", verbose=false )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom atol=1e-5
end

@testset "H atom SCF Rhoe mix pulay" begin
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="pulay", verbose=false )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom atol=1e-5
end

@testset "H atom SCF Rhoe mix ppulay" begin
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="ppulay", verbose=false )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom atol=1e-5
end

@testset "H atom SCF Rhoe mix rpulay" begin
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="rpulay", verbose=true )
    Etot = sum(Ham.energies)
    @test Etot ≈ Etot_H_atom atol=1e-5
end