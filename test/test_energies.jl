# XXX These tests are assuming that Random.seed! is stable accross Julia versions

function test_energy_01()
    Ham = create_Ham_Si_fcc_oncv();
    
    Random.seed!(1234);
    psiks = rand_BlochWavefunc(Ham);
    calc_rhoe!(Ham, psiks, Ham.rhoe);
    update_from_rhoe!(Ham, psiks, Ham.rhoe);
    energies = calc_energies(Ham, psiks);

    @test energies.Kinetic ≈ 95.5765255700
    @test energies.Ps_loc ≈ -4.6284911767
    @test energies.Ps_nloc ≈ 0.8114754475
    @test energies.Hartree ≈ 1.1818870819
    @test energies.XC ≈ -3.8283308065

    E_NN = calc_E_NN(Ham.atoms)
    @test E_NN ≈ -8.397927400714764
    return
end


