function test_ElecVars( Ham::Hamiltonian )
    evars = ElecVars(Ham)
    println(evars)

    E_fermi, mTS = update_occ!( Ham, evars, 0.01 )
    
    println("E_fermi = ", E_fermi)

    println(Ham.electrons)

    println(Ham.electrons.ebands[:,1])

    println("Pass here in test_ElecVars")
end