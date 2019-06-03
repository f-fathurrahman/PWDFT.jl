function calc_forces_finite_diff(
    atoms::Atoms,
    pspfiles::Array{String,1},
    ecutwfc::Float64,
    meshk::Array{Int64,1}
)

    pos_orig = copy(atoms.positions)
    
    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs

    forces = zeros(3,Natoms)
    F_NN = zeros(3,Natoms)
    F_Hartree = zeros(3,Natoms)
    F_Ps_loc = zeros(3,Natoms)
    F_XC = zeros(3,Natoms)
    F_Kinetic = zeros(3,Natoms)
    F_Ps_nloc = zeros(3,Natoms)

    Random.seed!(1234)

    Δ = 0.001
    for ia = 1:Natoms
    for ii = 1:3
        
        # set to original positions
        atoms.positions[:,:] = pos_orig[:,:]

        atoms.positions[ii,ia] = pos_orig[ii,ia] + 0.5*Δ
        println(atoms)
        
        Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk )
        
        KS_solve_SCF!( Ham, mix_method="rpulay", etot_conv_thr=1e-6 )
        Eplus = sum(Ham.energies)
        energies1 = deepcopy(Ham.energies)
        
        atoms.positions[ii,ia] = pos_orig[ii,ia] - 0.5*Δ
        println(atoms)

        Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk )
        #
        KS_solve_SCF!( Ham, mix_method="rpulay", etot_conv_thr=1e-6 )
        Eminus = sum(Ham.energies)
        energies2 = deepcopy(Ham.energies)
        #
        forces[ii,ia] = -(Eplus - Eminus)/Δ

        # in Ry/au, for comparison with PWSCF
        F_Kinetic[ii,ia] = -2*(energies1.Kinetic - energies2.Kinetic)/Δ
        F_NN[ii,ia] = -2*(energies1.NN - energies2.NN)/Δ
        F_Ps_loc[ii,ia] = -2*(energies1.Ps_loc - energies2.Ps_loc)/Δ
        F_Ps_nloc[ii,ia] = -2*(energies1.Ps_nloc - energies2.Ps_nloc)/Δ
        F_Hartree[ii,ia] = -2*(energies1.Hartree - energies2.Hartree)/Δ
        F_XC[ii,ia] = -2*(energies1.XC - energies2.XC)/Δ

        println("")
        @printf("ia = %d, idir = %d\n", ia, ii)
        @printf("Eplus  = %18.10f\n", Eplus)
        @printf("Eminus = %18.10f\n", Eminus)
        @printf("diff   = %18.10f\n", Eplus - Eminus)
        @printf("F      = %18.10f\n", forces[ii,ia])
    end
    end

    println("Kinetic forces:")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Kinetic[1,ia], F_Kinetic[2,ia], F_Kinetic[3,ia] )
    end

    println("XC forces:")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_XC[1,ia], F_XC[2,ia], F_XC[3,ia] )
    end

    println("NN forces:")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_NN[1,ia], F_NN[2,ia], F_NN[3,ia] )
    end

    println("Ps loc forces:")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_loc[1,ia], F_Ps_loc[2,ia], F_Ps_loc[3,ia] )
    end

    println("Ps nloc forces:")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_nloc[1,ia], F_Ps_nloc[2,ia], F_Ps_nloc[3,ia] )
    end

    return forces

end
