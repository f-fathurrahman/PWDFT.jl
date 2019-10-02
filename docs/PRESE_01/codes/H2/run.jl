using Random
using PWDFT

function main()
    
    Random.seed!(1234) # to ensure reproducibility
    
    atoms = Atoms( xyz_file="H2.xyz",
                   LatVecs = gen_lattice_sc(16.0) )
    
    write_xsf("H2.xsf", atoms)
    
    pspfiles = ["H-q1.gth"]
    
    ecutwfc = 15.0
    
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )

    #KS_solve_SCF!( Ham, betamix=0.5, mix_method="rpulay" )
    #KS_solve_Emin_PCG!( Ham )
    KS_solve_SCF_potmix!( Ham, betamix=0.5 )
end

@time main()
@time main()