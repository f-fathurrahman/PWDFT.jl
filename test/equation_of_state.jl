using PWDFT

include("data_crystal.jl")
include("all_pade_psp.jl")
include("all_pbe_psp.jl")

function test_main()
    atsymb = "Ag"

    a = FCC_DATA[atsymb]*ANG2BOHR

    scal0 = 0.95
    Ndata = 6
    Δ = 0.02
    
    ecutwfc = 15.0
    pspfiles = [ "../pseudopotentials/pade_gth/"*ALL_PADE_PSP[atsymb][1] ]
    #pspfiles = [ "../pseudopotentials/pbe_gth/"*ALL_PBE_PSP[atsymb][1] ]

    for i = 1:Ndata

        scal = scal0 + (i-1)*Δ

        println("Lattice parameter = ", a*scal)
    
        atoms = Atoms( xyz_string_frac="""
        1

        $atsymb 0.0 0.0 0.0
        """, LatVecs=gen_lattice_fcc(a*scal) )
    
        Ham = Hamiltonian(
           atoms, pspfiles, ecutwfc, meshk=[5,5,5], extra_states=4,
        )
        println(Ham)

        @time KS_solve_SCF!( Ham, use_smearing=true, mix_method="rpulay" )
    end

end

test_main()


