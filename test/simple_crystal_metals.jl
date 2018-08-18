using PWDFT

include("data_crystal.jl")
include("all_pade_psp.jl")

function test_bcc()
    for atsymb in keys(BCC_DATA)
        atoms = Atoms( xyz_string_frac="""
        1

        $atsymb 0.0 0.0 0.0
        """, LatVecs=gen_lattice_bcc( BCC_DATA[atsymb]*ANG2BOHR ) )
        pspfiles = [ "../pseudopotentials/pade_gth/"*ALL_PADE_PSP[atsymb][1] ]
        ecutwfc = 15.0

        Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], extra_states=4 )
        println(Ham)

        @time KS_solve_SCF!( Ham, use_smearing=true, mix_method="rpulay" )
    end
end

function test_fcc()
    for atsymb in keys(FCC_DATA)
        atoms = Atoms( xyz_string_frac="""
        1

        $atsymb 0.0 0.0 0.0
        """, LatVecs=gen_lattice_fcc( FCC_DATA[atsymb]*ANG2BOHR ) )
        pspfiles = [ "../pseudopotentials/pade_gth/"*ALL_PADE_PSP[atsymb][1] ]
        ecutwfc = 15.0

        Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], extra_states=4 )
        println(Ham)

        @time KS_solve_SCF!( Ham, use_smearing=true, mix_method="rpulay" )
    end
end

function test_hcp()
    
    for atsymb in keys(HCP_DATA)

        a = HCP_DATA[atsymb][1]*ANG2BOHR
        c = HCP_DATA[atsymb][2]*ANG2BOHR

        atoms = Atoms( xyz_string_frac="""
        1

        $atsymb 0.0 0.0 0.0
        """, LatVecs=gen_lattice_hexagonal(a,c) )
        pspfiles = [ "../pseudopotentials/pade_gth/"*ALL_PADE_PSP[atsymb][1] ]
        ecutwfc = 15.0

        Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], extra_states=4 )
        println(Ham)

        @time KS_solve_SCF!( Ham, use_smearing=true, mix_method="rpulay" )
    end
end

function test_hcp_spinpol()
    
    for atsymb in keys(HCP_DATA)

        a = HCP_DATA[atsymb][1]*ANG2BOHR
        c = HCP_DATA[atsymb][2]*ANG2BOHR

        atoms = Atoms( xyz_string_frac="""
        1

        $atsymb 0.0 0.0 0.0
        """, LatVecs=gen_lattice_hexagonal(a,c) )
        pspfiles = [ "../pseudopotentials/pade_gth/"*ALL_PADE_PSP[atsymb][1] ]
        ecutwfc = 15.0

        Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], extra_states=4, Nspin=2 )
        println(Ham)

        @time KS_solve_SCF!( Ham, use_smearing=true, mix_method="rpulay" )
    end
end

#test_hcp()
#test_hcp_spinpol()
test_fcc()
#test_bcc()


