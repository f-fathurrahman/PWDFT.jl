using PWDFT

using Random

Random.seed!(1234)

include("data_crystal.jl")
include("all_pade_psp.jl")

function test_conv_kpts( atsymb::String, lattice::String )
    
    if lattice == "hcp"
        a = HCP_DATA[atsymb][1]*ANG2BOHR
        c = HCP_DATA[atsymb][2]*ANG2BOHR
        LatVecs = gen_lattice_hexagonal(a,c)
    
    elseif lattice == "fcc"
        a = FCC_DATA[atsymb][1]*ANG2BOHR
        LatVecs = gen_lattice_fcc(a)
    
    elseif lattice == "bcc"
        a = BCC_DATA[atsymb][1]*ANG2BOHR
        LatVecs = gen_lattice_bcc(a)
    
    else
        error("Unknown lattice: ", lattice)
    end

    pspfiles = [ "../pseudopotentials/pade_gth/"*ALL_PADE_PSP[atsymb][1] ]

    atoms = Atoms( xyz_string_frac=
        """
        1

        $atsymb 0.0 0.0 0.0
        """, LatVecs=LatVecs )

    ecutwfc = 15.0
    Nk_start = 3
    Ndata = 10
    for i = 1:Ndata
        Nk = Nk_start + (i-1)
        Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[Nk,Nk,Nk], extra_states=4 )
        println(Ham)
        @time KS_solve_SCF!( Ham, use_smearing=true, kT=0.01, mix_method="rpulay", update_psi="LOBPCG" )
    end
end

test_conv_kpts("Mo", "bcc")
