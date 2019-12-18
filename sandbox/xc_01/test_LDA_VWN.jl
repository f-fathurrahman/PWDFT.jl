using Random
using LinearAlgebra
using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("XCCalculator.jl")

include("LDA_VWN_pwdft_xc.jl")

function main_nospin()

    @printf("---------------\n")
    @printf("Testing LDA VWN\n")
    @printf("---------------\n")

    Random.seed!(1234)
    
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))
    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth", "Si-q4.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )

    pw = Ham.pw
    Npoints = prod(pw.Ns)
    dVol = pw.CellVolume/Npoints

    Rhoe = guess_rhoe_atomic(Ham)
    @printf("Integrated rhoe = %18.10f\n", sum(Rhoe)*dVol)

    epsxc = calc_epsxc_VWN( Ham.xc_calc, Rhoe )
    E_xc = dot( Rhoe, epsxc ) * dVol
    Vxc = calc_Vxc_VWN( Ham.xc_calc, Rhoe )

    epsxc = calc_epsxc_VWN( XCCalculator(), Rhoe )
    E_xc_v2 = dot( Rhoe, epsxc ) * dVol
    Vxc_v2 = calc_Vxc_VWN( XCCalculator(), Rhoe )

    @printf("E_xc v1 = %18.10f\n", E_xc)
    @printf("E_xc v2 = %18.10f\n", E_xc_v2)
    
    @printf("sum Vxc v1 = %18.10f\n", sum(Vxc))
    @printf("sum Vxc v2 = %18.10f\n", sum(Vxc_v2))

end

function main_spin()

    @printf("-----------------------\n")
    @printf("Testing LDA VWN spinpol\n")
    @printf("-----------------------\n")

    Random.seed!(1234)
    
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))
    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth", "Si-q4.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], Nspin=2 )

    pw = Ham.pw
    Npoints = prod(pw.Ns)
    dVol = pw.CellVolume/Npoints

    Rhoe = guess_rhoe_atomic(Ham, starting_magnetization=[0.1])
    @printf("Integrated rhoe = %18.10f\n", sum(Rhoe)*dVol)

    Rhoe_total = Rhoe[:,1] + Rhoe[:,2]

    epsxc = calc_epsxc_VWN( Ham.xc_calc, Rhoe )
    E_xc = dot( Rhoe_total, epsxc ) * dVol
    Vxc = calc_Vxc_VWN( Ham.xc_calc, Rhoe )

    epsxc = calc_epsxc_VWN( XCCalculator(), Rhoe )
    E_xc_v2 = dot( Rhoe_total, epsxc ) * dVol
    Vxc_v2 = calc_Vxc_VWN( XCCalculator(), Rhoe )

    @printf("E_xc v1 = %18.10f\n", E_xc)
    @printf("E_xc v2 = %18.10f\n", E_xc_v2)

    for ispin = 1:2
        @printf("sum Vxc v1 = %18.10f\n", sum(Vxc[:,ispin]))
        @printf("sum Vxc v2 = %18.10f\n", sum(Vxc_v2[:,ispin]))        
    end

end

main_nospin()
main_spin()