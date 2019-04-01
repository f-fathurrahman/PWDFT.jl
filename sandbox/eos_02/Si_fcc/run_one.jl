using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include(joinpath(DIR_PWDFT, "sandbox", "ABINIT.jl"))
include(joinpath(DIR_PWDFT, "sandbox", "PWSCF.jl"))
include(joinpath(DIR_PWDFT, "sandbox", "KS_solve_SCF_potmix.jl"))

function init_Ham_Si_fcc( a::Float64, meshk::Array{Int64,1} )
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(a))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk )
end

function run_one( idata::Int64 )

    LATCONST = 10.2631

    scal_i = 0.8
    scal_f = 1.2
    Ndata = 15
    Δ = (scal_f - scal_i)/(Ndata-1)
    
    scal = scal_i + (idata-1)*Δ
    a = LATCONST*scal    
    Ham = init_Ham_Si_fcc( a, [3,3,3] )

    println(Ham)

    KS_solve_SCF_potmix!( Ham, betamix=0.2, print_final_ebands=false, startingrhoe=:random )

end

run_one(15)

