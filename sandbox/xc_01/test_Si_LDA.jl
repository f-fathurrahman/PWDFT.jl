using LinearAlgebra
using Random
using Printf
using Libxc
using PWDFT

import PWDFT: guess_rhoe_atomic

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("Libxc_LDA.jl")

include("XC_x_slater.jl")
include("XC_c_pw.jl")
include("XC_c_vwn.jl")

include("XC_x_slater_spin.jl")
include("XC_c_pw_spin.jl")
include("XC_c_vwn_spin.jl")

include("XC_x_pbe.jl")
include("XC_c_pbe.jl")

include("XC_c_pbe_spin.jl")

function test_Si()

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

    third = 1.0/3.0
    pi34 = 0.6203504908994  # pi34=(3/4pi)^(1/3)    
    etxc = 0.0
    for ip in 1:Npoints
        
        rs = pi34/Rhoe[ip,1]^third
        
        ss_xc, vxc = XC_c_pw( rs )
        
        etxc = etxc + ss_xc*Rhoe[ip]
    end

    epsxc = calc_epsxc_LDA( Rhoe[:,1], Libxc.LDA_C_PW )
    E_xc_v1 = dot( Rhoe, epsxc ) * dVol
    @printf("Exc v1 = %18.10f\n", E_xc_v1)

    E_xc_v2 = etxc*dVol
    @printf("Exc v2 = %18.10f\n", E_xc_v2)
end

test_Si()
