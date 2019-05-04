using Printf
using LinearAlgebra
using Random
using SpecialFunctions: erf

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("Rhoe_aux.jl")

function init_Ham_Si_fcc( a::Float64, meshk::Array{Int64,1} )
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(a))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    #return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk, Ns_=(32,32,32) )
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk )
end



function init_Ham_H_in_box( a::Float64 )
    atoms = Atoms(xyz_string_frac=
        """
        1

        H  0.5  0.5  0.5
    """, in_bohr=true, LatVecs=gen_lattice_sc(a))

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end



function main()

    LATCONST = 10.2631
    Ham = init_Ham_Si_fcc( LATCONST, [3,3,3] )

    #Ham = init_Ham_H_in_box(16.0)

    @assert Ham.electrons.Nspin == 1

    pw = Ham.pw
    dVol = pw.CellVolume/prod(pw.Ns)
    println("dVol = ", dVol)

    psiks = rand_BlochWavefunc( Ham )
    Rhoe = calc_rhoe( Ham, psiks )
    println("integ of Rhoe = ", sum(Rhoe)*dVol)

    Rhoe_tot = Rhoe[:,1]
    Rhoe_tot_G = R_to_G(pw, Rhoe_tot)
    println("Rhoe_tot_G(G=0) = ", Rhoe_tot_G[1]*dVol)

    Rhoe_aux, E_self = gen_Rhoe_aux_G( Ham.atoms, Ham.pw, Ham.pspots )
    println("integ of Rhoe_aux = ", sum(Rhoe_aux)*dVol)
    println("E self = ", E_self)

    exit()

    V_aux_G = Poisson_solve( pw, Rhoe_aux )
    V_aux = real( G_to_R( pw, V_aux_G ) )
    println("integ of V_aux = ", sum(V_aux)*dVol)

    Rhoe_aux_R = gen_Rhoe_aux_R( Ham )
    V_aux_R = gen_V_aux_R( Ham )

    V_aux = gen_V_aux_G( Ham.atoms, Ham.pw )
    println("Some V_aux:")
    for i = 1:5
        @printf("%3d %18.10f\n", i, V_aux[i])
    end

    Rhoe_from_V_aux = -op_nabla_dot( pw, op_nabla(pw, V_aux) )/(4*pi)

    println("Some Rhoe_from_V_aux:")
    for i = 1:5
        @printf("%3d %18.10f %18.10f\n", i, Rhoe_from_V_aux[i], Rhoe_aux[i])
    end

    println
    println("sum = ", sum(Rhoe_from_V_aux)*dVol)

end

main()

