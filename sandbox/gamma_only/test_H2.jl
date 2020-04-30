using LinearAlgebra
using Random
using Printf

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

function main()

    Random.seed!(1234)

    # Atoms
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
                   LatVecs = gen_lattice_sc(16.0) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )

    println(Ham)
    exit()

    #psiks = rand_BlochWavefunc(Ham)
    #KS_solve_SCF!( Ham, psiks, mix_method="pulay" )

    pw = Ham.pw
    idx_gw2r = Ham.pw.gvecw.idx_gw2r[1]
    Ngw = Ham.pw.gvecw.Ngw[1]
    Npoints = prod(pw.Ns)
    
    psi_R = rand(Float64,Npoints)
    
    ctmp = R_to_G(pw, psi_R)
    
    psi_G = zeros(ComplexF64,Ngw)
    for igw in 1:Ngw
        ip = idx_gw2r[igw]
        psi_G[igw] = ctmp[ip]
    end

    ctmp2 = zeros(ComplexF64,Npoints)
    for igw in 1:Ngw
        ip = idx_gw2r[igw]
        ctmp2[ip] = psi_G[igw]
    end

    psi_R_v2 = G_to_R(pw, ctmp2)

    psi_R_v3 = G_to_R(pw, ctmp)

    println(psi_R[2])
    println(psi_R_v2[2])
    println(psi_R_v3[2])


    println("Pass here")
end

main()
