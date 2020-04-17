using Printf
using Random
using LinearAlgebra

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include(joinpath(DIR_PWDFT, "utilities", "PWSCF.jl"))

function init_Ham_H2O()
    # Atoms
    atoms = Atoms( ext_xyz_file="H2O.xyz" )

    # Initialize Hamiltonian
    ecutwfc = 15.0
    pspfiles = [joinpath(DIR_PSP, "O-q6.gth"),
                joinpath(DIR_PSP, "H-q1.gth")]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )

    return Ham
end

function initial_hessian( Natoms )
    H = diagm( 0 => 70*ones(3*Natoms) )
    return H
end
    
function update_hessian( H_old, r, f, r0, f0 )
    Natoms = size(r,2)
    H_new = H_
    dr = r - r0
    if maximum( abs(dr) ) < 1e-7
        return diagm( 0 => 70*ones(3*Natoms) )
    end
    df = f - f0
    a = dot(dr, df)
    dg = H_old*dr
    b = dot(dr, dg)
    return H_old - (df * df')/a + (dg*dg')/ b
end

function run_pwscf( Ham )
    run(`rm -rfv TEMP_pwscf/\*`)
    write_pwscf( Ham, prefix_dir="TEMP_pwscf" )
    cd("./TEMP_pwscf")
    run(pipeline(`pw.x`, stdin="PWINPUT", stdout="LOG1"))
    cd("../")

    pwscf_energies = read_pwscf_etotal("TEMP_pwscf/LOG1")
    pwscf_forces = read_pwscf_forces("TEMP_pwscf/LOG1")
    return pwscf_energies, pwscf_forces
end

function main()
    Ham = init_Ham_H2O()
    Natoms = Ham.atoms.Natoms

    r = vec( copy(Ham.atoms.positions) )
    energies, forces = run_pwscf(Ham)
    display(forces); println()
    
    H = initial_hessian(Natoms)

    MAXSTEP = 0.04*ANG2BOHR

    f = vec(forces)
    omega, V = eigen(H)
    dr = V * (V*f ./ abs.(omega))
    steplengths = sqrt.(sum( dr.^2, dims=1 ))
    maxsteplength = maximum(steplengths)
    if maxsteplength >= MAXSTEP:
        dr = dr * MAXSTEP / maxsteplength

    println("dr = ")
    display(reshape(dr,(3,Natoms))); println()

    println("Before update")
    display(Ham.atoms.positions); println()
    #display(r); println()

    r = r + dr
    println("After update")
    Ham.atoms.positions = copy(reshape(r,(3,Natoms)))
    display(Ham.atoms.positions); println()
    #display(r); println()

    println("Pass here")
end

main()