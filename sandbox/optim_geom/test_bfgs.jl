using Printf
using Random
using LinearAlgebra

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include(joinpath(DIR_PWDFT, "utilities", "PWSCF.jl"))

include("../NLopt_v3/calc_energies_grad.jl")
include("../NLopt_v3/KS_solve_Emin_PCG_new.jl")
include("../NLopt_v3/KS_solve_Emin_PCG_vec.jl")
include("update_positions.jl")

Random.seed!(1234)

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

function init_Ham_H2()
    # Atoms
    atoms = Atoms( xyz_file="H2.xyz", LatVecs=gen_lattice_sc(16.0) )
    # Initialize Hamiltonian
    ecutwfc = 15.0
    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )

    return Ham
end

function initial_hessian( Natoms )
    h = 70/(2*Ry2eV) * (ANG2BOHR^2)
    #h = 70.0
    H = diagm( 0 => h*ones(3*Natoms) )
    return H
end
    
function update_hessian( H_old, r, f, r0, f0 )
    Natoms = size(r,2)
    dr = r - r0
    # FIXME: need this?
    if maximum( abs.(dr) ) < 1e-7
        return diagm( 0 => 70*ones(3*Natoms) )
    end
    df = f - f0
    a = dot(dr, df)
    dg = H_old*dr
    b = dot(dr, dg)
    return H_old - (df*df')/a - (dg*dg')/b
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

function run_pwdft_jl( Ham )
    psiks = rand_BlochWavefunc(Ham)
    #KS_solve_Emin_PCG!( Ham, psiks )
    KS_solve_Emin_PCG_vec!( Ham, psiks )
    forces = calc_forces( Ham, psiks )
    return sum(Ham.energies), forces
end

function run_pwdft_jl!( Ham, psiks )
    KS_solve_Emin_PCG_vec!( Ham, psiks, skip_initial_diag=true )
    forces = calc_forces( Ham, psiks )
    return sum(Ham.energies), forces
end

function main()
    
    Ham = init_Ham_H2O()
    #Ham = init_Ham_H2()

    Natoms = Ham.atoms.Natoms
    
    #energies, forces = run_pwscf(Ham)
    #energies, forces = run_pwdft_jl(Ham)

    psiks = rand_BlochWavefunc(Ham)
    energies, forces = run_pwdft_jl!(Ham, psiks)

    println("Initial r  =")
    display(Ham.atoms.positions'); println()
    println("Initial forces = ")
    display(forces'); println()
    
    MAXSTEP = 0.04*ANG2BOHR

    f = vec(copy(forces))
    r = vec(copy(Ham.atoms.positions))
    r0 = zeros(size(r))
    f0 = zeros(size(f))

    H = initial_hessian(Natoms)

    NiterMax = 15
    for iter = 1:NiterMax

        println("Hessian = ")
        display(H); println()

        omega, V = eigen( Symmetric(H) )
        dr = V * (V'*f ./ abs.(omega))
        steplengths = sqrt.(sum( dr.^2, dims=1 ))
        maxsteplength = maximum(steplengths)
        if maxsteplength >= MAXSTEP
            println("Scaling dr")
            dr = dr * MAXSTEP / maxsteplength
        end

        r0 = copy(r)
        f0 = copy(f)
        H_old = copy(H)

        r[:] = r[:] + dr[:]
        #Ham.atoms.positions = copy(reshape(r,(3,Natoms)))
        update_positions!( Ham, reshape(r,3,Natoms) )

        energies_old = energies
        #energies, forces = run_pwscf(Ham)
        #energies, forces = run_pwdft_jl(Ham)
        energies, forces = run_pwdft_jl!(Ham, psiks)
        
        @printf("\nIter = %3d, Etot = %18.10f\n", iter, sum(energies))
        println("Forces = ")
        display(forces'); println()
        println("dr = ")
        display(reshape(dr,(3,Natoms))'); println()
        println("r  =")
        display(Ham.atoms.positions'); println()
        println("r (in Angstrom) =")
        display(Ham.atoms.positions' / ANG2BOHR); println()

        f = vec(copy(forces))
        H = update_hessian( H_old, r, f, r0, f0 )

    end

end

main()