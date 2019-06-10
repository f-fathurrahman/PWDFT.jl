using Printf
using Random
using LinearAlgebra

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include(joinpath(DIR_PWDFT, "sandbox", "scf_01", "alt1_KS_solve_SCF.jl"))

function create_Hamiltonian_CH4()
    # Atoms
    atoms = init_atoms_xyz(joinpath(DIR_STRUCTURES,"CH4.xyz"))
    atoms.LatVecs = gen_lattice_cubic(16.0)

    # Initialize Hamiltonian
    ecutwfc_Ry = 30.0
    pspfiles = [joinpath(DIR_PSP,"C-q4.gth"),
                joinpath(DIR_PSP,"H-q1.gth")]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, extra_states=2 )

    return Ham
end

function do_calc()

    Ham = create_Hamiltonian_CH4()
    println(Ham)

    if isfile("WFC_ikspin_1.data")
        println("\nWave function is read from file\n")
        psiks = read_wfc(Ham)
    else
        psiks = nothing
    end

    # Solve the KS problem
    @time KS_solve_SCF!(
        Ham, etot_conv_thr=1e-6, NiterMax=50, betamix=0.5, update_psi="davidson",
        savewfc=true, startingwfc=psiks
    )

    psiks = read_wfc(Ham)

    psi = psiks[1]

    Nstates = Ham.electrons.Nstates
    Npoints = prod(Ham.pw.Ns)

    cpsi = zeros(ComplexF64, Npoints)

    ik = 1
    ist = 4
    pw = Ham.pw

    cpsi[:] .= 0.0 + im*0.0
    # Transform to real space
    idx = pw.gvecw.idx_gw2r[ik]
    cpsi[idx] = psi[:,ist]
    psiR_real = real(G_to_R(pw, cpsi))
    psiR_imag = imag(G_to_R(pw, cpsi))
    #
    psiR_real = sqrt(Npoints/pw.CellVolume)*psiR_real # normalize
    psiR_imag = sqrt(Npoints/pw.CellVolume)*psiR_imag # normalize
    #
    write_xsf( "TEMP_psi_real.xsf", Ham.atoms )
    write_xsf_data3d_crystal( "TEMP_psi_real.xsf", Ham.atoms, Ham.pw.Ns, psiR_real )
    #
    write_xsf( "TEMP_psi_imag.xsf", Ham.atoms )
    write_xsf_data3d_crystal( "TEMP_psi_imag.xsf", Ham.atoms, Ham.pw.Ns, psiR_imag )

end

function read_wfc(Ham)

    Nstates = Ham.electrons.Nstates
    Ngw = Ham.pw.gvecw.Ngw
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    psiks = Array{Array{ComplexF64,2},1}(undef,Nkpt)

    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        # Don't forget to use read mode
        wfc_file = open("WFC_ikspin_"*string(ikspin)*".data","r")
        psiks[ikspin] = Array{ComplexF64}(undef,Ngw[ik],Nstates)
        psiks[ikspin] = read!( wfc_file, psiks[ikspin] )
        close( wfc_file )
    end
    end

    return psiks
end


do_calc()