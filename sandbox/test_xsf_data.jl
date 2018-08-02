using Printf
using Random
using LinearAlgebra

using PWDFT

include("alt1_KS_solve_SCF.jl")
include("../src/mix_rpulay.jl")

function create_Hamiltonian_CH4()
    # Atoms
    atoms = init_atoms_xyz("../structures/CH4.xyz")
    atoms.LatVecs = gen_lattice_cubic(16.0)
    println(atoms)

    # Initialize Hamiltonian
    ecutwfc_Ry = 30.0
    pspfiles = ["../pseudopotentials/pade_gth/C-q4.gth",
                "../pseudopotentials/pade_gth/H-q1.gth"]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, extra_states=2 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    return Ham
end

function do_calc()

    Ham = create_Hamiltonian_CH4()

    psiks = read_wfc(Ham)

    # Solve the KS problem
    @time KS_solve_SCF!(
        Ham, ETOT_CONV_THR=1e-6, NiterMax=50, betamix=0.5, update_psi="davidson",
        savewfc=true, startingwfc=psiks
    )
    
    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands
    
    println("\nBand energies:")
    for ist = 1:Nstates
        @printf("%8d  %18.10f = %18.10f eV\n", ist, ebands[ist], ebands[ist]*Ry2eV*2)
    end
    
    println("\nTotal energy components")
    println(Ham.energies)

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