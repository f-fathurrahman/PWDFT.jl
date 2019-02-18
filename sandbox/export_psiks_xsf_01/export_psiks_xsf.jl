using Printf
using Random
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("../forces_nonlocal_01/read_psiks.jl")

function export_psiks_xsf(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    ik_range,
    ist_range
)

    pw = Ham.pw
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    @assert Nspin == 1

    cpsi = zeros(ComplexF64, Npoints)

    for ik in ik_range
    for ist in ist_range

        ikspin = ik + (Nspin-1)*Nkpt
        psi = psiks[ikspin]

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
        fname = "TEMP_psi_"*string(ik)*"_"*string(ist)*"_real.xsf"
        write_xsf( fname, Ham.atoms )
        write_xsf_data3d_crystal( fname, Ham.atoms, Ham.pw.Ns, psiR_real )
        #
        fname = "TEMP_psi_"*string(ik)*"_"*string(ist)*"_imag.xsf"
        write_xsf( "TEMP_psi_imag.xsf", Ham.atoms )
        write_xsf_data3d_crystal( "TEMP_psi_imag.xsf", Ham.atoms, Ham.pw.Ns, psiR_imag )
    end
    end

end


function main()
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
    println(Ham)

    Random.seed!(1234)
    
    #KS_solve_SCF!( Ham, mix_method="rpulay", savewfc=true )

    psiks = read_psiks(Ham)

    export_psiks_xsf( Ham, psiks, 1:4, 1:4 )

    println("Program ended normally")

end

main()
