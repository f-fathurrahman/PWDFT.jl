using Printf
using Random
using LinearAlgebra

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

function init_Ham_Si_fcc( a::Float64, meshk::Array{Int64,1} )
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(a))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk, Ns_=(32,32,32) )
end


function do_calc()

    LATCONST = 10.2631
    Ham = init_Ham_Si_fcc( LATCONST, [3,3,3] )
    println(Ham)

    # Solve the KS problem
    #KS_solve_SCF!(
    #    Ham, etot_conv_thr=1e-6, NiterMax=50, betamix=0.5, mix_method="rpulay",
    #    savewfc=true
    #)

    psiks = read_wfc(Ham)

    Nstates = Ham.electrons.Nstates
    Npoints = prod(Ham.pw.Ns)

    cpsi = zeros(ComplexF64, Npoints)

    ik = 4
    ist = 4
    psi = psiks[ik]
    pw = Ham.pw

    cpsi[:] .= 0.0 + im*0.0
    # Transform to real space
    idx = pw.gvecw.idx_gw2r[ik]
    cpsi[idx] = psi[:,ist]

    G_to_R!(pw, cpsi)
    ss = norm(cpsi)
    cpsi = sqrt(Npoints/pw.CellVolume).*cpsi/ss
    psiR_real = real(cpsi)
    psiR_imag = imag(cpsi)

    println("integ cpsi = ", sum(abs.(cpsi).^2)*pw.CellVolume/Npoints)
    println("integ real = ", sum(psiR_real.^2)*pw.CellVolume/Npoints)
    println("integ imag = ", sum(psiR_imag.^2)*pw.CellVolume/Npoints)
    #
    filname = "TEMP_psi_real_"*string(ik)*"_"*string(ist)*".xsf"
    write_xsf( filname, Ham.atoms )
    write_xsf_data3d_crystal( filname, Ham.atoms, Ham.pw.Ns, psiR_real )
    #
    filname = "TEMP_psi_imag_"*string(ik)*"_"*string(ist)*".xsf"
    write_xsf( filname, Ham.atoms )
    write_xsf_data3d_crystal( filname, Ham.atoms, Ham.pw.Ns, psiR_imag )

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