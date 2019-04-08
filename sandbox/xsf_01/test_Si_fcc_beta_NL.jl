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
        1

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

    Npoints = prod(Ham.pw.Ns)

    ik = 2
    ibetaNL = 2
    psi = Ham.pspotNL.betaNL[:,ibetaNL,ik]
    cpsi = zeros(ComplexF64, Npoints)

    pw = Ham.pw
    Ngw = pw.gvecw.Ngw

    cpsi[:] .= 0.0 + im*0.0
    # Transform to real space
    idx = pw.gvecw.idx_gw2r[ik]
    cpsi[idx] = psi[1:Ngw[ik]]

    G_to_R!(pw, cpsi)
    ss = norm(cpsi)
    cpsi = sqrt(Npoints/pw.CellVolume).*cpsi/ss
    psiR_real = real(cpsi)
    psiR_imag = imag(cpsi)

    println("integ cpsi = ", sum(abs.(cpsi).^2)*pw.CellVolume/Npoints)
    println("integ real = ", sum(psiR_real.^2)*pw.CellVolume/Npoints)
    println("integ imag = ", sum(psiR_imag.^2)*pw.CellVolume/Npoints)
    #
    filname = "TEMP_betaNL_real.xsf"
    write_xsf( filname, Ham.atoms )
    write_xsf_data3d_crystal( filname, Ham.atoms, Ham.pw.Ns, psiR_real )
    #
    filname = "TEMP_betaNL_imag.xsf"
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