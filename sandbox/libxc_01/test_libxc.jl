using LinearAlgebra
using Random
using Printf
using PWDFT

import PWDFT: guess_rhoe_atomic

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

function test_LDA_VWN()

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
    dVol = pw.CellVolume/prod(pw.Ns)

    Rhoe = guess_rhoe_atomic(Ham)
    @printf("Integrated rhoe = %18.10f\n", sum(Rhoe)*dVol)

    Vxc = calc_Vxc_VWN( Rhoe )
    @printf("sum Vxc = %18.10f\n", sum(Vxc))

    epsxc = calc_epsxc_VWN( Rhoe )
    E_xc = dot( Rhoe, epsxc ) * dVol
    @printf("E_xc = %18.10f\n", E_xc)

    psiks = rand_BlochWavefunc(Ham)
    update!( Ham, Rhoe )

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Ngw = Ham.pw.gvecw.Ngw
    Nstates = Ham.electrons.Nstates

    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin-1)*Nkpt
        psiks[ikspin] = ones(Ngw[ik], Nstates)/Ngw[ik]
        psiks[ikspin] = op_H(Ham, psiks[ikspin])
        ss = sum(psiks[ikspin])
        @printf("sum psiks[ikspin] = [%18.10f,%18.10f]\n", real(ss), imag(ss))
    end
    end

end


function test_LDA_VWN_spinpol()

    @printf("-----------------------\n")
    @printf("Testing LDA VWN spinpol\n")
    @printf("-----------------------\n")

    Random.seed!(1234)
    
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))
    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth", "Si-q4.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], Nspin=2 )

    pw = Ham.pw
    dVol = pw.CellVolume/prod(pw.Ns)

    Rhoe = guess_rhoe_atomic(Ham, starting_magnetization=[0.0])
    @printf("Integrated rhoe = %18.10f\n", sum(Rhoe)*dVol)

    Vxc = calc_Vxc_VWN( Rhoe )
    for ispin = 1:2
        @printf("sum Vxc = %18.10f\n", sum(Vxc[:,ispin]))
    end

    epsxc = calc_epsxc_VWN( Rhoe )
    Rhoe_total = Rhoe[:,1] + Rhoe[:,2]
    E_xc = dot( Rhoe_total, epsxc ) * dVol
    @printf("E_xc = %18.10f\n", E_xc)

    psiks = rand_BlochWavefunc(Ham)
    update!( Ham, Rhoe )

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Ngw = Ham.pw.gvecw.Ngw
    Nstates = Ham.electrons.Nstates

    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin-1)*Nkpt
        psiks[ikspin] = ones(Ngw[ik], Nstates)/Ngw[ik]
        #
        psiks[ikspin] = op_H(Ham, psiks[ikspin])
        ss = sum(psiks[ikspin])
        @printf("sum psiks[ikspin] = [%18.10f,%18.10f]\n", real(ss), imag(ss))
    end
    end

end



function test_GGA_PBE()

    @printf("---------------\n")
    @printf("Testing GGA PBE\n")
    @printf("---------------\n")

    Random.seed!(1234)
    
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))
    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pbe_gth", "Si-q4.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[1,1,1], xcfunc="PBE" )

    pw = Ham.pw
    dVol = pw.CellVolume/prod(pw.Ns)

    Rhoe = guess_rhoe_atomic(Ham)
    @printf("Integrated rhoe = %18.10f\n", sum(Rhoe)*dVol)

    Vxc = calc_Vxc_PBE( pw, Rhoe )

    epsxc = calc_epsxc_PBE( pw, Rhoe )
    E_xc = dot( Rhoe, epsxc ) * dVol
    @printf("E_xc = %18.10f\n", E_xc)

    psiks = rand_BlochWavefunc(Ham)
    update!( Ham, Rhoe )

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Ngw = Ham.pw.gvecw.Ngw
    Nstates = Ham.electrons.Nstates

    for ispin = 1:Nspin
        for ik = 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin-1)*Nkpt
            psiks[ikspin] = ones(Ngw[ik], Nstates)/Ngw[ik]
            #
            psiks[ikspin] = op_V_loc(Ham, psiks[ikspin])
            ss = sum(psiks[ikspin])
            @printf("sum psiks[ikspin] = [%18.10f,%18.10f]\n", real(ss), imag(ss))
        end
    end

end



function test_GGA_PBE_spinpol()
    
    @printf("-----------------------\n")
    @printf("Testing GGA PBE spinpol\n")
    @printf("-----------------------\n")

    Random.seed!(1234)
    
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pbe_gth", "Si-q4.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[1,1,1], xcfunc="PBE", Nspin=2 )

    pw = Ham.pw
    dVol = pw.CellVolume/prod(pw.Ns)

    Rhoe = guess_rhoe_atomic(Ham, starting_magnetization=[0.0])

    @printf("Integrated rhoe = %18.10f\n", sum(Rhoe)*dVol)

    Vxc = calc_Vxc_PBE( pw, Rhoe )

    Npoints = prod(Ham.pw.Ns)
    f = open("TEMP_Vxc_spinpol.dat", "w")
    for ip = 1:Npoints
        @printf(f, "%d %18.10f %18.10f\n", ip, Vxc[ip,1], Vxc[ip,2])
    end
    close(f)

    epsxc = calc_epsxc_PBE( pw, Rhoe )
    Rhoe_total = Rhoe[:,1] + Rhoe[:,2]
    E_xc = dot( Rhoe_total, epsxc ) * dVol
    @printf("E_xc = %18.10f\n", E_xc)

    psiks = rand_BlochWavefunc(Ham)
    update!( Ham, Rhoe )

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Ngw = Ham.pw.gvecw.Ngw
    Nstates = Ham.electrons.Nstates

    for ispin = 1:Nspin
        for ik = 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin-1)*Nkpt
            psiks[ikspin] = ones(Ngw[ik], Nstates)/Ngw[ik]
            psiks[ikspin] = op_V_loc(Ham, psiks[ikspin])
            ss = sum(psiks[ikspin])
            @printf("sum psiks[ikspin] = [%18.10f,%18.10f]\n", real(ss), imag(ss))
        end
    end

end


#test_LDA_VWN()
#test_LDA_VWN_spinpol()
test_GGA_PBE()
test_GGA_PBE_spinpol()

#test_spinpol(xc="VWN", Nspin=1)
#test_spinpol(xc="VWN", Nspin=2)
#test_spinpol(xc="PBE", Nspin=1)
#test_spinpol(xc="PBE", Nspin=2)
