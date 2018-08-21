using PWDFT

function test_main()

    atoms = init_atoms_xyz("../structures/H.xyz")
    atoms.LatVecs = gen_lattice_sc(16.0)
    
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5 )

    KS_solve_Emin_PCG!( Ham, NiterMax=5, verbose=true, savewfc=true )

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

    KS_solve_DCM!( Ham, NiterMax=15, startingwfc=psiks )

end

test_main()