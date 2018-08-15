using Printf
using PWDFT

function test_main()

    atoms = init_atoms_xyz_string("""
    1

    H  0.0  0.0  0.0
    """)
    atoms.LatVecs = gen_lattice_cubic(16.0)
    
    ecutwfc_Ry = 30.0
    kpoints = KPoints( atoms, [2,2,2], [0,0,0] )
    pw = PWGrid( ecutwfc_Ry*0.5, atoms.LatVecs, kpoints=kpoints )
    Ngw = pw.gvecw.Ngw
    Nkpt = pw.gvecw.kpoints.Nkpt

    Nstates = 4
    Focc = 2.0*ones(Nstates,Nkpt)

    psik = Array{Array{ComplexF64,2},1}(undef,Nkpt)
    for ik = 1:Nkpt
        psi = rand( ComplexF64, Ngw[ik], Nstates )
        psik[ik] = ortho_gram_schmidt(psi)  # orthogonalize in G-space
    end

    rhoe = calc_rhoe( pw, Focc, psik )
    dVol = pw.CellVolume/prod(pw.Ns)
    @printf("Integrated rhoe = %18.10f\n", sum(rhoe)*dVol)

end

test_main()

