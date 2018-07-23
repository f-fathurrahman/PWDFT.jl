using Printf
using Random
using PWDFT

function analytic_free_band( pw::PWGrid, Nstates::Int64, ik=1 )
    Ngw_ik = pw.gvecw.Ngw[ik]
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    k = pw.gvecw.kpoints.k

    Gw = zeros(3,Ngw_ik)
    Gw2 = zeros(Ngw_ik)
    G = pw.gvec.G
    for igk = 1:Ngw_ik
        ig = idx_gw2g[igk]
        Gw[:,igk] = G[:,ig] + k[:,ik]
        Gw2[igk] = Gw[1,igk]^2 + Gw[2,igk]^2 + Gw[3,igk]^2
    end
    Gw2_sorted = sort(Gw2)
    evals = zeros(Nstates)
    for ist = 1:Nstates
        evals[ist] = 0.5*Gw2_sorted[ist]
    end
    return evals
end

# Test free_electron_Hamiltonian
function test_01()

    atoms = init_atoms_xyz_string(
        """
        1

        H  0.0   0.0   0.0
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(6.0)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)

    srand(1234)
    kpoints = KPoints(atoms)
    kpoints.k = kpoints.RecVecs*rand( 0.0:0.01:0.5, (3,1) )
    println(kpoints)

    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc = 50.0
    Ham = free_electron_Hamiltonian(
            atoms, pspfiles, ecutwfc, kpoints=kpoints, extra_states=4,
            verbose=true
          )

    Ngwx = Ham.pw.gvecw.Ngwx
    Nstates = Ham.electrons.Nstates
    psi = ortho_sqrt( rand(ComplexF64,Ngwx,Nstates) )

    #evals, psi = diag_lobpcg(Ham, psi, verbose=true, Nstates_conv=1)
    #evals, psi = diag_davidson(Ham, psi, verbose=true)
    evals, psi = diag_Emin_PCG(Ham, psi, verbose=true, NiterMax=400)

    evals_analytic = analytic_free_band( Ham.pw, Nstates )

    for ist = 1:Nstates
        @printf("%4d %18.10f %18.10f\n", ist, evals[ist], evals_analytic[ist])
    end

end


# Test with proper pseudopotentials
function test_02( diag_method::String )

    atoms = init_atoms_xyz_string(
        """
        1

        H  0.0   0.0   0.0
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(6.0)
    atoms.positions = atoms.LatVecs*atoms.positions

    srand(1234)
    kpoints = KPoints(atoms)
    kpoints.k = kpoints.RecVecs*rand( 0.0:0.01:0.5, (3,1) )

    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc = 50.0
    Ham = Hamiltonian(
            atoms, pspfiles, ecutwfc, kpoints=kpoints, extra_states=4,
          )

    Ngwx = Ham.pw.gvecw.Ngwx
    Nstates = Ham.electrons.Nstates
    psi = ortho_sqrt( rand(ComplexF64,Ngwx,Nstates) )

    if diag_method == "LOBPCG"
        evals, psi = diag_lobpcg(Ham, psi, verbose=true, Nstates_conv=1)
    elseif diag_method == "davidson"
        evals, psi = diag_davidson(Ham, psi, verbose=true)
    elseif diag_method == "PCG"
        evals, psi = diag_Emin_PCG(Ham, psi, verbose=true, NiterMax=400)
    else
        println("ERROR: unknown diag_method = ", diag_method)
        error("STOPPED")
    end

    for ist = 1:Nstates
        @printf("%4d %18.10f\n", ist, evals[ist])
    end

end

#test_01()
test_02("LOBPCG")
test_02("davidson")
test_02("PCG")