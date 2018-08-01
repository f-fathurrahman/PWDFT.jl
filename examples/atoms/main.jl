using PWDFT

const ALL_PS_PADE_GTH = """
Ag-q11.gth  Bi-q5.gth   Cs-q1.gth   Ge-q4.gth   Kr-q8.gth   Nb-q5.gth   P-q5.gth    Sc-q11.gth  Te-q6.gth   Zn-q12.gth
Ag-q19.gth  B-q3.gth    Cs-q9.gth   He-q2.gth   La-q11.gth  Nd-q14.gth  Pr-q13.gth  Sc-q3.gth   Ti-q12.gth  Zn-q20.gth
Ag-q1.gth   Br-q7.gth   Cu-q11.gth  Hf-q12.gth  Li-q1.gth   Ne-q8.gth   Pt-q10.gth  Se-q6.gth   Ti-q4.gth   Zn-q2.gth
Al-q3.gth   Ca-q10.gth  Cu-q19.gth  Hg-q12.gth  Li-q3.gth   Ni-q10.gth  Pt-q18.gth  Si-q4.gth   Tl-q13.gth  Zr-q12.gth
Ar-q8.gth   Ca-q2.gth   Cu-q1.gth   Hg-q2.gth   Lu-q25.gth  Ni-q18.gth  Rb-q1.gth   Sm-q16.gth  Tl-q3.gth   Zr-q4.gth
As-q5.gth   Cd-q12.gth  Dy-q20.gth  Ho-q21.gth  Mg-q10.gth  N-q5.gth    Rb-q9.gth   Sn-q4.gth   Tm-q23.gth  
At-q7.gth   Cd-q2.gth   Er-q22.gth  H-q1.gth    Mg-q2.gth   O-q6.gth    Re-q15.gth  S-q6.gth    V-q13.gth   
Au-q11.gth  Ce-q12.gth  Eu-q17.gth  In-q13.gth  Mn-q15.gth  Os-q16.gth  Re-q7.gth   Sr-q10.gth  V-q5.gth    
Au-q19.gth  Cl-q7.gth   Fe-q16.gth  In-q3.gth   Mn-q7.gth   Os-q8.gth   Rh-q17.gth  Sr-q2.gth   W-q14.gth   
Au-q1.gth   Co-q17.gth  Fe-q8.gth   I-q7.gth    Mo-q14.gth  Pb-q4.gth   Rh-q9.gth   Ta-q13.gth  W-q6.gth    
Ba-q10.gth  Co-q9.gth   F-q7.gth    Ir-q17.gth  Mo-q6.gth   Pd-q10.gth  Rn-q8.gth   Ta-q5.gth   Xe-q8.gth   
Ba-q2.gth   C-q4.gth    Ga-q13.gth  Ir-q9.gth   Na-q1.gth   Pd-q18.gth  Ru-q16.gth  Tb-q19.gth  Yb-q24.gth  
Be-q2.gth   Cr-q14.gth  Ga-q3.gth   K-q1.gth    Na-q9.gth   Pm-q15.gth  Ru-q8.gth   Tc-q15.gth  Y-q11.gth   
Be-q4.gth   Cr-q6.gth   Gd-q18.gth  K-q9.gth    Nb-q13.gth  Po-q6.gth   Sb-q5.gth   Tc-q7.gth   Y-q3.gth
"""

function do_run( psp_filename; ecutwfc_Ry=30.0 )

    # psp_filename is assumed to be a GTH pspot name (without full path)
    atsymb = split(psp_filename,"-")[1]

    # Atoms
    atoms = init_atoms_xyz_string(
    """
    1
    
    $atsymb   0.0   0.0   0.0
    """)
 
    atoms.LatVecs = gen_lattice_cubic(16.0)
    println(atoms)

    #
    # Initialize Hamiltonian
    #
    pspfiles = ["../../pseudopotentials/pade_gth/"*psp_filename]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, verbose=true,
                         Nspin=1 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    #KS_solve_SCF!( Ham, mix_method="simple", betamix=0.1 )
    KS_solve_Emin_PCG!( Ham, I_CG_BETA=4 )

    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands
    
    println("\nBand energies:")
    for ist = 1:Nstates
        @printf("%8d  %18.10f = %18.10f eV\n", ist, ebands[ist], ebands[ist]*Ry2eV*2)
    end
    
    println("\nTotal energy components")
    println(Ham.energies)

    pspcore_ene = calc_PspCore_ene(atoms, Ham.pspots, Ham.pw.Ω)
    
    @printf("\nPspCore ene = %18.10e\n", pspcore_ene)
    @printf("\nTotEne + PspCore = %18.10e\n", pspcore_ene + Ham.energies.Total)

end

function main()
    ps_name = split(ALL_PS_PADE_GTH)
    for p in ps_name
        do_run(p)
    end
end

#@time do_run("H-q1.gth")
#@time do_run("He-q2.gth")
#@time do_run("Li-q1.gth")
#@time do_run("Li-q3.gth")
#@time do_run("C-q4.gth", ecutwfc_Ry=30.0)
@time do_run("Si-q4.gth")
#@time do_run("Pd-q10.gth")

#main()

