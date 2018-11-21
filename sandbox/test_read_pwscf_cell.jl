import Glob
filelist = Glob.glob("*.in")

stdout_orig = stdout

for f in filelist
    atsymb = replace(f, ".in" => "")
    
    fileout = open("LOG_"*atsymb, "w")
    
    redirect_stdout(fileout)
    
    atoms, meshk = read_pwscf_cell(f)

    pspfiles = [ "../pseudopotentials/pade_gth/"*ALL_PADE_PSP[atsymb][1] ]
    ecutwfc = 15.0
    Ham = Hamiltonian(atoms, pspfiles, ecutwfc, meshk=meshk, extra_states=4)
    println(Ham)

    KS_solve_SCF!( Ham, mix_method="rpulay", betamix=0.1, use_smearing=true )
    
    close(fileout)
end
