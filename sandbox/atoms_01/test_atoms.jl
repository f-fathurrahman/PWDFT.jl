using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

function test01()
    filename = joinpath(DIR_STRUCTURES,"CuSO4.xyz")
    atoms = init_atoms_xyz(filename, verbose=true)
    println(atoms)
    
    dummy_atoms = Atoms()
    println(dummy_atoms)
end

function test02()
    filename = joinpath(DIR_STRUCTURES, "CuSO4.xyz")
    atoms = Atoms(xyz_file=filename)
    println(atoms)
end


function test03()
    atoms = Atoms(xyz_string="""
    2

    H  0.0  0.0  0.0
    Cl 1.5  0.0  0.0
    """, in_bohr=true)
    println(atoms)
end


function test04()
    atoms = Atoms(xyz_string_frac="""
    2

    Si  0.0   0.0   0.0
    Si  0.25  0.25  0.25
    """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))
    println(atoms)
    write_xsf("TEMP_Si_fcc.xsf", atoms)
end

function test05()
    atoms = Atoms(xyz_string_frac="""
    2

    Ga  0.0   0.0   0.0
    As  0.25  0.25  0.25
    """, in_bohr=true, LatVecs=gen_lattice_fcc(5.6536*ANG2BOHR))
    println(atoms)
    write_xsf("TEMP_GaAs_fcc.xsf", atoms)
end

test01()
test02()
test03()
test04()
test05()

