using Printf
using PWDFT

function test_Si_diamond()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(5.431*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)

    kpoints = KPoints( atoms, [3,3,3], [0,0,0], verbose=true )
    println(kpoints)
end

function test_Fe_bcc()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        1

        Fe  0.0   0.0   0.0
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_bcc(5.0*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)

    kpoints = KPoints( atoms, [3,3,3], [0,0,0], verbose=true )
    println(kpoints)
end

function test_kpath()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(5.431*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)

    kdict = get_special_kpoints("fcc")
    println(kdict)
end


function kpoints_from_file(atoms::Atoms, filename::String)
    file = open(filename)
    str = readline(file)
    Nkpt = parse( Int, str )
    kred = zeros( Float64, 3,Nkpt )
    for ik = 1:Nkpt
        str = split(readline(file))
        kred[1,ik] = parse( Float64, str[1] )
        kred[2,ik] = parse( Float64, str[2] )
        kred[3,ik] = parse( Float64, str[3] )
    end
    close(file)
    # kpts need to be converted to Cartesian form
    RecVecs = 2*pi*PWDFT.invTrans_m3x3(atoms.LatVecs)
    kpts = RecVecs*kred
    #
    wk = ones(Nkpt) # not used for non-scf calculations
    #
    return KPoints(Nkpt, kred, wk, RecVecs)
end

function test_kpoints_from_file()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(5.431*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)
    #
    kpoints = kpoints_from_file(atoms, "KPATH_FCC_60")
    println(kpoints)
    #
    pw = PWGrid(15.0, atoms.LatVecs, kpoints=kpoints)
    println(pw)
end

test_kpoints_from_file()

#test_Si_diamond()

#test_Fe_bcc()

#test_kpath()
