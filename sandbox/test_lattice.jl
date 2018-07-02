using PWDFT

function test_fcc()
    LatVecs = gen_lattice_fcc(16.0)
    ecutwfc_Ry = 1.0

    pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs)
    println(pw)

    write_xsf("TEMP_fcc.xsf", LatVecs, pw.r)
end

function test_bcc()
    LatVecs = gen_lattice_bcc(16.0)
    ecutwfc_Ry = 1.0

    pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs)
    println(pw)

    write_xsf("TEMP_bcc.xsf", LatVecs, pw.r)
end

function test_bcc_v2()
    LatVecs = gen_lattice_bcc_v2(16.0)
    ecutwfc_Ry = 1.0

    pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs)
    println(pw)

    write_xsf("TEMP_bcc_v2.xsf", LatVecs, pw.r)
end

function test_hexagonal()
    LatVecs = gen_lattice_hexagonal(16.0,5.0)
    ecutwfc_Ry = 1.0

    pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs)
    println(pw)

    write_xsf("TEMP_hexagonal.xsf", LatVecs, pw.r)
end

function test_trigonal()
    #
    LatVecs = gen_lattice_trigonal(16.0, 30.0)
    ecutwfc_Ry = 1.0
    pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs)
    println(pw)
    write_xsf("TEMP_trigonal_30.xsf", LatVecs, pw.r)    
    #
    LatVecs = gen_lattice_trigonal(16.0, 45.0)
    ecutwfc_Ry = 1.0
    pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs)
    println(pw)
    write_xsf("TEMP_trigonal_45.xsf", LatVecs, pw.r)
    #
    LatVecs = gen_lattice_trigonal(16.0, 60.0)
    ecutwfc_Ry = 1.0
    pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs)
    println(pw)
    write_xsf("TEMP_trigonal_60.xsf", LatVecs, pw.r)
    #
    LatVecs = gen_lattice_trigonal(16.0, 90.0)
    ecutwfc_Ry = 1.0
    pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs)
    println(pw)
    write_xsf("TEMP_trigonal_90.xsf", LatVecs, pw.r)
end

function test_tetragonal_P()
    LatVecs = gen_lattice_tetragonal_P(16.0, 5.0)
    ecutwfc_Ry = 1.0

    pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs )
    println(pw)

    write_xsf("TEMP_tetragonal_P.xsf", LatVecs, pw.r)
end

function test_tetragonal_I()
    LatVecs = gen_lattice_tetragonal_I(16.0, 5.0)
    ecutwfc_Ry = 1.0

    pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs )
    println(pw)

    write_xsf("TEMP_tetragonal_I.xsf", LatVecs, pw.r)
end

function test_orthorhombic()
    LatVecs = gen_lattice_orthorhombic(16.0, 5.0, 7.0)
    ecutwfc_Ry = 1.0

    pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs )
    println(pw)

    write_xsf("TEMP_orthorhombic.xsf", LatVecs, pw.r)
end

function test_monoclinic()
    LatVecs = gen_lattice_monoclinic(11.0, 10.0, 5.0, 60.0)
    ecutwfc_Ry = 1.0

    pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs )
    println(pw)

    write_xsf("TEMP_monoclinic.xsf", LatVecs, pw.r)
end

function test_triclinic()
    LatVecs = gen_lattice_triclinic(11.0, 10.0, 5.0, 60.0, 90.0, 45.0)
    ecutwfc_Ry = 1.0

    pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs )
    println(pw)

    write_xsf("TEMP_triclinic.xsf", LatVecs, pw.r)
end

test_fcc()
test_bcc()
test_bcc_v2()
test_hexagonal()
test_trigonal()
test_tetragonal_P()
test_tetragonal_I()
test_orthorhombic()
test_monoclinic()
test_triclinic()
