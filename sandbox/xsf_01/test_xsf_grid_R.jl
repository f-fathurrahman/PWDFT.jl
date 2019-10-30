using PWDFT

function main()

    LatVecs = gen_lattice_fcc(5.0)
    pw = PWGrid(2.0, LatVecs)
    atpos = zeros(3,1)

    println(pw)

    filename = "TEMP_fcc_grid_R.xsf"
    write_xsf( filename, LatVecs, pw.r )

end

main()