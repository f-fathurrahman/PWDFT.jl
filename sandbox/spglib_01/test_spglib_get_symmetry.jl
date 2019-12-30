using PWDFT

import LibSymspg

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

function spg_get_symmetry( atoms::Atoms; symprec=1e-5 )

    lattice = Matrix(atoms.LatVecs')
    positions = Matrix(inv(atoms.LatVecs))*atoms.positions # convert to fractional coordinates
    ctypes = Base.cconvert( Array{Int32,1}, atoms.atm2species)
    num_atom = Base.cconvert( Int32, atoms.Natoms )

    max_size = 48*num_atom
    cmax_size = Base.cconvert(Int32, max_size)
    out_rot = zeros(Int32,3,3,max_size)
    out_translations = zeros(Float64,3,max_size)

    Nsyms =
    ccall((:spg_get_symmetry, LibSymspg.libsymspg), Int32,
          (Ptr{Int32}, Ptr{Float64}, Int32,
           Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Int32, Float64),
           out_rot, out_translations, cmax_size,
           lattice, positions, ctypes, num_atom, symprec)

    return Nsyms, out_rot[:,:,1:Nsyms], out_translations[:,1:Nsyms]
end


function test_Si_fcc()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true,
        LatVecs=gen_lattice_fcc(5.431*ANG2BOHR))
    
    println(atoms)
    Nsyms, rots, trans = spg_get_symmetry(atoms)
    println("Nsyms = ", Nsyms)
end


function test_Si8()
    atoms = Atoms(xyz_string=
        """
        8

        Ga     0.00000000000000     0.10000000000000     0.00000000000000
        Si     1.35750000000000     1.35750000000000     1.35750000000000
        Si     0.00000000000000     2.71500000000000     2.71500000000000
        Si     1.35750000000000     4.07250000000000     4.07250000000000
        Si     2.71500000000000     0.10000000000000     2.71500000000000
        Si     4.07250000000000     1.35750000000000     4.07250000000000
        Si     2.71500000000000     2.71500000000000     0.10000000000000
        Si     4.07250000000000     4.07250000000000     1.35750000000000
        """,
        LatVecs=gen_lattice_sc(5.430*ANG2BOHR))

    write_xsf("TEMP_Si8_pwdft.xsf", atoms)

    println(atoms)
    Nsyms, rots, trans = spg_get_symmetry(atoms)
    println("Nsyms = ", Nsyms)
end


function test_GaAs()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Ga  0.0   0.0   0.0
        As  0.25  0.25  0.25
        """, in_bohr=true,
        LatVecs=gen_lattice_fcc(10.6839444516))
    
    println(atoms)
    Nsyms, rots, trans = spg_get_symmetry(atoms)
    println("Nsyms = ", Nsyms)

    Nsyms2, rots2, trans2 = PWDFT.spg_get_symmetry(atoms)
    println("Nsyms2 = ", Nsyms2)

end


function test_CH4()
    # Atoms
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "CH4.xyz"),
                   LatVecs=gen_lattice_cubic(16.0))
    println(atoms)
    
    Nsyms, rots, trans = spg_get_symmetry(atoms)
    Nsyms = size(trans,2)

end


test_Si8()
#test_Si_fcc()
#test_GaAs()
#test_CH4()