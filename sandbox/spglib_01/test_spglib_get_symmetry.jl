using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

function spg_get_symmetry( atoms::Atoms; symprec=1e-5 )

    lattice = Matrix(atoms.LatVecs')
    positions = Matrix(inv(atoms.LatVecs))*atoms.positions # convert to fractional coordinates
    ctypes = Base.cconvert( Array{Int32,1}, atoms.atm2species)
    num_atom = Base.cconvert( Int32, atoms.Natoms )

    max_size = 50
    cmax_size = Base.cconvert(Int32, max_size)
    out_rot = zeros(Int32,3,3,max_size)
    out_translations = zeros(Float64,3,max_size)

    Nsym =
    ccall((:spg_get_symmetry, PWDFT.LIBSYMSPG), Int32,
          (Ptr{Int32}, Ptr{Float64}, Int32,
           Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Int32, Float64),
           out_rot, out_translations, cmax_size,
           lattice, positions, ctypes, num_atom, symprec)
    
    for i = 1:Nsym
        println("out_translations = ", out_translations[:,i])
    end

    for i = 1:Nsym
        println("rotation = ", i)
        println(out_rot[:,:,i])
    end

    println("Nsym = ", Nsym)

    return
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


    spg_get_symmetry(atoms)
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

    spg_get_symmetry(atoms)
end


function test_CH4()
    # Atoms
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "CH4.xyz"),
                   LatVecs=gen_lattice_cubic(16.0))
    println(atoms)
    spg_get_symmetry(atoms)
end


test_Si_fcc()
#test_GaAs()
#test_CH4()