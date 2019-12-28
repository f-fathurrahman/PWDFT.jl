using CUDAdrv

using PWDFT
using PWDFT_cuda

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

function main()

    atoms = Atoms(xyz_string=
    """
    3

    H   0.0  0.0  0.0
    H   0.0  1.51  0.0
    H   1.5  0.0  0.0
    """, LatVecs=gen_lattice_sc(10.0))

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 50.0

    Ham = CuHamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )

    Npoints = prod(Ham.pw.Ns)

    psiks = rand_CuBlochWavefunc( Ham )
    Rhoe = calc_rhoe( Ham, psiks )

    update!( Ham, Rhoe )
    
    Vpsi = op_V_loc( Ham, psiks[1] )

    println("Pass here")
end

CUDAdrv.@profile main()