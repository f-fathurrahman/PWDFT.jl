using Test

using PWDFT
using PWDFT_cuda

function main()
    LatVecs = gen_lattice_sc(10.0)
    atoms = Atoms( xyz_string=
    """
    1

    H   0.0  0.0  0.0
    """, LatVecs=LatVecs )

    pw = CuPWGrid( 15.0, LatVecs )
    Sf = calc_strfact( atoms, pw )

    pw_cpu = PWGrid( 15.0, LatVecs )
    Sf_cpu = calc_strfact( atoms, pw_cpu )

    Sf_gpu = collect( Sf )
    @test Sf_gpu ≈ Sf_cpu

    println("Pass here")
end

main()
