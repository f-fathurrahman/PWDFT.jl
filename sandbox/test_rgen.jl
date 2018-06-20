using Printf
using LinearAlgebra
using PWDFT

#=
rgen is presently not used now. It is remnant of required functions for calc_Ewald_qe
which is adapted from QE.
It is still included here because it might be used again to test Ewald force calculation
againsts QE.
=#
include("../src/rgen.jl")

function test_rgen()
    at = gen_lattice_sc(16.0)

    bg = zeros(Float64,3,3)

    TPI = 2*pi
    bg[1,1] = TPI/at[1,1]
    bg[2,2] = TPI/at[2,2]
    bg[3,3] = TPI/at[3,3]

    alat = 1.0
    alpha = 0.015
    rmax = 4.0/sqrt(alpha)/alat
    @printf("rmax = %18.10f\n", rmax)
    dtau = [0.1, 0.1, 0.2]
    mxr = 50

    nrm, r, r2 = rgen( dtau, rmax, mxr, at, bg )
    @printf("nrm = %d\n", nrm)

    for ir = 1:nrm
        @printf("%4d%18.10f%18.10f%18.10f%18.10f\n", ir, r[1,ir], r[2,ir], r[3,ir], r2[ir])
    end

end

test_rgen()
