using Printf
using PWDFT

include("../test/all_pade_psp.jl")

function test_atmlength()
    prefix = "../pseudopotentials/pade_gth/"
    atsymbs = keys(ALL_PADE_PSP)
    for atsymb in atsymbs
        pspfiles = ALL_PADE_PSP[atsymb]
        for fil in pspfiles
            pspot = PsPot_GTH(prefix*fil)
            zval = Float64(pspot.zval)
            znucl = Float64(PWDFT.ZATOMS[pspot.atsymb])
            Length = atmlength( 0.0, zval, znucl )
            @printf("%s %f %f %f\n", atsymb, zval, znucl, Length)
        end
    end 
end
#test_atmlength()


function test_guess_rhoe_H2()
    atoms = Atoms( xyz_file="../structures/H2.xyz", LatVecs=gen_lattice_sc(16.0) )
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, 15.0 )
    Rhoe = guess_rhoe( Ham )

    write_xsf( "TEMP_H2_guess_rhoe.xsf", atoms )
    write_xsf_data3d_crystal( "TEMP_H2_guess_rhoe.xsf", atoms, Ham.pw.Ns, Rhoe )
end
test_guess_rhoe_H2()

function test_guess_rhoe_PtO()
    atoms = Atoms( xyz_file="../structures/PtO.xyz", LatVecs=gen_lattice_sc(16.0) )
    pspfiles = ["../pseudopotentials/pade_gth/Pt-q18.gth",
                "../pseudopotentials/pade_gth/O-q6.gth"]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, 15.0 )
    Rhoe = guess_rhoe( Ham )

    write_xsf( "TEMP_PtO_guess_rhoe.xsf", atoms )
    write_xsf_data3d_crystal( "TEMP_PtO_guess_rhoe.xsf", atoms, Ham.pw.Ns, Rhoe )
end
test_guess_rhoe_PtO()
