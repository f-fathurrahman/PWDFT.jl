using Printf
using PWDFT

function write_psp_pwscf( psp::PsPot_GTH )
    atsymb = psp.atsymb
    zval = psp.zval
    zatom = PWDFT.ZATOMS[atsymb]
    lmax = psp.lmax
    
    f = open("TEMP_"*atsymb*".gth", "w")
    @printf(f, "Hartwigsen-Goedecker-Hutter psp for %s, from PRB58, 3641 (1998)\n", atsymb)
    @printf(f, "%d %d 010605\n", Int(zatom), Int(zval))
    @printf(f, "10 1 %d 0 2001 0\n", lmax)
    
    # local pseudopotential
    rlocal = psp.rlocal
    c = psp.c
    ncoef = count(abs.(c) .> 0.0)
    @printf(f, "%18.10f %d", rlocal, ncoef)
    for i = 1:ncoef
        @printf(f, "%18.10f ", c[i])
    end
    @printf(f, "\n")
    close(f)

    @printf("%d\n", lmax+1)
    for l = 0:lmax
    end
    return
end

function test_main()
    psp = PsPot_GTH("../pseudopotentials/pade_gth/Pd-q10.gth")
    println(psp)

    write_psp_pwscf(psp)
end

test_main()
