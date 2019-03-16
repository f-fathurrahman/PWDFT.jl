# Display information about GTH pseudopotential
import Base: println
function println( psp::PsPot_GTH; header=true )

    ANGMOM = ["s", "p", "d", "f"]

    Nproj_l = psp.Nproj_l
    c = psp.c
    rlocal = psp.rlocal
    rc = psp.rc
    h = psp.h

    if header
        @printf("\n")
        @printf("                                   ---------\n")
        @printf("                                   PsPot_GTH\n")
        @printf("                                   ---------\n")
        @printf("\n")
    end
    @printf("Species: %s\n\n", psp.atsymb)
    @printf("zval: %d\n", psp.zval)
    @printf("File: %s\n", psp.pspfile)
    @printf("\nLocal pseudopotential parameters\n\n")
    @printf("rloc = %18.10f\n", rlocal)
    @printf("c[1] = %18.10f\n", c[1])
    @printf("c[2] = %18.10f\n", c[2])
    @printf("c[3] = %18.10f\n", c[3])
    @printf("c[4] = %18.10f\n", c[4])
    @printf("\n")
    if psp.lmax > -1
        @printf("Nonlocal pseudopotential parameters:\n\n")
    else
        @printf("No non-local pseudopotential.\n\n")
    end
    for i=1:psp.lmax+1
        @printf("Angular momentum: %s, rc = %f\n", ANGMOM[i], rc[i])
        if Nproj_l[i] >= 1
            @printf("h = \n")
        end
        for pi = 1:Nproj_l[i]
            for pj = 1:Nproj_l[i]
                @printf("%18.10f ", h[i,pi,pj])
            end
            @printf("\n")
        end
        @printf("\n")
    end

end
