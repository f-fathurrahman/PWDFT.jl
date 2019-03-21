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


"""
Write GTH/HGH pseudopotential in the format that can be read by ABINIT and PWSCF.
"""
function write_psp10( psp::PsPot_GTH; prefix="", extension="psp10" )
    atsymb = psp.atsymb
    zval = psp.zval
    zatom = PWDFT.ZATOMS[atsymb]
    lmax = psp.lmax
    
    if extension == ".gth"
        f = open( prefix*atsymb*".gth", "w" )
    else
        f = open( prefix*atsymb*".psp10", "w" )
    end
    @printf(f, "Hartwigsen-Goedecker-Hutter psp for %s, from PRB58, 3641 (1998)\n", atsymb)
    @printf(f, "%d %d 010605\n", Int(zatom), Int(zval))
    @printf(f, "10 1 %d 0 2001 0\n", lmax)
    
    # local pseudopotential
    rlocal = psp.rlocal
    c = psp.c
    ncoef = count(abs.(c) .> 0.0)
    @printf(f, "%18.10f %5d", rlocal, ncoef)
    for i = 1:ncoef
        @printf(f, "%18.10f ", c[i])
    end
    @printf(f, "\n")


    # Nonlocal part
    Nproj_l = psp.Nproj_l
    rc = psp.rc
    h = psp.h
    @printf(f, "%5d\n", lmax+1)
    for l = 0:lmax
        @printf(f, "%18.10f %5d", rc[l+1], Nproj_l[l+1])
        for i = 1:Nproj_l[l+1]
            for j = i:Nproj_l[l+1]
                @printf(f, "%18.10f ", h[l+1,i,j])
            end
            @printf(f, "               \n")
        end
        if l > 0
            for i = 1:Nproj_l[l+1]
                for j = i:Nproj_l[l+1]
                    @printf(f, "%18.10f ", 0.0)
                end
                @printf(f, "               \n")
            end
        end
    end

    close(f)

    return
end