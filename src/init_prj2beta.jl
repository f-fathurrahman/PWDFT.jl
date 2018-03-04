function init_prj2beta( Pspots::Array{PsPot_GTH,1}, atoms::Atoms )

    Natoms = atoms.Natoms

    prj2beta = Array{Int64}( 3, Natoms, 4, 7 )
    prj2beta[:,:,:,:] = -1   # set to invalid index

    # 4: indexed from 0:3
    # 0, 1, 2, 3  -> l indexed
    # 1, 2, 3, 4  -> l + 1

    # -3:3
    # -3, -2, -1, 0, 1, 2, 3  -> m
    #  1,  2,  3, 4, 5, 6, 7  -> 4 + m, lmax = 3 + 1

    # -2, -1, 0, 1, 2  -> m
    #  1,  2, 3, 4, 5  -> 3 + m, lmax = 2 + 1

    atm2species = atoms.atm2species
    atpos = atoms.positions

    NbetaNL = 0
    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = PsPot_GTH[isp]
        for l = 0:psp.lmax
            for iprj = 1:psp.Nproj_l[l+1]
                for m = -l:l
                    NbetaNL = NbetaNL + 1
                    prj2beta[iprj,ia,l+1,m+psp.lmax+1] = NbetaNL
                    @printf("NbetaNL, l, m: %3d %3d %3d\n", NbetaNL, l, m)
                end
            end
        end
    end

    return NbetaNL, prj2beta
end

