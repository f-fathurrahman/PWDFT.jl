import Base: println

"""
Display some information about an instance of `KPoints`.
"""
function println( kpoints::KPoints; header=true )

    if header
        @printf("\n")
        @printf("                                     -------\n")
        @printf("                                     KPoints\n")
        @printf("                                     -------\n")
        @printf("\n")
    end

    @printf("\n")
    if kpoints.mesh != (0,0,0)
        @printf("Mesh: (%4d,%4d,%4d)\n", kpoints.mesh[1], kpoints.mesh[2], kpoints.mesh[3])
    end
    @printf("Total number of kpoints = %d\n", kpoints.Nkpt )

    @printf("\n")
    @printf("kpoints in Cartesian coordinate (unscaled)\n")
    @printf("\n")
    kcart = copy(kpoints.k)
    for ik = 1:kpoints.Nkpt
        @printf("%4d [%14.10f %14.10f %14.10f] %14.10f\n",
                ik, kcart[1,ik], kcart[2,ik], kcart[3,ik], kpoints.wk[ik])
    end

    RecVecs = kpoints.RecVecs
    LatVecs = inv(Matrix(RecVecs'))
    alat = norm(LatVecs[:,1])
    ss = 1.0/alat

    # This is useful for comparison with pwscf
    @printf("\n")
    @printf("kpoints in Cartesian coordinate (scale: %f)\n", ss)
    @printf("\n")
    kcart = kcart/ss

    for ik = 1:kpoints.Nkpt
        @printf("%4d [%14.10f %14.10f %14.10f] %14.10f\n",
                ik, kcart[1,ik], kcart[2,ik], kcart[3,ik], kpoints.wk[ik])
    end
    
    @printf("\n")
    @printf("sum wk = %f\n", sum(kpoints.wk))
end
