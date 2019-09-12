import Base: show

"""
Display some information about an instance of `KPoints`.
"""
function show( io::IO, kpoints::KPoints; header=true )

    if header
        @printf(io, "\n")
        @printf(io, "                                     -------\n")
        @printf(io, "                                     KPoints\n")
        @printf(io, "                                     -------\n")
        @printf(io, "\n")
    end

    @printf("\n")
    if kpoints.mesh != (0,0,0)
        @printf(io, "Mesh: (%4d,%4d,%4d)\n", kpoints.mesh[1], kpoints.mesh[2], kpoints.mesh[3])
    end
    @printf(io, "Total number of kpoints = %d\n", kpoints.Nkpt )

    @printf(io, "\n")
    @printf(io, "kpoints in Cartesian coordinate (unscaled)\n")
    @printf(io, "\n")
    kcart = copy(kpoints.k)
    for ik = 1:kpoints.Nkpt
        @printf(io, "%4d [%14.10f %14.10f %14.10f] %14.10f\n",
                ik, kcart[1,ik], kcart[2,ik], kcart[3,ik], kpoints.wk[ik])
    end

    RecVecs = kpoints.RecVecs
    LatVecs = inv(Matrix(RecVecs'))
    alat = norm(LatVecs[:,1])
    ss = 1.0/alat

    # This is useful for comparison with pwscf
    @printf(io, "\n")
    @printf(io, "kpoints in Cartesian coordinate (scale: %f)\n", ss)
    @printf(io, "\n")
    kcart = kcart/ss

    for ik = 1:kpoints.Nkpt
        @printf(io, "%4d [%14.10f %14.10f %14.10f] %14.10f\n",
                ik, kcart[1,ik], kcart[2,ik], kcart[3,ik], kpoints.wk[ik])
    end
    
    @printf(io, "\n")
    @printf(io, "sum wk = %f\n", sum(kpoints.wk))
end
show( kpoints::KPoints; header=true ) = show( stdout, kpoints, header=header )