function write_KPoints( f::IOStream, kpoints::KPoints )
    write(f, kpoints.Nkpt)
    write(f, kpoints.mesh[1])
    write(f, kpoints.mesh[2])
    write(f, kpoints.mesh[3])
    write(f, kpoints.k)
    write(f, kpoints.wk)
    write(f, kpoints.RecVecs)
end

function read_KPoints( f::IOStream )
    
    tmpInt = Array{Int64}(undef,1)
    
    read!(f, tmpInt)
    Nkpt = tmpInt[1]
    
    k = Array{Float64}(undef,3,Nkpt)
    wk = Array{Float64}(undef,Nkpt)
    RecVecs = Array{Float64}(undef,3,3)

    read!(f, tmpInt)
    mesh1 = tmpInt[1]
    
    read!(f, tmpInt)
    mesh2 = tmpInt[1]
    
    read!(f, tmpInt)
    mesh3 = tmpInt[1]
    
    mesh = (mesh1, mesh2, mesh3)

    read!(f, k)
    read!(f, wk)
    read!(f, RecVecs)

    return KPoints(Nkpt, mesh, k, wk, RecVecs)
end

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
    ss = maximum(abs.(RecVecs))

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
