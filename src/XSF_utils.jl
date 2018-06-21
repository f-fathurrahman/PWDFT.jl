
# XSF spec: http://www.xcrysden.org/doc/XSF.html


function write_xsf( filenam::String, atoms::Atoms )
    write_xsf( filenam, atoms.LatVecs/ANG2BOHR, atoms.positions/ANG2BOHR;
               atsymbs=atoms.atsymbs )
end


function write_xsf( filnam::String, LL::Array{Float64,2}, atpos::Array{Float64,2};
                    atsymbs=nothing, molecule=false )
    #
    f = open(filnam, "w")
    Natoms = size(atpos)[2]
    #
    if molecule
        @printf(f, "MOLECULE\n")
    else
        @printf(f, "CRYSTAL\n")
    end
    v1 = LL[:,1]
    v2 = LL[:,2]
    v3 = LL[:,3]
    @printf(f, "PRIMVEC\n")
    @printf(f, "%18.10f %18.10f %18.10f\n", v1[1], v1[2], v1[3])
    @printf(f, "%18.10f %18.10f %18.10f\n", v2[1], v2[2], v2[3])
    @printf(f, "%18.10f %18.10f %18.10f\n", v3[1], v3[2], v3[3])
    @printf(f, "PRIMCOORD\n")
    @printf(f, "%8d %8d\n", Natoms, 1)
    #
    if atsymbs == nothing
        for ia = 1:Natoms
            @printf(f, "X  %18.10f %18.10f %18.10f\n", atpos[1,ia], atpos[2,ia], atpos[3,ia])
        end
    else
        for ia = 1:Natoms
            @printf(f, "%s  %18.10f %18.10f %18.10f\n", atsymbs[ia], atpos[1,ia], atpos[2,ia], atpos[3,ia])
        end
    end

    close(f)
end


function write_xsf_data3d_crystal(
        filnam::String, Ns::Tuple{Int64,Int64,Int64}, LL::Array{Float64,2},
        data3d::Array{Float64,1};
        center=zeros(3) )
    #
    f = open(filnam, "a")  # FIXME: What if filnam is not exist?
    #
    @printf(f, "BEGIN_BLOCK_DATAGRID_3D\n")
    @printf(f, "made_by_ffr\n")
    @printf(f, "DATAGRID_3D_UNKNOWN\n")
    @printf(f, "%8d %8d %8d\n", Ns[1]+1, Ns[2]+1, Ns[3]+1 )
    @printf(f, "%18.10f %18.10f %18.10f\n", center[1], center[2], center[3])
    v1 = LL[:,1]
    v2 = LL[:,2]
    v3 = LL[:,3]
    @printf(f, "%18.10f %18.10f %18.10f\n", v1[1], v1[2], v1[3])
    @printf(f, "%18.10f %18.10f %18.10f\n", v2[1], v2[2], v2[3])
    @printf(f, "%18.10f %18.10f %18.10f\n", v3[1], v3[2], v3[3])
    #
    rDat3d = reshape( data3d, (Ns[1],Ns[2],Ns[3]) )
    for k = 1:Ns[3]+1
        for j = 1:Ns[2]+1
            for i = 1:Ns[1]+1
                ii = ( i == Ns[1]+1 ? 1 : i )
                jj = ( j == Ns[2]+1 ? 1 : j )
                kk = ( k == Ns[3]+1 ? 1 : k )
                #
                @printf(f, "%18.10f\n", rDat3d[ii,jj,kk])
            end
        end
    end
    @printf(f, "END_DATAGRID_3D\n")
    @printf(f, "END_BLOCK_DATAGRID_3D\n")
    close(f)
end
