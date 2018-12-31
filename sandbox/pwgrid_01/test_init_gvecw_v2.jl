#
# Test using k-points dependent Hamiltonian
#

using Printf
using PWDFT

function read_kpts( filname )
    file = open(filname)
    str = readline(file)
    Nkpts = parse( Int, str )
    kpts = zeros( Float64, 3,Nkpts )
    for ik = 1:Nkpts
        str = split(readline(file))
        kpts[1,ik] = parse( Float64, str[1] )
        kpts[2,ik] = parse( Float64, str[2] )
        kpts[3,ik] = parse( Float64, str[3] )
    end
    close(file)
    return kpts
end

function init_gvecw_kpts( ecutwfc, gvec::GVectors, kpts::Array{Float64,2} )
    G = gvec.G
    Ng = gvec.Ng
    idx_g2r = gvec.idx_g2r
    Nkpts = size(kpts)[2]
    #
    Gk2 = zeros(Float64,Ng)
    Gk = zeros(Float64,3)
    idx_gw2g = Array{Array{Int64,1},1}(undef,Nkpts)
    idx_gw2r = Array{Array{Int64,1},1}(undef,Nkpts)
    Ngk = Array{Int64,1}(undef,Nkpts)
    #
    for ik = 1:Nkpts
        for ig = 1:Ng
            Gk[:] = G[:,ig] .+ kpts[:,ik]
            Gk2[ig] = Gk[1]^2 + Gk[2]^2 + Gk[3]^2
        end
        idx_gw2g[ik] = findall( 0.5*Gk2 .<= ecutwfc )
        idx_gw2r[ik] = idx_g2r[idx_gw2g[ik]]
        Ngk[ik] = length(idx_gw2g[ik])
        @printf("ik = %8d, k = [%10.5f,%10.5f,%10.5f]\n",
                ik, kpts[1,ik], kpts[2,ik], kpts[3,ik])
        @printf("Ngk = %8d\n\n", Ngk[ik])
    end
    Nelems = sum(Ngk)
    memMiB = 2*Nelems*sizeof(Int64)/1024.0/1024.0
    @printf("mem = %f MiB\n\n", memMiB)
end

function test_kpath()
    LatVecs = gen_lattice_fcc(16.0)
    ecutwfc = 30.0
    pw = PWGrid( ecutwfc, LatVecs )
    println(pw)

    kpts = read_kpts("KPATH_FCC_60")
    #kpts =  (kpts_red' * pw.RecVecs)'
    kpts = pw.RecVecs*kpts
    init_gvecw_kpts( ecutwfc, pw.gvec, kpts)
end


function test_kgrid()

    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(5.431*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions

    ecutwfc = 30.0
    pw = PWGrid( ecutwfc, atoms.LatVecs )
    println(pw)

    kmesh  = [3,3,3]
    kshift = [0,0,0]
    kpts, wk = gen_kgrid_reduced( atoms, kmesh, kshift )
    kpts = pw.RecVecs*kpts
    init_gvecw_kpts( ecutwfc, pw.gvec, kpts )
end

test_kpath()
test_kgrid()


