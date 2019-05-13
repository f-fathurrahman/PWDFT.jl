using Printf
using LinearAlgebra
using Random

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("Gshells.jl")

function init_Ham_Si_fcc()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true,
        LatVecs=gen_lattice_fcc(5.431*ANG2BOHR))

    ecutwfc = 15.0
    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
end


function init_Ham_GaAs()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Ga  0.0  0.0  0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.6839444516))

    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
end


function init_Ham_atom_H()
    atoms = Atoms( xyz_string="""
            1

            H   0.1   0.0   0.0
            """, LatVecs=gen_lattice_sc(16.0) )

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end

function init_Ham_H2()
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
                   LatVecs = gen_lattice_sc(16.0) )

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end

function init_Ham_H2_v2()
    atoms = Atoms( xyz_string=
        """
        2
            
        H  0.0  0.0  0.0
        H  1.5  0.0  0.0
        """, in_bohr=true, LatVecs = gen_lattice_sc(16.0) )

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end

function main()

    #Ham = init_Ham_Si_fcc()
    #Ham = init_Ham_GaAs()
    #Ham = init_Ham_atom_H()
    #Ham = init_Ham_H2()
    Ham = init_Ham_H2_v2()

    pw = Ham.pw
    atoms = Ham.atoms
    sym_info = Ham.sym_info
    println(sym_info)

    #G2_gshells, idx_gshells = init_Gshells(pw.gvec)
    #Ngshells = length(G2_gshells)
    #println("Ngshells = ", Ngshells)

    sym_gshell = Array{Int64,2}(undef,3,48)

    LatVecs = atoms.LatVecs
    G = pw.gvec.G

    G_crystal = zeros(Int64,3,pw.gvec.Ng)
    is_done = zeros(Bool,pw.gvec.Ng)
    for ig = 1:pw.gvec.Ng
        xx = ( LatVecs[1,:]*G[1,ig] + LatVecs[2,:]*G[2,ig] + LatVecs[3,:]*G[3,ig] )/(2*pi)
        G_crystal[:,ig] = round.(Int64,xx)
        is_done[ig] = false
        #println(G_crystal[:,ig])
    end


    Nsyms = sym_info.Nsyms
    s = sym_info.s
    inv_s = sym_info.inv_s

    println("Nsyms = ", Nsyms)

    Ngs = 0
    sn = zeros(Int64,3)

    shell_G_sym = Array{Array{Int64,1},1}(undef,pw.gvec.Ng)


    Ngtot = 0  # for verification
    for ig = 1:pw.gvec.Ng
        if is_done[ig]
            continue
        end
        Ngs = Ngs + 1
        ng = 0
        for isym = 1:Nsyms
            sn[:] = s[:,1,isym]*G_crystal[1,ig] + s[:,2,isym]*G_crystal[2,ig] + s[:,3,isym]*G_crystal[3,ig]
            found = false
            for i = 1:ng
                found = ( sn[:] == sym_gshell[:,i] )
                if found
                    break
                end
            end
            if !found
                #println("Adding to sym_gshell")
                ng = ng + 1
                if( ng > 48)
                    error("This should not happen: ng > 48")
                end
                sym_gshell[:,ng] = sn[:]
            end
        end
        #println("Ngs = ", Ngs, " ng = ", ng)
        Ngtot = Ngtot + ng
        shell_G_sym[Ngs] = Array{Int64,1}(undef,ng)
        for i = 1:ng
            found = false
            for j = ig:pw.gvec.Ng
                if is_done[j]
                    continue
                end
                found = ( sym_gshell[:,i] == G_crystal[:,j] )
                if found
                    is_done[j] = true
                    shell_G_sym[Ngs][i] = j
                    break
                end
            end
            if !found
                error("Lone G-vector")
            end
        end
    end

    println("Ngtot = ", Ngtot)
    println("pw.gvec.Ng = ", pw.gvec.Ng)
    println("Ngs = ", Ngs)

    Random.seed!(1234)


    dVol = pw.CellVolume/prod(pw.Ns)
    psiks = rand_BlochWavefunc(Ham)
    Rhoe = calc_rhoe(Ham, psiks)
    println("integ Rhoe = ", sum(Rhoe)*dVol)

    # assuming Nspin = 1
    RhoeG = R_to_G(pw, Rhoe[:,1])

    SMALL = 1e-10
    ft = sym_info.ft
    ft_ = zeros(3,Nsyms)
    non_symmorphic = Array{Bool,1}(undef,Nsyms)
    # convert fractional translations to Cartesian
    for isym = 1:Nsyms
        non_symmorphic[isym] = ( (abs(ft[1,isym]) >= SMALL) &&
                                 (abs(ft[2,isym]) >= SMALL) &&
                                 (abs(ft[3,isym]) >= SMALL) )
        #println("non_symmorphic ", isym, " ", non_symmorphic[isym])
        if non_symmorphic[isym]
            ft_[:,isym] = LatVecs[:,1]*ft[1,isym] + LatVecs[:,2]*ft[2,isym] + LatVecs[:,3]*ft[3,isym]
        end
    end

    sg = zeros(Float64,3)
    irot = zeros(Int64,48)

    for igl = 1:Ngs
        ng = length(shell_G_sym[igl])
        g0 = zeros(3,ng)
        is_done_shell = zeros(Bool,ng)
        for ig = 1:ng
            g0[:,ig] = pw.gvec.G[:,shell_G_sym[igl][ig]]
        end
        cryst_to_cart!(g0, LatVecs/(2*pi), -1)
        for ig = 1:ng
            
            if !is_done_shell[ig]
                
                rhosum = 0.0 + im*0.0

                for isym = 1:Nsyms
                    sg[:] = inv_s[:,1,isym]*g0[1,ig] + inv_s[:,2,isym]*g0[2,ig] + inv_s[:,3,isym]*g0[3,ig]
                
                    irot[isym] = 0
                    for isg = 1:ng
                        if all( abs.(sg[:] - g0[:,isg]) .< 1e-5 )
                            irot[isym] = isg
                            break
                        end
                    end

                    if (irot[isym] < 1) || (irot[isym] > ng)
                        error("Error in determining irot")
                    end

                    # isg is the index of rotated G-vector
                    isg = shell_G_sym[igl][irot[isym]]
                    ip = pw.gvec.idx_g2r[isg]

                    if non_symmorphic[isym]
                        arg = pw.gvec.G[1,isg]*ft_[1,isym] + pw.gvec.G[2,isg]*ft_[2,isym] + pw.gvec.G[3,isg]*ft_[3,isym]
                        fact = cos(arg) - im*sin(arg)
                        # assume Nspin == 1
                        rhosum = rhosum + RhoeG[ip]*fact
                    else
                        rhosum = rhosum + RhoeG[ip]
                    end
                
                end # isym

                rhosum = rhosum/Nsyms

                for isym = 1:Nsyms
                    isg = shell_G_sym[igl][irot[isym]]
                    ip = pw.gvec.idx_g2r[isg]
                    if non_symmorphic[isym]
                        arg = pw.gvec.G[1,isg]*ft_[1,isym] + pw.gvec.G[2,isg]*ft_[2,isym] + pw.gvec.G[3,isg]*ft_[3,isym]
                        fact = cos(arg) + im*sin(arg)
                        RhoeG[ip] = rhosum*fact
                    else
                        RhoeG[ip] = rhosum
                    end
                    is_done_shell[irot[isym]] = true
                end

            end  # is_done_shell
        end # ng
    end # ngl

    Rhoe_sym = real(G_to_R(pw, RhoeG))
    println("integ Rhoe_sym = ", sum(Rhoe_sym)*dVol)

    println("diff sum abs = ", sum(abs.(Rhoe_sym - Rhoe)))

end

function cryst_to_cart!( vec, trmat, iflag )
#subroutine cryst_to_cart (nvec, vec, trmat, iflag)
#  !-----------------------------------------------------------------------
#  !
#  !     This routine transforms the atomic positions or the k-point
#  !     components from crystallographic to cartesian coordinates 
#  !     ( iflag=1 ) and viceversa ( iflag=-1 ).
#  !     Output cartesian coordinates are stored in the input ('vec') array
#  !
#  !
  #
  #     Compute the cartesian coordinates of each vectors
  #     (atomic positions or k-points components)
  #
    nvec = size(vec)[2]
    vau = zeros(3)
    for nv=1:nvec
        if iflag==1
            for kpol=1:3
                vau[kpol] = trmat[kpol,1]*vec[1,nv] + trmat[kpol,2]*vec[2,nv] + trmat[kpol,3]*vec[3,nv]
            end
        else
            for kpol=1:3
                vau[kpol] = trmat[1,kpol]*vec[1,nv] + trmat[2,kpol]*vec[2,nv] + trmat[3,kpol]*vec[3,nv]
            end
        end
        for kpol=1:3
            vec[kpol,nv] = vau[kpol]
        end
    end
    return
end


main()
