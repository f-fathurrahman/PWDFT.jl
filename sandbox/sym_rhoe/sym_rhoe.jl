function init_sym_rhoe( Ham::Hamiltonian )

    atoms = Ham.atoms
    pw = Ham.pw
    sym_info = Ham.sym_info
    LatVecs = atoms.LatVecs
    G = pw.gvec.G

    sym_gshell = Array{Int64,2}(undef,3,48)
    G_crystal = zeros(Int64,3,pw.gvec.Ng)
    is_done = zeros(Bool,pw.gvec.Ng)
    shell_G_sym = Array{Array{Int64,1},1}(undef,pw.gvec.Ng)

    for ig = 1:pw.gvec.Ng
        xx = ( LatVecs[1,:]*G[1,ig] + LatVecs[2,:]*G[2,ig] + LatVecs[3,:]*G[3,ig] )/(2*pi)
        G_crystal[:,ig] = round.(Int64,xx)
        is_done[ig] = false
    end

    Nsyms = sym_info.Nsyms
    s = sym_info.s
    inv_s = sym_info.inv_s
    ft = sym_info.ft
    non_symmorphic = sym_info.non_symmorphic

    Ngs = 0
    sn = zeros(Int64,3)

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
                ng = ng + 1
                if( ng > 48)
                    error("This should not happen: ng > 48")
                end
                sym_gshell[:,ng] = sn[:]
            end
        end

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

    ft_ = zeros(3,Nsyms)
    # convert fractional translations to Cartesian
    alat = norm(LatVecs[:,1])
    for isym = 1:Nsyms
        if non_symmorphic[isym]
            for ii = 1:3
                ft_[ii,isym] = LatVecs[ii,1]*ft[1,isym] + LatVecs[ii,2]*ft[2,isym] + LatVecs[ii,3]*ft[3,isym]
            end
        end
    end

    return Ngs, shell_G_sym, ft_

end


function sym_rhoe!( Ham, Rhoe, Ngs, shell_G_sym, ft_ )

    pw = Ham.pw
    sym_info = Ham.sym_info

    LatVecs = pw.LatVecs
    Nsyms = sym_info.Nsyms
    inv_s = sym_info.inv_s
    ft = sym_info.ft
    non_symmorphic = sym_info.non_symmorphic

    Npoints = prod(pw.Ns)
    Nspin = Ham.electrons.Nspin

    RhoeG = zeros(ComplexF64,Npoints,Nspin)
    for ispin = 1:Nspin
        RhoeG[:,ispin] = R_to_G(pw, Rhoe[:,ispin])
    end

    sg = zeros(Float64,3)
    irot = zeros(Int64,48)

    g0 = zeros(3,48)
    is_done_shell = zeros(Bool,48)

    trmat = LatVecs'/(2*pi)

    rhosum = zeros(ComplexF64,Nspin)

    for igl = 1:Ngs
        
        ng = length(shell_G_sym[igl])

        for ig = 1:ng
            g0[1,ig] = pw.gvec.G[1,shell_G_sym[igl][ig]]
            g0[2,ig] = pw.gvec.G[2,shell_G_sym[igl][ig]]
            g0[3,ig] = pw.gvec.G[3,shell_G_sym[igl][ig]]
            is_done_shell[ig] = false
        end
        
        #_cryst_to_cart!(g0, LatVecs/(2*pi), -1)
        g0 = trmat*g0
        #for ig = 1:ng
        #    for ii=1:3
        #        g0[ii,ig] = (LatVecs[1,ii]*g0[1,ig] + LatVecs[2,ii]*g0[2,ig] + LatVecs[3,ii]*g0[3,ig])/(2*pi)
        #    end
        #end

        for ig = 1:ng
            
            if !is_done_shell[ig]
                
                for ispin = 1:Nspin
                    rhosum[ispin] = 0.0 + im*0.0
                end

                for isym = 1:Nsyms
                    for ii = 1:3
                        sg[ii] = inv_s[ii,1,isym]*g0[1,ig] + inv_s[ii,2,isym]*g0[2,ig] + inv_s[ii,3,isym]*g0[3,ig]
                    end

                    irot[isym] = 0

                    for isg = 1:ng
                        if (abs(sg[1] - g0[1,isg]) < 1e-5) &&
                           (abs(sg[2] - g0[2,isg]) < 1e-5) &&
                           (abs(sg[3] - g0[3,isg]) < 1e-5)
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
                        fact = cos(arg) + im*sin(arg)
                        for ispin = 1:Nspin
                            rhosum[ispin] = rhosum[ispin] + RhoeG[ip,ispin]*fact
                        end
                    else
                        for ispin = 1:Nspin
                            rhosum[ispin] = rhosum[ispin] + RhoeG[ip,ispin]
                        end
                    end
                
                end # isym

                for ispin = 1:Nspin
                    rhosum[ispin] = rhosum[ispin]/Nsyms
                end

                for isym = 1:Nsyms
                    isg = shell_G_sym[igl][irot[isym]]
                    ip = pw.gvec.idx_g2r[isg]
                    if non_symmorphic[isym]
                        arg = pw.gvec.G[1,isg]*ft_[1,isym] + pw.gvec.G[2,isg]*ft_[2,isym] + pw.gvec.G[3,isg]*ft_[3,isym]
                        fact = cos(arg) - im*sin(arg)
                        for ispin = 1:Nspin
                            RhoeG[ip,ispin] = rhosum[ispin]*fact
                        end
                    else
                        for ispin = 1:Nspin
                            RhoeG[ip,ispin] = rhosum[ispin]
                        end
                    end
                    is_done_shell[irot[isym]] = true
                end

            end  # is_done_shell
        end # ng
    end # ngl

    for ispin = 1:Nspin
        Rhoe[:,ispin] = real(G_to_R(pw, RhoeG[:,ispin]))
    end

    return

end

function _cryst_to_cart!( vec, trmat, iflag )
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