struct RhoeSymmetrizer
    Ngs::Int64
    shell_G_sym::Array{Array{Int64,1},1}
    fts::Array{Float64,2}
end

# dummy RhoeSymmetrizer
function RhoeSymmetrizer()
    shell_G_sym = Array{Array{Int64,1},1}(undef,1)
    fts = zeros(1,1)
    return RhoeSymmetrizer(0, shell_G_sym, fts)
end

function RhoeSymmetrizer( atoms::Atoms, pw::PWGrid, sym_info::SymmetryInfo )

    LatVecs = atoms.LatVecs
    G = pw.gvec.G

    Nsyms = sym_info.Nsyms
    s = sym_info.s
    inv_s = sym_info.inv_s
    ft = sym_info.ft
    non_symmorphic = sym_info.non_symmorphic

    sym_gshell = Array{Int64,2}(undef,3,Nsyms)
    G_crystal = zeros(Int64,3,pw.gvec.Ng)
    is_done = zeros(Bool,pw.gvec.Ng)
    shell_G_sym = Array{Array{Int64,1},1}(undef,pw.gvec.Ng)

    for ig = 1:pw.gvec.Ng
        for ii = 1:3
            xx = ( LatVecs[1,ii]*G[1,ig] + LatVecs[2,ii]*G[2,ig] + LatVecs[3,ii]*G[3,ig] )/(2*pi)
            G_crystal[ii,ig] = round(Int64,xx)
        end
    end

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
            for ii = 1:3
                sn[ii] = s[ii,1,isym]*G_crystal[1,ig] + s[ii,2,isym]*G_crystal[2,ig] + s[ii,3,isym]*G_crystal[3,ig]
            end
            found = false
            for i = 1:ng
                found = ( (sn[1] == sym_gshell[1,i]) &&
                          (sn[2] == sym_gshell[2,i]) &&
                          (sn[3] == sym_gshell[3,i]) )
                if found
                    break
                end
            end
            if !found
                ng = ng + 1
                if( ng > Nsyms)
                    error("This should not happen: ng > Nsyms")
                end
                for ii = 1:3
                    sym_gshell[ii,ng] = sn[ii]
                end
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
                found = ( (sym_gshell[1,i] == G_crystal[1,j]) &&
                          (sym_gshell[2,i] == G_crystal[2,j]) &&
                          (sym_gshell[3,i] == G_crystal[3,j]) )
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

    fts = zeros(3,Nsyms)
    # convert fractional translations to Cartesian
    for isym = 1:Nsyms
        if non_symmorphic[isym]
            for ii = 1:3
                fts[ii,isym] = LatVecs[ii,1]*ft[1,isym] + LatVecs[ii,2]*ft[2,isym] + LatVecs[ii,3]*ft[3,isym]
            end
        end
    end

    return RhoeSymmetrizer(Ngs, shell_G_sym, fts)

end



function symmetrize_rhoe!(
    pw::PWGrid,
    sym_info::SymmetryInfo,
    rhoe_symmetrizer::RhoeSymmetrizer,
    Rhoe::Array{Float64,2}
)

    Ngs = rhoe_symmetrizer.Ngs
    shell_G_sym = rhoe_symmetrizer.shell_G_sym
    fts = rhoe_symmetrizer.fts

    LatVecs = pw.LatVecs
    Nsyms = sym_info.Nsyms
    inv_s = sym_info.inv_s
    ft = sym_info.ft
    non_symmorphic = sym_info.non_symmorphic

    Npoints = prod(pw.Ns)
    Nspin = size(Rhoe)[2]

    RhoeG = zeros(ComplexF64,Npoints,Nspin)
    for ispin = 1:Nspin
        @views RhoeG[:,ispin] = R_to_G(pw, Rhoe[:,ispin])
    end

    sg = zeros(Float64,3)
    irot = zeros(Int64,Nsyms)

    g0 = zeros(3,Nsyms)
    is_done_shell = zeros(Bool,Nsyms)

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

        g0 = trmat*g0

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
                        arg = pw.gvec.G[1,isg]*fts[1,isym] + pw.gvec.G[2,isg]*fts[2,isym] + pw.gvec.G[3,isg]*fts[3,isym]
                        fact = cos(arg) - im*sin(arg)
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
                        arg = pw.gvec.G[1,isg]*fts[1,isym] + pw.gvec.G[2,isg]*fts[2,isym] + pw.gvec.G[3,isg]*fts[3,isym]
                        fact = cos(arg) + im*sin(arg)
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
        @views Rhoe[:,ispin] = real(G_to_R(pw, RhoeG[:,ispin]))
    end

    return

end

function symmetrize_rhoe_v2!(
    pw::PWGrid,
    sym_info::SymmetryInfo,
    rhoe_symmetrizer::RhoeSymmetrizer,
    Rhoe::Array{Float64,2}
)

    Ngs = rhoe_symmetrizer.Ngs
    shell_G_sym = rhoe_symmetrizer.shell_G_sym
    fts = rhoe_symmetrizer.fts

    LatVecs = pw.LatVecs
    RecVecs = pw.RecVecs
    G = pw.gvec.G
    idx_g2r = pw.gvec.idx_g2r

    Nsyms = sym_info.Nsyms
    inv_s = sym_info.inv_s
    idx_inv_s = sym_info.idx_inv_s
    sname = sym_info.sname
    s = sym_info.s
    ft = sym_info.ft
    t_rev = sym_info.t_rev
    non_symmorphic = sym_info.non_symmorphic

    Npoints = prod(pw.Ns)
    Nspin_dens = size(Rhoe, 2)
    if Nspin_dens == 4
        Nspin_lsda = 1
    elseif Nspin_dens in [1,2]
        Nspin_lsda = Nspin_dens
    else
        error("Invalid value for Nspin_dens = $Nspin_dens")
    end

    RhoeG = zeros(ComplexF64, Npoints, Nspin_dens)
    for ispin in 1:Nspin_dens
        @views RhoeG[:,ispin] = R_to_G(pw, Rhoe[:,ispin])
    end

    sg = zeros(Float64, 3)
    irot = zeros(Int64, Nsyms)

    g0 = zeros(Float64, 3, Nsyms)
    is_done_shell = zeros(Bool, Nsyms)

    trmat = LatVecs'/(2*pi)

    rhosum = zeros(ComplexF64, 2)
    magsum = zeros(ComplexF64, 3)
    mag = zeros(ComplexF64, 3)
    magrot = zeros(ComplexF64, 3)
    
    non_symmorphic = zeros(Bool, 48)
    #
    # convert fractional translations to cartesian, in a0 units
    for isym in 1:Nsyms
        if non_symmorphic[isym]
            ft_[:,isym] = at[:,1]*ft[1,isym] + at[:,2]*ft(2,ns) + at[:,3]*ft[3,isym]
        end
    end
    
    if Nspin_dens == 4
        Nspin_lsda = 1
    elseif Nspin_dens in [1,2]
        Nspin_lsda = Nspin_dens
    else
        @error("Not valid number of Nspin_dens=$(Nspin_dens)")
    end
    #
    # scan shells of G-vectors
    #
    for igl in 1:Ngs
        
        # symmetrize: \rho_sym(G) = \sum_S rho(SG) for all G-vectors in the star
        #
        ng = length(shell_G_sym[igl])
        @assert ng >= 1
        for ig in 1:ng
            idx_shell = shell_G_sym[igl][ig]
            g0[1,ig] = G[1,idx_shell]
            g0[2,ig] = G[2,idx_shell]
            g0[3,ig] = G[3,idx_shell]
            is_done_shell[ig] = false
        end
        g0[:,:] = trmat*g0[:,:]
        #
        #  rotate G-vectors
        #
        for ig in 1:ng
            #
            if is_done_shell[ig]
                continue
            end
            
            fill!(rhosum, 0.0 + im*0.0)
            fill!(magsum, 0.0 + im*0.0)

            for isym in 1:Nsyms
                for ii in 1:3
                    sg[ii] = inv_s[ii,1,isym]*g0[1,ig] + inv_s[ii,2,isym]*g0[2,ig] + inv_s[ii,3,isym]*g0[3,ig]
                end
                #
                irot[isym] = 0
                for isg in 1:ng
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
                ip = idx_g2r[isg]
                #
                # non-spin-polarized case: component 1 is the charge
                # LSDA case: components 1,2 are spin-up and spin-down charge
                # non colinear case: component  1 is the charge density,
                #                    components 2,3,4 are the magnetization
                # non colinear case: components 2,3,4 are the magnetization
                #
                if Nspin_dens == 4
                    # bring magnetization to crystal axis
                    mag[:] = RhoeG[ip,2] * RecVecs[1,:] +
                             RhoeG[ip,3] * RecVecs[2,:] +
                             RhoeG[ip,4] * RecVecs[3,:]
                    # rotate and add magnetization
                    magrot[:] = inv_s[1,:,isym] * mag[1] +
                                inv_s[2,:,isym] * mag[2] +
                                inv_s[3,:,isym] * mag[3]
                    if sname[idx_inv_s][1:3] == "inv"
                        magrot[:] = -magrot[:]
                    end
                    if t_rev[idx_inv_s[isym]]
                        magrot[:] = -magrot[:]
                    end
                end

                if non_symmorphic[isym]
                    arg = G[1,isg]*fts[1,isym] + G[2,isg]*fts[2,isym] + G[3,isg]*fts[3,isym]
                    fact = cos(arg) - im*sin(arg)
                    for ispin in 1:Nspin_lsda
                        rhosum[ispin] += RhoeG[ip,ispin]*fact
                    end
                    if Nspin_dens == 4
                        magsum[:] .+= magrot[:] * fact
                    end
                else
                    for ispin in 1:Nspin_lsda
                        rhosum[ispin] += RhoeG[ip,ispin]
                    end
                    if Nspin_dens == 4
                        magsum[:] .+= magrot[:]
                    end
                end
                
            end # isym
            #
            for ispin in 1:Nspin_lsda
                rhosum[ispin] = rhosum[ispin]/Nsyms
            end

            for isym in 1:Nsyms
                isg = shell_G_sym[igl][irot[isym]]
                ip = idx_g2r[isg]
                if non_symmorphic[isym]
                    arg = G[1,isg]*fts[1,isym] + G[2,isg]*fts[2,isym] + G[3,isg]*fts[3,isym]
                    fact = cos(arg) + im*sin(arg)
                    for ispin in 1:Nspin_lsda
                        RhoeG[ip,ispin] = rhosum[ispin]*fact
                    end
                else
                    for ispin in 1:Nspin_lsda
                        RhoeG[ip,ispin] = rhosum[ispin]
                    end
                end
                is_done_shell[irot[isym]] = true
            end

        end # ng
    end # ngl

    for ispin in 1:Nspin_dens
        @views Rhoe[:,ispin] = real(G_to_R(pw, RhoeG[:,ispin]))
    end

    return
end
