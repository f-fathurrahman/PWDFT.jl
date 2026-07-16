mutable struct EXXVariables
    is_active::Bool
    exx_alpha::Float64
    #
    eps_occ::Float64
    eps_qdiv::Float64
    use_coulomb_vcut_ws::Bool
    use_coulomb_vcut_spheric::Bool
    use_regularization::Bool
    exxdiv::Float64
    gau_scrlen::Float64
    erf_scrlen::Float64
    erfc_scrlen::Float64
    yukawa::Float64
    x_gamma_extrapolation::Bool
    exxdiv_treatment::String
    grid_factor::Float64
    use_ace::Bool
    #
    ecutfock::Float64
    Nq1::Int64
    Nq2::Int64
    Nq3::Int64
    nkqs::Int64
    xkq::Matrix{Float64}
    nqs::Int64
    index_xkq::Matrix{Int64} # index_xkq(nks,nqs)
    index_xk::Vector{Int64} # index_xk(nkqs)
    index_sym::Vector{Int64} # index_sym(nkqs)
    rir::Matrix{Int64}
    Ns::Tuple{Int64,Int64,Int64}
    gvec::GVectors
    planfw::PWDFT.PLANFW_TYPE
    planbw::PWDFT.PLANBW_TYPE
    idx_gw2r::Vector{Vector{Int64}}
    x_occupation::Matrix{Float64}
    buff::Array{ComplexF64,3}
end


function _init_exx_grid(Ham, pwinput)
    
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nsyms = Ham.sym_info.Nsyms
    nq1 = pwinput.nqx1
    nq2 = pwinput.nqx2
    nq3 = pwinput.nqx3

    max_nk = Nkpt * min(48, 2*Nsyms)
    temp_index_xk = zeros(Int64, max_nk)
    temp_index_sym = zeros(Int64, max_nk)
    temp_index_ikq = zeros(Int64, max_nk)
    temp_xkq = zeros(Float64, 3, max_nk)
    
    pw = Ham.pw

    # find all k-points equivalent by symmetry to the points in the k-list
    temp_nkqs = 0
    xk_cryst = zeros(Float64, 3)
    sxk = zeros(Float64, 3)
    s = Ham.sym_info.s
    SMALL_Q = 1e-6
    noncolin = Ham.electrons.noncollinear
    domag = Ham.electrons.domag
    dxk = zeros(Float64, 3)

    invRecVecs = inv(pw.RecVecs)
    for isym in 1:Nsyms
        println("\nisym = ", isym)
        for ik in 1:Nkpt
            # go to crystalline coordinates
            xk_cryst[:] = pw.gvecw.kpoints.k[:,ik]
            println()
            @printf("%4d in bg   = [%18.10f %18.10f %18.10f]\n", ik, xk_cryst...)
            #cryst_to_cart!( 1, xk_cryst, pw.LatVecs/(2π), -1)
            xk_cryst[:] = invRecVecs*xk_cryst[:] # This is also can be used, but using inverse
            @printf("%4d in cart = [%18.10f %18.10f %18.10f]\n", ik, xk_cryst...)
            # rotate with this sym.op.
            sxk[:] = s[:,1,isym]*xk_cryst[1] + s[:,2,isym]*xk_cryst[2] + s[:,3,isym]*xk_cryst[3]
            # add sxk to the auxiliary list IF it is not already present
            xk_not_found = true
            # do-loop skipped the first time because temp_nkqs == 0
            for ikq in 1:temp_nkqs
                if xk_not_found
                    dxk[:] = sxk[:] - temp_xkq[:,ikq] - round.(Int64, sxk[:] - temp_xkq[:,ikq])
                    if all( abs.(dxk) .<= SMALL_Q )
                        xk_not_found = false
                    end
                end
            end
            #
            if xk_not_found
                temp_nkqs += 1
                temp_xkq[:,temp_nkqs] = sxk[:]
                temp_index_xk[temp_nkqs] = ik
                temp_index_sym[temp_nkqs] = isym
            end
            #
            sxk[:] = -sxk[:]
            xk_not_found = true
            for ikq in 1:temp_nkqs
                if xk_not_found
                    dxk[:] = sxk[:] - temp_xkq[:,ikq] - round.(Int64, sxk[:] - temp_xkq[:,ikq])
                    if all( abs.(dxk) .<= SMALL_Q )
                        xk_not_found = false
                    end
                end
            end
            if xk_not_found && !(noncolin && domag)
                temp_nkqs += 1
                temp_xkq[:,temp_nkqs] = sxk[:]
                temp_index_xk[temp_nkqs] = ik
                temp_index_sym[temp_nkqs] = -isym
            end
        end
    end
    println("temp_nkqs = ", temp_nkqs)
    
    # Find good q-point grid. Decrease the nqX until a good grid is found or
    # until it is 1 x 1 x 1 (always good)
    idx = 1
    sign_ = -1
    nqx = [nq1, nq2, nq3]
    nk1 = pw.gvecw.kpoints.mesh[1]
    nk2 = pw.gvecw.kpoints.mesh[2]
    nk3 = pw.gvecw.kpoints.mesh[3]
    #
    nqs = nq1*nq2*nq3
    nkqs, index_xkq = _exx_qgrid_init!(
        pw, nqs, nq1, nq2, nq3,
        max_nk, temp_nkqs, temp_xkq, temp_index_ikq, dxk
    )
    println("after exx_qgrid_init: nkqs = ", nkqs)
    println("dxk = ", dxk)
    #
    # Good q-point mesh
    is_good_q_pts = all( abs.(dxk) .<= SMALL_Q )
    #
    while !is_good_q_pts
        #
        # Try q-points around the input mesh, prioritizing smaller mesh
        nq1 = nqx[1] + idx * sign_
        nq2 = nqx[2] + idx * sign_
        nq3 = nqx[3] + idx * sign_
        #
        # Ensure no values smaller than 1
        if nq1 < 1 nq1 = 1 end
        if nq2 < 1 nq2 = 1 end
        if nq3 < 1 nq3 = 1 end
        #
        # Enforce nqX <= nkX. This is important for surfaces to keep the Z q-point 1.
        #
        if nq1 > nk1 nq1 = nk1 end
        if nq2 > nk2 nq2 = nk2 end
        if nq3 > nk3 nq3 = nk3 end
        #
        nqs = nq1 * nq2 * nq3
        #
        sign_ = -1 * sign_ # reduce of increase q grid
        #
        # Increase idx every other time sign is changed
        if sign_ < 0
            idx += 1
        end

        nkqs, index_xkq = exx_qgrid_init!(
            pw, nqs, nq1, nq2, nq3,
            max_nk, temp_nkqs, temp_xkq, temp_index_ikq, dxk
        )
        println("after exx_qgrid_init: nkqs = ", nkqs)
        println("dxk = ", dxk)
        #
        # Good q-point mesh
        is_good_q_pts = all( abs.(dxk) .<= SMALL_Q )
        if is_good_q_pts
            if idx > 1
                println("EXX: WARNING: q-point mesh has been updated!")
            end
            println("Found good qpoint mesh for EXX: q-point mesh: $nq1 $nq2 $nq3")
        end
    end

    nspin_lsda = 1 # XXX hardcoded

    xkq_collect = zeros(Float64, 3, nspin_lsda*nkqs)
    index_xk = zeros(Int64, nspin_lsda*nkqs)
    index_sym = zeros(Int64, nspin_lsda*nkqs)
    for ik in 1:nkqs
        ikq = temp_index_ikq[ik]
        xkq_collect[:,ik] = temp_xkq[:,ikq]
        index_xk[ik] = temp_index_xk[ikq]
        index_sym[ik] = temp_index_sym[ikq]
    end

    xkq_collect[:,:] = pw.RecVecs*xkq_collect[:,:]
    println("xkq_collect in Cartesian:")
    for ik in 1:nkqs
        @printf("[%8.5f %8.5f %8.5f] %4d %4d\n",
            xkq_collect[1,ik], xkq_collect[2,ik], xkq_collect[3,ik], 
            index_xk[ik], index_sym[ik])
    end

    _exx_grid_check(
        pw, Ham.sym_info,
        nq1, nq2, nq3, index_xk, index_xkq, index_sym
    )

    qnrm = 0.0
    for iq in 1:nkqs, ik in 1:Nkpt
        qnrm = max(qnrm, sqrt( sum( (pw.gvecw.kpoints.k[:,ik] - xkq_collect[:,iq]).^2) ))
    end
    println("qnrm = ", qnrm)

    return nq1, nq2, nq3, nkqs, xkq_collect, index_xkq, index_xk, index_sym
end



function _calc_fft_grid_size(LatVecs, ecut)
    alat = sqrt(LatVecs[1,1]^2 + LatVecs[2,1]^2 + LatVecs[3,1]^2)
    tpiba  = 2π/alat
    tpiba2 = tpiba^2
    gcutm = ecut/tpiba2

    LatVecsLen1 = norm(LatVecs[:,1])
    LatVecsLen2 = norm(LatVecs[:,2])
    LatVecsLen3 = norm(LatVecs[:,3])
    #=
    # This is the old
    Ns1 = 2*round( Int64, sqrt(ecut/2)*LatVecsLen1/pi ) + 1
    Ns2 = 2*round( Int64, sqrt(ecut/2)*LatVecsLen2/pi ) + 1
    Ns3 = 2*round( Int64, sqrt(ecut/2)*LatVecsLen3/pi ) + 1
    #
    =#
    Ns1 = round(Int64, sqrt(gcutm)*LatVecsLen1/pi ) + 1
    Ns2 = round(Int64, sqrt(gcutm)*LatVecsLen2/pi ) + 1
    Ns3 = round(Int64, sqrt(gcutm)*LatVecsLen3/pi ) + 1

    Ns1 = PWDFT.good_fft_order(Ns1)
    Ns2 = PWDFT.good_fft_order(Ns2)
    Ns3 = PWDFT.good_fft_order(Ns3)
    #
    return Ns1, Ns2, Ns3
end

function _init_g_to_r_map(gvec, kpoints, ecutwfc)
    G = gvec.G
    Ng = gvec.Ng
    idx_g2r = gvec.idx_g2r
    #
    kpts = kpoints.k
    Nkpt = kpoints.Nkpt
    #
    Gk2 = zeros(Float64, Ng)
    Gk = zeros(Float64, 3)
    idx_gw2g = Array{Array{Int64,1},1}(undef, Nkpt)
    idx_gw2r = Array{Array{Int64,1},1}(undef, Nkpt)
    #
    for ik in 1:Nkpt
        for ig in 1:Ng
            Gk[1] = G[1,ig] + kpts[1,ik]
            Gk[2] = G[2,ig] + kpts[2,ik]
            Gk[3] = G[3,ig] + kpts[3,ik]
            Gk2[ig] = Gk[1]^2 + Gk[2]^2 + Gk[3]^2
        end
        idx_gw2g[ik] = findall( 0.5*Gk2 .<= ecutwfc )
        idx_gw2r[ik] = idx_g2r[idx_gw2g[ik]]
    end
    #XXX we should only need idx_gw2r
    #XXX I think we can use idx_gw2g from pw.gvecw because they should be the same
    return idx_gw2r
end


function EXXVariables(Ham, pwinput)
        
    Nq1, Nq2, Nq3, nkqs, xkq, index_xkq, index_xk, index_sym = _init_exx_grid(Ham, pwinput)
    nqs = Nq1*Nq2*Nq3
    
    pw = Ham.pw
    ecutfock = pwinput.ecutfock
    kpoints = pw.gvecw.kpoints
    Nkpt = kpoints.Nkpt
    wk = kpoints.wk
    Ngw = pw.gvecw.Ngw
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc

    Ns = _calc_fft_grid_size(pw.LatVecs, ecutfock)
    Npoints_exx = prod(Ns)
    gvec = PWDFT.init_gvec(Ns, pw.RecVecs, ecutfock)
    planfw = plan_fft!(zeros(ComplexF64, Ns), flags = FFTW.MEASURE)
    planbw = plan_ifft!(zeros(ComplexF64, Ns), flags = FFTW.MEASURE)

    idx_gw2r = _init_g_to_r_map(gvec, kpoints, pw.ecutwfc)

    rir = exx_set_symm(Ham.sym_info, Ns)

    # set occupations of wavefunctions used in the calculation of exchange term
    x_occupation = zeros(Float64, Nstates, Nkpt)
    # do not include orbitals with small occupations (?)
    SMALL_OCC = 1e-8 # XXX Harcoded
    for ik in 1:Nkpt
        if abs(wk[ik]) > SMALL_OCC # XXX why only check wk?
            x_occupation[1:Nstates,ik] = Focc[1:Nstates,ik]*0.5  # Set to [0,1]
            #FIXME: Need to account for spin?
        else
            x_occupation[1:Nstates,ik] = 0.0
        end
    end

    buff = zeros(ComplexF64, Npoints_exx, Nstates, nkqs)

    is_active = false
    # Some variables from pwinput
    x_gamma_extrapolation = pwinput.x_gamma_extrapolation
    exxdiv_treatment = pwinput.exxdiv_treatment
    use_ace = false

    eps_occ = 1e-8
    eps_qdiv = 1e-8
    use_coulomb_vcut_ws =  false
    use_coulomb_vcut_spheric =  false
    
    gau_scrlen = 0.0
    erf_scrlen = 0.0
    erfc_scrlen = 0.0
    yukawa = 0.0
    
    grid_factor = 1.0
    if x_gamma_extrapolation
        grid_factor = 8/7
    end

    exxdiv_treatment = pwinput.exxdiv_treatment
    use_regularization = !(exxdiv_treatment == "none")
    # exxdiv has energy unit
    exxdiv = calc_exx_divergence(
        pw, Nq1, Nq2, Nq3,
        x_gamma_extrapolation,
        erf_scrlen, erfc_scrlen, yukawa, grid_factor, nqs,
        use_regularization = use_regularization
    )

    exx_alpha = 1.0

    return EXXVariables(
        is_active,
        exx_alpha,
        eps_occ, eps_qdiv,
        use_coulomb_vcut_ws,
        use_coulomb_vcut_spheric,
        use_regularization,
        exxdiv,
        gau_scrlen, erf_scrlen, erfc_scrlen,
        yukawa,
        x_gamma_extrapolation,
        exxdiv_treatment,
        grid_factor,
        use_ace,
        ecutfock,
        Nq1, Nq2, Nq3,
        nkqs,
        xkq,
        nqs,
        index_xkq,
        index_xk,
        index_sym,
        rir,
        Ns, gvec, planfw, planbw, idx_gw2r,
        x_occupation,
        buff
    )

end


