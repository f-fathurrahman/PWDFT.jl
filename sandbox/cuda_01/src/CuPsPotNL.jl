struct CuPsPotNL
    NbetaNL::Int64
    prj2beta::Array{Int64,4}
    betaNL::Array{CuArray{ComplexF64,2},1}
    #
    # Several copy used in calc_E_Ps_nloc
    #
    prj2beta_gpu::CuArray{Int64,4}
    # Atom
    atm2species::CuArray{Int64,1}
    # PsPot_GTH (array)
    psp_lmax::CuArray{Int64,1}
    psp_Nproj_l::CuArray{Int64,2}
    psp_h::CuArray{Float64,4}
    # from KPoints
    wk::CuArray{Float64,1}
    # From Electrons
    Focc::CuArray{Float64,2}
end

function CuPsPotNL(
    atoms::Atoms,
    pspots::Array{PsPot_GTH,1},
    kpoints::KPoints,
    electrons::CuElectrons,
    pspNL::PsPotNL,
)
    Nkpt = length( pspNL.betaNL )
    betaNL = Array{CuArray{ComplexF64,2},1}(undef, Nkpt)
    for ik = 1:Nkpt
        betaNL[ik] = CuArray( pspNL.betaNL[ik] )
    end

    prj2beta_gpu = CuArray( pspNL.prj2beta )

    atm2species = CuArray( atoms.atm2species )

    Nspecies = atoms.Nspecies
    #
    psp_lmax_cpu = zeros( Int64, Nspecies )
    psp_Nproj_l_cpu = zeros( Int64, 4, Nspecies )
    psp_h_cpu = zeros( Float64, 4,3,3, Nspecies )

    for isp = 1:Nspecies
        psp_lmax_cpu[isp] = pspots[isp].lmax
        psp_Nproj_l_cpu[:,isp] = pspots[isp].Nproj_l
        psp_h_cpu[:,:,:,isp] = pspots[isp].h
    end

    psp_lmax = CuArray( psp_lmax_cpu )
    psp_Nproj_l = CuArray( psp_Nproj_l_cpu )
    psp_h = CuArray( psp_h_cpu )

    wk = CuArray( kpoints.wk )

    Focc = copy( electrons.Focc_gpu )

    return CuPsPotNL(
        pspNL.NbetaNL,
        pspNL.prj2beta,
        betaNL,
        prj2beta_gpu,
        atm2species,
        psp_lmax,
        psp_Nproj_l,
        psp_h,
        wk,
        Focc
    )
end


import PWDFT: calc_betaNL_psi

function calc_betaNL_psi(
    ik::Int64,
    betaNL::Array{CuArray{ComplexF64,2},1},
    psi::CuArray{ComplexF64,2}
)

    Nstates = size(psi)[2]
    NbetaNL = size(betaNL[1],2)

    betaNL_psi = CuArrays.zeros( ComplexF64, Nstates, NbetaNL )
    betaNL_psi[:,:] = conj( psi' * betaNL[ik] )
    return betaNL_psi
end

function calc_betaNL_psi(
    ik::Int64,
    betaNL::Array{CuArray{ComplexF64,2},1},
    psi::Array{ComplexF64,1}
)
    NbetaNL = size(betaNL[1],2)

    betaNL_psi = CuArrays.zeros( ComplexF64, NbetaNL )
    betaNL_psi[:] = conj( psi' * betaNL[ik] )
    return betaNL_psi
end

