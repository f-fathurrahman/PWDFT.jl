mutable struct CuHamiltonian
    pw::CuPWGrid
    potentials::CuPotentials
    energies::Energies
    rhoe::CuArray{Float64,2}
    electrons::CuElectrons
    atoms::Atoms
    sym_info::SymmetryInfo
    rhoe_symmetrizer::RhoeSymmetrizer
    pspots::Array{PsPot_GTH,1}
    pspotNL::CuPsPotNL
    xcfunc::String
    xc_calc::AbstractXCCalculator
    ik::Int64
    ispin::Int64
end


function CuHamiltonian( atoms::Atoms, pspfiles::Array{String,1},
                        ecutwfc::Float64 ;
                        Nspin=1,
                        meshk=[1,1,1], shiftk=[0,0,0], time_reversal=true,
                        Ns_=(0,0,0),
                        kpoints=nothing,
                        kpts_str="",
                        xcfunc="VWN",
                        use_xc_internal=true,
                        extra_states=0,
                        use_symmetry=true )

    if use_symmetry == false
        sym_info = SymmetryInfo()
    else
        sym_info = SymmetryInfo(atoms)
    end

    if sym_info.Nsyms != 1
        println("-------------------------------------------------------------------")
        println("Your system has high symmetry: Nsyms = ", sym_info.Nsyms)
        println("Please break the symmetry of your system or set use_symmetry=false.")
        println("-------------------------------------------------------------------")
        error("We only support Nsyms = 1 for CUDA-enabled calculations")
    end

    # kpoints
    if kpoints == nothing
        if kpts_str == ""
            # automatic generation of kpoints
            kpoints = KPoints( atoms, meshk, shiftk, sym_info.s, time_reversal=time_reversal )
        else
            # use the given kpoints
            kpoints = kpoints_from_string( atoms, kpts_str )
        end
    else
        @assert typeof(kpoints) == KPoints
    end

    # Initialize plane wave grids
    pw = CuPWGrid( ecutwfc, atoms.LatVecs, kpoints=kpoints, Ns_=Ns_ )

    # XXX Need to calculate PWGrid again because it is needed in several places
    pw_ = PWGrid( ecutwfc, atoms.LatVecs, kpoints=kpoints, Ns_=Ns_ )

    Nspecies = atoms.Nspecies
    if Nspecies != size(pspfiles)[1]
        error( @sprintf("Length of pspfiles is not equal to %d\n", Nspecies) )
    end

    Npoints = prod(pw.Ns)
    CellVolume = pw_.CellVolume
    G2 = pw_.gvec.G2
    Ng = pw_.gvec.Ng
    idx_g2r = pw_.gvec.idx_g2r

    strf = calc_strfact( atoms, pw_ )

    #
    # Initialize pseudopotentials and local potentials
    #
    Pspots = Array{PsPot_GTH}(undef,Nspecies)

    G2_shells = pw_.gvec.G2_shells
    idx_g2shells = pw_.gvec.idx_g2shells

    Ngl = length(G2_shells)
    Vgl = zeros(Float64, Ngl)
    V_Ps_loc = zeros(Float64, Npoints)
    Vg = zeros(ComplexF64, Npoints)

    # TODO: using kernel (CUDAnative.erf is available)
    for isp = 1:Nspecies
        Pspots[isp] = PsPot_GTH( pspfiles[isp] )
        psp = Pspots[isp]
        #
        for igl = 1:Ngl
            Vgl[igl] = eval_Vloc_G( psp, G2_shells[igl] )
        end
        #
        for ig = 1:Ng
            ip = idx_g2r[ig]
            igl = idx_g2shells[ig]
            Vg[ip] = strf[ig,isp] * Vgl[igl]
        end
        #
        V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(pw_, Vg) ) * Npoints / CellVolume
    end

    # other potential terms are set to zero
    V_Hartree = CuArrays.zeros( Float64, Npoints )
    V_xc = CuArrays.zeros( Float64, Npoints, Nspin )
    V_loc_tot = CuArrays.zeros( Float64, Npoints, Nspin )
    potentials = CuPotentials( CuArray(V_Ps_loc), V_Hartree, V_xc, V_loc_tot )
    #
    energies = Energies()
    #
    rhoe = zeros( Float64, Npoints, Nspin )

    electrons_ = Electrons( atoms, Pspots, Nspin=Nspin, Nkpt=kpoints.Nkpt,
                            Nstates_empty=extra_states )
    electrons = CuElectrons( electrons_ )

    # NL pseudopotentials
    pspotNL_ = PsPotNL( atoms, pw_, Pspots, check_norm=false ) # XXX should be done on GPU?
    
    pspotNL = CuPsPotNL( atoms, Pspots, kpoints, electrons, pspotNL_ ) # copy

    atoms.Zvals = get_Zvals( Pspots )

    ik = 1
    ispin = 1

    if sym_info.Nsyms > 1
        rhoe_symmetrizer = RhoeSymmetrizer( atoms, pw_, sym_info )
    else
        rhoe_symmetrizer = RhoeSymmetrizer() # dummy rhoe_symmetrizer
    end

    @assert use_xc_internal == true
    xc_calc = XCCalculator()

    return CuHamiltonian( pw, potentials, energies, rhoe,
                          electrons, atoms, sym_info, rhoe_symmetrizer,
                          Pspots, pspotNL, xcfunc, xc_calc, ik, ispin )
end


import PWDFT: update!

function update!( Ham::CuHamiltonian, rhoe::CuArray{Float64,2} )
    
    Nspin = size(rhoe)[2]
    
    Ham.rhoe = rhoe[:,:]
    if Nspin == 2
        Rhoe_total = Ham.rhoe[:,1] + Ham.rhoe[:,2] # Nspin is 2
    else
        Rhoe_total = Ham.rhoe[:,1]
    end
    Ham.potentials.Hartree = real( G_to_R( Ham.pw, Poisson_solve(Ham.pw, Rhoe_total) ) )
    
    if Ham.xcfunc == "PBE"
        Ham.potentials.XC = calc_Vxc_PBE( Ham.xc_calc, Ham.pw, rhoe )
    else
        Ham.potentials.XC = calc_Vxc_VWN( Ham.xc_calc, rhoe )
    end
    
    Npoints = prod(Ham.pw.Ns)
    
    for ispin = 1:Nspin
        Ham.potentials.Total[:,ispin] = Ham.potentials.Ps_loc[:] + Ham.potentials.Hartree[:] +
                                        Ham.potentials.XC[:,ispin]  
    end
    
    return
end