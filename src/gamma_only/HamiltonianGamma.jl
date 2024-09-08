mutable struct HamiltonianGamma
    pw::PWGridGamma
    potentials::Potentials
    energies::Energies
    rhoe::Array{Float64,2} # spin dependent, spin is not implemented yet
    electrons::Electrons
    atoms::Atoms
    pspots::Array{PsPot_GTH,1}
    pspotNL::PsPotNLGamma
    xcfunc::String
    xc_calc::AbstractXCCalculator
    ispin::Int64 # No spin index is implemented yet
end

function HamiltonianGamma(
    atoms::Atoms, pspfiles::Array{String,1},
    ecutwfc::Float64 ;
    Nspin=1,
    Ns_=(0,0,0),
    xcfunc="VWN",
    use_xc_internal=false,
    extra_states=0,
)

    # Initialize plane wave grids
    pw = PWGridGamma( ecutwfc, atoms.LatVecs, Ns_=Ns_ )

    Nspecies = atoms.Nspecies
    if Nspecies != size(pspfiles)[1]
        error( @sprintf("Length of pspfiles is not equal to %d\n", Nspecies) )
    end

    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    idx_g2rm = pw.gvec.idx_g2rm

    strf = calc_strfact( atoms.positions, atoms.Nspecies, atoms.atm2species, pw.gvec.G )

    #
    # Initialize pseudopotentials and local potentials
    #
    pspots = Array{PsPot_GTH}(undef,Nspecies)
    G2_shells = pw.gvec.G2_shells
    Ngl = length(G2_shells)
    Vgl = zeros(Float64, Ngl)
    idx_g2shells = pw.gvec.idx_g2shells
    V_Ps_loc = zeros(Float64, Npoints)
    Vg = zeros(ComplexF64, Npoints)

    for isp = 1:Nspecies
        pspots[isp] = PsPot_GTH( pspfiles[isp] )
        psp = pspots[isp]
        #
        for igl = 1:Ngl
            Vgl[igl] = eval_Vloc_G( psp, G2_shells[igl] )
        end
        #
        Vg[1] = strf[1,isp]*Vgl[1] # G=(0,0,0)
        for ig = 2:Ng
            #
            igl = idx_g2shells[ig]
            #
            ip = idx_g2r[ig]
            Vg[ip] = strf[ig,isp] * Vgl[igl]
            #
            ipm = idx_g2rm[ig]
            Vg[ipm] = conj(strf[ig,isp]) * Vgl[igl]
        end
        #
        V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(pw, Vg) ) * Npoints / CellVolume
    end

    # other potential terms are set to zero
    V_Hartree = zeros( Float64, Npoints )
    V_xc = zeros( Float64, Npoints, Nspin )
    V_loc_tot = zeros( Float64, Npoints, Nspin )
    potentials = Potentials( V_Ps_loc, V_Hartree, V_xc, V_loc_tot )
    #
    energies = Energies()
    #
    rhoe = zeros( Float64, Npoints, Nspin )

    electrons = Electrons( atoms, pspots, Nspin=Nspin, Nkpt=1,
                           Nstates_empty=extra_states )

    # NL pseudopotentials
    pspotNL = PsPotNLGamma( atoms, pw, pspots )

    atoms.Zvals = get_Zvals( pspots )

    ispin = 1

    if use_xc_internal
        xc_calc = XCCalculator()
    else
        # Using Libxc is the default
        xc_calc = LibxcXCCalculator()
    end

    return HamiltonianGamma(
        pw, potentials, energies, rhoe,
        electrons, atoms,
        pspots, pspotNL, xcfunc, xc_calc, ispin )
end


function update!( Ham::HamiltonianGamma, Rhoe::Array{Float64,2} )
    
    Nspin = size(Rhoe,2)
    
    # Copy Rhoe
    Ham.rhoe[:] = Rhoe[:]

    Poisson_solve!( Ham.pw, Ham.rhoe, Ham.potentials.Hartree )
    
    if Ham.xcfunc == "PBE"
        Ham.potentials.XC[:] = calc_Vxc_PBE( Ham.xc_calc, Ham.pw, Rhoe )
    else  # VWN is the default
        Ham.potentials.XC[:] = calc_Vxc_VWN( Ham.xc_calc, Rhoe )
    end
    
    Npoints = prod(Ham.pw.Ns)
    for ispin in 1:Nspin, ip in 1:Npoints
        Ham.potentials.Total[ip,ispin] = Ham.potentials.Ps_loc[ip] + Ham.potentials.Hartree[ip] +
                                         Ham.potentials.XC[ip,ispin]  
    end
    
    return
end


import Base: print
function print( io::IO, Ham::HamiltonianGamma; header=true )
    if header
        @printf("\n")
        @printf("                           --------------------------\n")
        @printf("                           Hamiltonian (Γ-point only)\n")
        @printf("                           --------------------------\n")
        @printf("\n")
    end
    @printf(io, "size (MiB) = %18.5f\n", Base.summarysize(Ham)/1024/1024)
    println(io, "")
    println(io, "xcfunc     = ", Ham.xcfunc)
    println(io, "xc_calc    = ", Ham.xc_calc)
    println(io, "")
    print(io, Ham.atoms)
    #print(io, Ham.pw)
    print(io, Ham.electrons)
    for isp = 1:Ham.atoms.Nspecies
        print(io, Ham.pspots[isp])
    end
end
print( Ham::HamiltonianGamma; header=true ) = print( stdout, Ham, header=header )