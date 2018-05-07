include("Energies.jl")

include("Potentials.jl")

include("PsPotNL.jl")

mutable struct PWHamiltonian
    pw::PWGrid
    potentials::Potentials
    energies::Energies
    rhoe::Array{Float64,1}
    electrons::Electrons
    atoms::Atoms
    pspots::Array{PsPot_GTH,1}
    pspotNL::PsPotNL
    xcfunc::String
    ik::Int64   # current kpoint index
end


function PWHamiltonian( atoms::Atoms, pspfiles::Array{String,1},
                        ecutwfc::Float64 ;
                        meshk = [1,1,1], shiftk=[0,0,0],
                        xcfunc = "VWN", verbose=false, extra_states=0 )

    # kpoints
    kpoints = KPoints( atoms, meshk, shiftk )

    # Initialize plane wave grids
    pw = PWGrid( ecutwfc, atoms.LatVecs, kpoints=kpoints )
    if verbose
        println(pw)
    end

    Nspecies = atoms.Nspecies
    if Nspecies != size(pspfiles)[1]
        @printf("ERROR length of pspfiles is not equal to %d\n", Nspecies)
        exit()
    end

    Npoints = prod(pw.Ns)
    Ω = pw.Ω
    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r

    strf = calc_strfact( atoms, pw )

    #
    # Initialize pseudopotentials and local potentials
    #
    Vg = zeros(Complex128, Npoints)
    V_Ps_loc = zeros(Float64, Npoints)

    Pspots = Array{PsPot_GTH}(Nspecies)

    for isp = 1:Nspecies
        Pspots[isp] = PsPot_GTH( pspfiles[isp] )
        psp = Pspots[isp]
        if verbose
            println(psp)
        end
        for ig = 1:Ng
            ip = idx_g2r[ig]
            Vg[ip] = strf[ig,isp] * eval_Vloc_G( psp, G2[ig] ) / Ω
        end
        #
        V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(pw, Vg) ) * Npoints
    end

    # other potential terms are set to zero
    V_Hartree = zeros( Float64, Npoints )
    V_XC = zeros( Float64, Npoints )
    potentials = Potentials( V_Ps_loc, V_Hartree, V_XC )
    #
    energies = Energies()
    #
    rhoe = zeros( Float64, Npoints )

    electrons = Electrons( atoms, Pspots, Nkpt=kpoints.Nkpt, Nstates_empty=extra_states )
    if verbose
        println(electrons)
    end

    # NL pseudopotentials
    pspotNL = PsPotNL( pw, atoms, Pspots, kpoints, check_norm=false )

    atoms.Zvals = get_Zvals( Pspots )

    ik = 1
    return PWHamiltonian( pw, potentials, energies, rhoe,
                          electrons, atoms, Pspots, pspotNL, xcfunc, ik )
end


#
# No pspfiles given. Use Coulomb potential (all electrons)
#
function PWHamiltonian( atoms::Atoms, ecutwfc::Float64;
                        meshk = [1,1,1], shiftk=[0,0,0],    
                        xcfunc="VWN", verbose=false )
    
    # kpoints
    kpoints = KPoints( atoms, meshk, shiftk )

    # Initialize plane wave grids
    pw = PWGrid( ecutwfc, atoms.LatVecs, kpoints=kpoints )
    if verbose
        println(pw)
    end

    Nspecies = atoms.Nspecies

    Npoints = prod(pw.Ns)

    # dummy pspots
    Pspots = Array{PsPot_GTH}(1)
    Pspots[1] = PsPot_GTH()

    #
    # Coulomb potential as local potential
    #
    Zvals = get_Zatoms(atoms)
    strf = calc_strfact( atoms, pw )
    V_Ps_loc = init_V_coulomb_G( pw, strf, Zvals )

    # other potentials are set to zero
    V_Hartree = zeros( Float64, Npoints )
    V_XC = zeros( Float64, Npoints )
    potentials = Potentials( V_Ps_loc, V_Hartree, V_XC )
    #
    energies = Energies()
    #
    rhoe = zeros( Float64, Npoints )

    electrons = Electrons( atoms, Zvals, Nkpt=kpoints.Nkpt )
    if verbose
        println(electrons)
    end

    rhoe = zeros(Float64,Npoints)

    pspotNL = PsPotNL()

    Zatoms = get_Zatoms( atoms )
    atoms.Zvals = Zatoms
    
    ik = 1
    return PWHamiltonian( pw, potentials, energies, rhoe,
                          electrons, atoms, Pspots, pspotNL, xcfunc, ik )
end



include("op_K.jl")
include("op_V_loc.jl")
include("op_V_Ps_nloc.jl")
include("op_H.jl")

include("Poisson_solve.jl")


"""
Given rhoe in real space, update Ham.rhoe, Hartree and XC potentials.
"""
function update!(Ham::PWHamiltonian, rhoe::Array{Float64,1})
    Ham.rhoe = rhoe
    Ham.potentials.Hartree = real( G_to_R( Ham.pw, Poisson_solve(Ham.pw, rhoe) ) )
    if Ham.xcfunc == "PBE"
        Ham.potentials.XC = calc_Vxc_PBE( Ham.pw, rhoe )
    else  # VWN is the default
        Ham.potentials.XC = calc_Vxc_VWN( rhoe )
    end
end


function update!( Ham::PWHamiltonian, atoms::Atoms,
    strf::Array{Complex128,2}, pspfiles::Array{String,1} )

    Nspecies = atoms.Nspecies
    if Nspecies != size(pspfiles)[1]
        @printf("ERROR length of pspfiles is not equal to %d\n", Nspecies)
        exit()
    end

    pw = Ham.pw

    Npoints = prod(pw.Ns)
    Ω = pw.Ω
    G2 = pw.gvec.G2

    Vg = zeros(Complex128, Npoints)
    V_Ps_loc = zeros(Float64, Npoints)

    pw = Ham.pw
    for isp = 1:Nspecies
        psp = PsPot_GTH( pspfiles[isp] )
        println(psp)
        for ig = 1:Npoints
            Vg[ig] = strf[ig,isp] * eval_Vloc_G( psp, G2[ig], Ω )
        end
        #
        V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(pw, Vg) ) * Npoints
    end

    # Update
    Ham.potentials.Ps_loc = V_Ps_loc
    return

end
