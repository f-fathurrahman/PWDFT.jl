mutable struct Hamiltonian
    pw::PWGrid
    potentials::Potentials
    energies::Energies
    rhoe::Array{Float64,2} # spin dependent
    electrons::Electrons
    atoms::Atoms
    pspots::Array{PsPot_GTH,1}
    pspotNL::PsPotNL
    xcfunc::String
    ik::Int64   # current kpoint index
    ispin::Int64 # current spin index
end

import Base: println
function println( Ham::Hamiltonian; header=true )
    if header
        @printf("\n")
        @printf("                                  -----------\n")
        @printf("                                  Hamiltonian\n")
        @printf("                                  -----------\n")
        @printf("\n")
    end
    println("xcfunc = ", Ham.xcfunc)
    println("")
    println(Ham.atoms)
    println(Ham.pw)
    println(Ham.pw.gvecw.kpoints)
    println(Ham.electrons)
    for isp = 1:Ham.atoms.Nspecies
        println(Ham.pspots[isp])
    end
end

function Hamiltonian( atoms::Atoms, pspfiles::Array{String,1},
                      ecutwfc::Float64 ;
                      Nspin = 1,
                      meshk = [1,1,1], shiftk = [0,0,0],
                      kpoints = nothing,
                      xcfunc = "VWN",
                      extra_states = 0 )

    # kpoints
    if kpoints == nothing
        kpoints = KPoints( atoms, meshk, shiftk )
    else
        @assert typeof(kpoints) == KPoints
    end

    # Initialize plane wave grids
    pw = PWGrid( ecutwfc, atoms.LatVecs, kpoints=kpoints )

    Nspecies = atoms.Nspecies
    if Nspecies != size(pspfiles)[1]
        @printf("ERROR length of pspfiles is not equal to %d\n", Nspecies)
        exit()
    end

    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r

    strf = calc_strfact( atoms, pw )

    #
    # Initialize pseudopotentials and local potentials
    #
    Vg = zeros(ComplexF64, Npoints)
    V_Ps_loc = zeros(Float64, Npoints)

    Pspots = Array{PsPot_GTH}(undef,Nspecies)

    for isp = 1:Nspecies
        Pspots[isp] = PsPot_GTH( pspfiles[isp] )
        psp = Pspots[isp]
        for ig = 1:Ng
            ip = idx_g2r[ig]
            Vg[ip] = strf[ig,isp] * eval_Vloc_G( psp, G2[ig] ) / CellVolume
        end
        #
        V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(pw, Vg) ) * Npoints
    end

    # other potential terms are set to zero
    V_Hartree = zeros( Float64, Npoints )
    V_XC = zeros( Float64, Npoints, Nspin )
    potentials = Potentials( V_Ps_loc, V_Hartree, V_XC )
    #
    energies = Energies()
    #
    rhoe = zeros( Float64, Npoints, Nspin )

    electrons = Electrons( atoms, Pspots, Nspin=Nspin, Nkpt=kpoints.Nkpt,
                           Nstates_empty=extra_states )

    # NL pseudopotentials
    pspotNL = PsPotNL( pw, atoms, Pspots, kpoints, check_norm=false )

    atoms.Zvals = get_Zvals( Pspots )

    ik = 1
    ispin = 1
    return Hamiltonian( pw, potentials, energies, rhoe,
                          electrons, atoms, Pspots, pspotNL, xcfunc, ik, ispin )
end


#
# No pspfiles given. Use Coulomb potential (all electrons)
#
function Hamiltonian( atoms::Atoms, ecutwfc::Float64;
                      Nspin = 1,
                      meshk = [1,1,1], shiftk=[0,0,0],    
                      xcfunc="VWN", extra_states=0 )
    
    # kpoints
    kpoints = KPoints( atoms, meshk, shiftk )

    # Initialize plane wave grids
    pw = PWGrid( ecutwfc, atoms.LatVecs, kpoints=kpoints )

    Nspecies = atoms.Nspecies

    Npoints = prod(pw.Ns)

    # dummy pspots
    Pspots = Array{PsPot_GTH}(undef,1)
    Pspots[1] = PsPot_GTH()

    #
    # Coulomb potential as local potential
    #
    Zvals = get_Zatoms(atoms)
    strf = calc_strfact( atoms, pw )
    V_Ps_loc = init_V_coulomb_G( pw, strf, Zvals )

    # other potentials are set to zero
    V_Hartree = zeros( Float64, Npoints )
    V_XC = zeros( Float64, Npoints, Nspin )
    potentials = Potentials( V_Ps_loc, V_Hartree, V_XC )
    #
    energies = Energies()

    Zatoms = get_Zatoms( atoms )
    atoms.Zvals = Zatoms
    
    # use Zatoms as Zvals
    electrons = Electrons( atoms, Zatoms, Nspin=Nspin, Nkpt=kpoints.Nkpt,
                           Nstates_empty=extra_states )

    rhoe = zeros(Float64,Npoints,Nspin)

    pspotNL = PsPotNL()
    
    ik = 1
    ispin = 1
    return Hamiltonian( pw, potentials, energies, rhoe,
                          electrons, atoms, Pspots, pspotNL, xcfunc, ik, ispin )
end


"""
Given rhoe in real space, update Ham.rhoe, Hartree and XC potentials.
"""
function update!(Ham::Hamiltonian, rhoe::Array{Float64,1})
    # assumption Nspin = 1
    Ham.rhoe[:,1] = rhoe
    Ham.potentials.Hartree = real( G_to_R( Ham.pw, Poisson_solve(Ham.pw, rhoe) ) )
    if Ham.xcfunc == "PBE"
        Ham.potentials.XC[:,1] = calc_Vxc_PBE( Ham.pw, rhoe )
    else  # VWN is the default
        Ham.potentials.XC[:,1] = calc_Vxc_VWN( rhoe )
    end
    return
end

function update!(Ham::Hamiltonian, rhoe::Array{Float64,2})
    Nspin = size(rhoe)[2]
    if Nspin == 1
        update!(Ham, rhoe[:,1])
        return
    end
    Ham.rhoe = rhoe[:,:]
    Rhoe_total = Ham.rhoe[:,1] + Ham.rhoe[:,2] # Nspin is 2
    Ham.potentials.Hartree = real( G_to_R( Ham.pw, Poisson_solve(Ham.pw, Rhoe_total) ) )
    if Ham.xcfunc == "PBE"
        Ham.potentials.XC = calc_Vxc_PBE( Ham.pw, rhoe )
    else  # VWN is the default
        Ham.potentials.XC = calc_Vxc_VWN( rhoe )
    end
    return
end


function update!( Ham::Hamiltonian, atoms::Atoms,
    strf::Array{ComplexF64,2}, pspfiles::Array{String,1} )

    Nspecies = atoms.Nspecies
    if Nspecies != size(pspfiles)[1]
        @printf("ERROR length of pspfiles is not equal to %d\n", Nspecies)
        exit()
    end

    pw = Ham.pw

    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    G2 = pw.gvec.G2

    Vg = zeros(ComplexF64, Npoints)
    V_Ps_loc = zeros(Float64, Npoints)

    pw = Ham.pw
    for isp = 1:Nspecies
        psp = PsPot_GTH( pspfiles[isp] )
        println(psp)
        for ig = 1:Npoints
            Vg[ig] = strf[ig,isp] * eval_Vloc_G( psp, G2[ig], CellVolume )
        end
        #
        V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(pw, Vg) ) * Npoints
    end

    # Update
    Ham.potentials.Ps_loc = V_Ps_loc
    return

end
