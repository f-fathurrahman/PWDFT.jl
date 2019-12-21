push!(LOAD_PATH, "./")

using LinearAlgebra
using PWDFT2
using Random

include("calc_E_NN_mod.jl")

function gen_Rhoe_aux_G(
    atoms::Atoms, Zvals::Array{Float64,1}, pw::PWGrid; TOL = 1e-8
)

    # determine eta
    Gcut = 2*pw.ecutwfc/(2*pi)
    eta = 0.5*Gcut^2/-log(TOL)

    # structure factor
    Sf = calc_strfact( atoms, pw )
    
    Ng = pw.gvec.Ng
    G2 = pw.gvec.G2
    Nspecies = atoms.Nspecies
    Npoints = prod(pw.Ns)
    idx_g2r = pw.gvec.idx_g2r
    CellVolume = pw.CellVolume
    
    Rhoe_aux_G = zeros(ComplexF64,Npoints)
    Rhoe_aux = zeros(Float64,Npoints)

    for ig = 1:Ng
        ip = idx_g2r[ig]
        for isp = 1:Nspecies
            Rhoe_aux_G[ip] = Rhoe_aux_G[ip] + Zvals[isp]*exp(-0.125*G2[ig]/eta^2)*Sf[ig,isp]/CellVolume
        end
    end
    Rhoe_aux = real( G_to_R(pw,Rhoe_aux_G) )
    return -Rhoe_aux*Npoints
end

function create_Ham_Si()
    atoms = Atoms(xyz_string_frac=
        """
        2
        
        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631) )
    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Si-q4.gth"]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
end

function create_Ham_GaAs()
    atoms = Atoms(xyz_string_frac=
        """
        2
        
        Ga  0.0  0.0  0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.6537*ANG2BOHR) )
    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Ga-q3.gth",
                "../pseudopotentials/pade_gth/As-q5.gth"]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
end

function create_Ham_H()
    atoms = Atoms(xyz_string=
        """
        1
        
        H  0.0  0.0  0.0
        """, in_bohr=true, LatVecs=gen_lattice_sc(16.0) )
    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end    


function test_main()

    Ham = create_Ham_Si()
    #Ham = create_Ham_H()
    #Ham = create_Ham_GaAs()

    atoms = Ham.atoms
    println(atoms)

    pw = Ham.pw
    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    electrons = Ham.electrons
    Nspin = electrons.Nspin
    Nkpt = pw.gvecw.kpoints.Nkpt
    Focc = electrons.Focc
    Npoints = prod(pw.Ns)
    Nspecies = atoms.Nspecies

    dVol = CellVolume/Npoints
    #
    Nels = PWDFT2.get_Nelectrons( atoms, Ham.pspots )
    #
    Zvals = get_Zvals(Ham.pspots)
    Rhoe_aux = gen_Rhoe_aux_G(atoms, Zvals, pw)
    #
    println("Nelectrons = ", Nels)
    println("G space: integ Rhoe_aux = ", sum(Rhoe_aux)*dVol)

    #
    # Random guess of wave function
    #
    Random.seed!(1234)
    psiks = rand_BlochWavefunc(pw, electrons)
    Rhoe = zeros(Float64,Npoints,Nspin)
    for ispin = 1:Nspin
        idxset = (Nkpt*(ispin-1)+1):(Nkpt*ispin)
        Rhoe[:,ispin] = calc_rhoe( Nels, pw, Focc[:,idxset], psiks[idxset], Nspin )
    end

    @assert Nspin == 1  # restrict to non spin-polarized

    # Hartree potential (real space, using plain Rhoe)
    Rhoe_ = Rhoe[:,1]
    VHartree = real( G_to_R( pw, Poisson_solve(pw, Rhoe_) ) )
    Ehartree = 0.5*dot(VHartree,Rhoe_)*dVol
    println("Ehartree = ", Ehartree)

    # Local ionic potential
    E_ps_loc = dot( Ham.potentials.Ps_loc, Rhoe_ )*dVol
    println("E_ps_loc = ", E_ps_loc)

    E_NN = calc_E_NN(atoms)
    println("E_NN = ", E_NN)

    E_pspcore = calc_PspCore_ene( Ham.atoms, Ham.pspots )
    println("E_pspcore = ", E_pspcore)

    println("Electrostatic ene (excluding pspcore) = ", Ehartree + E_ps_loc + E_NN)
    println("Electrostatic ene (including pspcore) = ", Ehartree + E_ps_loc + E_NN + E_pspcore)

    println("\nUsing G-space formula")

    Ng = pw.gvec.Ng
    G = pw.gvec.G
    G2 = pw.gvec.G2
    idx_g2r = pw.gvec.idx_g2r

    # Need to normalize by 1/Npoints 
    RhoeG = R_to_G(pw, Rhoe_)/Npoints
    #
    Rhoe_aux_G = R_to_G(pw, Rhoe_aux)/Npoints
    #
    RhoeG_T = RhoeG + Rhoe_aux_G
    #
    # Calculate Hartree potential from this total chg + gaussian charge
    VHartreeG = zeros(ComplexF64,Npoints)
    VHartreeG[1] = 0.0 + im*0.0
    for ig = 2:Ng
        ip = idx_g2r[ig]
        VHartreeG[ip] = RhoeG_T[ig]/G2[ig]
    end
    VHartreeG = 4*pi*VHartreeG
    VHartree2 = real( G_to_R( pw, VHartreeG ) )
    #
    pspots = Ham.pspots
    Vg = zeros(ComplexF64,Npoints)
    V_Ps_loc = zeros(Float64,Npoints)
    #
    TOL  = 1e-8
    Gcut = 2*pw.ecutwfc/(2*pi)
    gexp = -log(TOL)
    eta  = 0.5*Gcut^2/gexp
    #
    strf = calc_strfact(atoms, pw)
    #
    # Modified V_ps_loc (plus potential correspond to by Gaussian chg)
    #
    for isp = 1:Nspecies
        psp = pspots[isp]
        for ig = 2:Ng
            ip = idx_g2r[ig]
            Vg[ip] = strf[ig,isp] * ( eval_Vloc_G( psp, G2[ig] ) + 
                     4*pi/G2[ig] * Zvals[isp] * exp(-0.125*G2[ig]/eta^2 )
                     )/ CellVolume
        end
        #
        V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(pw, Vg) ) * Npoints
    end
    V_Ps_locG = R_to_G(pw, V_Ps_loc)/Npoints
    #
    EhartreeG = 0.0
    E_ps_locG = 0.0
    #
    for ig = 2:Ng
        ip = idx_g2r[ig]
        #EhartreeG = EhartreeG + 0.5*real(VHartreeG[ip]*conj(RhoeG[ip]))*CellVolume
        EhartreeG = EhartreeG + 2*pi/G2[ig] * real(RhoeG_T[ip]*conj(RhoeG_T[ip]))*CellVolume  # using chgden + gaussian
        #
        E_ps_locG = E_ps_locG + real(V_Ps_locG[ip]*conj(RhoeG[ip]))*CellVolume  # using only the chgden
    end
    
    println("atoms.Zvals = ", atoms.Zvals)
    
    #E_NN_2 = calc_E_NN_mod(pw, atoms)
    E_NN_2 = calc_E_NN(atoms)  # E_NN_2 contains self energy of Gaussian chgden

    println("EhartreeG = ", EhartreeG)
    println("E_ps_locG = ", E_ps_locG)
    println("E_NN_2    = ", E_NN_2)

    println("E ps loc by real space integ = ", sum(Ham.potentials.Ps_loc.*Rhoe_)*dVol)

    println("E Electrostatic v2 = ", EhartreeG + E_ps_locG + E_NN_2 + E_pspcore) # not correct ? gaussian chg contrib is counted two times

    println("E Electrostatic v3 = ", EhartreeG + sum(Ham.potentials.Ps_loc.*Rhoe_)*dVol + E_NN_2 + E_pspcore)


    # Comparing V_ps_loc + V_Hartree
    Vloc1 = Ham.potentials.Ps_loc + VHartree
    Vloc2 = V_Ps_loc + VHartree2

    println("diff Vloc = ", sum(Vloc1 - Vloc2)/Npoints)

end

test_main()