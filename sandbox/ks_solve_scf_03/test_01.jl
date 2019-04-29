using Printf
using LinearAlgebra
using Random
using QuadGK
using SpecialFunctions

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("init_aux_loc.jl")

function init_Ham_Si_fcc()

    LATCONST = 10.2631

    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(LATCONST))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
end


function main()
    Ham = init_Ham_Si_fcc()
    println(Ham)

    V_Ps_loc, Rhoe_aux, E_alphat, E_sisa = init_aux_loc( Ham )
    
    V_Ps_loc_G = R_to_G(Ham.pw, V_Ps_loc)
    Rhoe_aux_G = R_to_G(Ham.pw, Rhoe_aux)

    V_aux = real( G_to_R( Ham.pw, Poisson_solve(Ham.pw, Rhoe_aux) ) )

    println(sum(Ham.potentials.Ps_loc))
    println(sum(V_Ps_loc))

    pw = Ham.pw

    Ng = pw.gvec.Ng
    Npoints = prod(pw.Ns)
    idx_g2r = pw.gvec.idx_g2r
    G2 = pw.gvec.G2
    CellVolume = pw.CellVolume
    dVol = CellVolume/Npoints

    Ngw = pw.gvecw.Ngw
    idx_gw2g = pw.gvecw.idx_gw2g
    G = pw.gvec.G

    electrons = Ham.electrons
    Nelectrons = electrons.Nelectrons
    Focc = copy(electrons.Focc) # make sure to use the copy
    Nstates = electrons.Nstates
    Nstates_occ = electrons.Nstates_occ
    Nspin = electrons.Nspin
    
    kpoints = pw.gvecw.kpoints
    Nkpt = kpoints.Nkpt
    k = kpoints.k
    wk = kpoints.wk

    Nkspin = Nkpt*Nspin


    u_Ps_loc_G = zeros(ComplexF64,Npoints)
    for ig = 2:Ng
        ip = idx_g2r[ig]
        u_Ps_loc_G[ip] = V_Ps_loc_G[ip] - 4*pi*Rhoe_aux_G[ip]/G2[ig]
    end
    u_Ps_loc = real( G_to_R(Ham.pw, u_Ps_loc_G) )

    Ham.potentials.Ps_loc = u_Ps_loc
    #Ham.potentials.Ps_loc = V_Ps_loc - V_aux


    Random.seed!(1234)

    psiks = rand_BlochWavefunc( Ham )
    Rhoe = calc_rhoe(Ham, psiks)
    
    Rhoe_tot = Rhoe[:,1]
    Rhoe_neu = Rhoe_tot + Rhoe_aux
    
    println("integ Rhoe     = ", sum(Rhoe)*dVol)
    println("integ Rhoe_neu = ", sum(Rhoe_neu)*dVol)

    @assert( Nspin == 1 )

    Ham.potentials.Hartree = real( G_to_R(pw, Poisson_solve(pw, Rhoe_neu)) )
    if Ham.xcfunc == "PBE"
        Ham.potentials.XC[:,1] = calc_Vxc_PBE( Ham.pw, Rhoe[:,1] )
    else  # VWN is the default
        Ham.potentials.XC[:,1] = calc_Vxc_VWN( Rhoe[:,1] )
    end


    Etot_old = 0.0

    Rhoe_new = zeros(Float64,Npoints,Nspin)

    diffRhoe = zeros(Nspin)

    evals = zeros(Float64,Nstates,Nkspin)

    ETHR_EVALS_LAST = 1e-6

    ethr = 0.1

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    # calculate PspCore energy
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    CONVERGED = 0

    betamix = 0.5

    ETOT_CONV_THR = 1e-6


    for iterSCF = 1:30

        evals = diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false,
                          Nstates_conv=Nstates_occ )

        Rhoe_new[:,:] = calc_rhoe( Nelectrons, pw, Focc, psiks, Nspin, renormalize=true )
        for ispin = 1:Nspin
            diffRhoe[ispin] = norm(Rhoe_new[:,ispin] - Rhoe[:,ispin])
        end

        for ispin = 1:Nspin
            Rhoe[:,ispin] = betamix*Rhoe_new[:,ispin] + (1-betamix)*Rhoe[:,ispin]
        end

        # Update potentials
        Ham.rhoe = Rhoe[:,:]
        Rhoe_tot = Rhoe[:,1]
        Rhoe_neu = Rhoe_tot + Rhoe_aux
                
        V_HartreeG = Poisson_solve(pw, Rhoe_neu)
        Ham.potentials.Hartree = real( G_to_R(pw, V_HartreeG) )

        if Ham.xcfunc == "PBE"
            Ham.potentials.XC[:,1] = calc_Vxc_PBE( Ham.pw, Rhoe[:,1] )
        else  # VWN is the default
            Ham.potentials.XC[:,1] = calc_Vxc_VWN( Rhoe[:,1] )
        end

       #
        # Kinetic energy
        #
        E_kin = 0.0
        for ispin = 1:Nspin
        for ik = 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin - 1)*Nkpt
            psi = psiks[ikspin]
            for ist = 1:Nstates
                psiKpsi = 0.0
                for igk = 1:Ngw[ik]
                    ig = idx_gw2g[ik][igk]
                    Gw2 = (G[1,ig] + k[1,ik])^2 + (G[2,ig] + k[2,ik])^2 + (G[3,ig] + k[3,ik])^2
                    psiKpsi = psiKpsi + abs(psi[igk,ist])^2*Gw2
                end
                E_kin = E_kin + wk[ik]*Focc[ist,ikspin]*psiKpsi
            end
        end
        end
        E_kin = 0.5*E_kin

        cRhoeG = conj(R_to_G(pw, Rhoe_neu))/Npoints
        cRhoeG_tot = conj(R_to_G(pw, Rhoe_tot))/Npoints

        if Ham.xcfunc == "PBE"
            epsxc = calc_epsxc_PBE( Ham.pw, Ham.rhoe )
        else
            epsxc = calc_epsxc_VWN( Ham.rhoe )
        end
        E_xc = dot( epsxc, Rhoe_tot ) * dVol

        E_self = sum( V_aux .* Rhoe_aux )*dVol

        E_Hartree = 0.0
        E_Ps_loc = 0.0
        for ig = 2:pw.gvec.Ng
            ip = pw.gvec.idx_g2r[ig]
            E_Ps_loc = E_Ps_loc + real( u_Ps_loc_G[ip] * conj(cRhoeG_tot[ip]) )
            E_Hartree = E_Hartree + 4*pi/G2[ig] * real( V_HartreeG[ip] * conj(cRhoeG[ip]) )
        end
        E_Ps_loc = E_Ps_loc*dVol + E_alphat - E_sisa
        E_Hartree = 0.5*E_Hartree*dVol + E_self
        #E_Hartree = 0.5*sum( Ham.potentials.Hartree.*Rhoe_neu )*dVol

        if Ham.pspotNL.NbetaNL > 0
            E_Ps_nloc = calc_E_Ps_nloc( Ham, psiks )
        else
            E_Ps_nloc = 0.0
        end

        Ham.energies.Kinetic = E_kin
        Ham.energies.Ps_loc  = E_Ps_loc
        Ham.energies.Ps_nloc = E_Ps_nloc
        Ham.energies.Hartree = E_Hartree
        Ham.energies.XC      = E_xc

        Etot = sum(Ham.energies)
        diffE = abs( Etot - Etot_old )

        @printf("\nSCF: %8d %18.10f %18.10e %18.10e\n", iterSCF, Etot, diffE, diffRhoe[1] )
        @printf("integ Rhoe = %18.10f\n", sum(Rhoe)*dVol)
   
        println(Ham.energies)

        if diffE < ETOT_CONV_THR
            CONVERGED = CONVERGED + 1
        else  # reset CONVERGED
            CONVERGED = 0
        end

        if CONVERGED >= 2
            if verbose
                @printf("SCF is converged: iter: %d , diffE = %10.7e\n", iterSCF, diffE)
            end
            break
        end
        #
        Etot_old = Etot

        flush(stdout)



    end

end

main()