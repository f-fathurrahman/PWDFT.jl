using Printf
using Random
using PWDFT

function investigate( Ham, psiks )

    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    Natoms = Ham.atoms.Natoms
    atm2species = Ham.atoms.atm2species
    Pspots = Ham.pspots
    prj2beta = Ham.pspotNL.prj2beta
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    wk = Ham.pw.gvecw.kpoints.wk
    NbetaNL = Ham.pspotNL.NbetaNL
    Nspin = Ham.electrons.Nspin

    for ispin = 1:Nspin, ik = 1:Nkpt
        
        ikspin = ik + (ispin - 1)*Nkpt
        psi = psiks[ikspin]
        #betaNL_psi = calc_betaNL_psi( ik, Ham.pspotNL.betaNL, psi )
        
        for ist = 1:Nstates

            @printf("\nStates = %d\n", ist)

            idx_enl = 0
            for ia = 1:Natoms
                @printf("\n")
                isp = atm2species[ia]
                psp = Pspots[isp]
                for l = 0:psp.lmax
                for m = -l:l
                for iprj = 1:psp.Nproj_l[l+1]
                for jprj = 1:psp.Nproj_l[l+1]
                    ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                    jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                    hij = psp.h[l+1,iprj,jprj]
                    #enl1 = enl1 + hij*real(conj(betaNL_psi[ist,ibeta])*betaNL_psi[ist,jbeta])
                    idx_enl = idx_enl + 1
                    @printf("%5d [%5d] [%2d,%2d] %5d %5d %18.10f\n", idx_enl, ia, l, m, ibeta, jbeta, hij)
                end
                end
                end # m
                end # l
            end
            #E_Ps_nloc = E_Ps_nloc + wk[ik]*Focc[ist,ikspin]*enl1
        end

    end

end


function main()

    Random.seed!(1234)

    atoms = Atoms(xyz_string=
    """
    1

    Pt   0.0  0.0  0.0
    """, LatVecs=gen_lattice_sc(10.0))

    pspfiles = ["../../pseudopotentials/pade_gth/Pt-q10.gth"]
    ecutwfc = 15.0

    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
    psiks = rand_BlochWavefunc( Ham )

    println("")
    println("NbetaNL = ", Ham.pspotNL.NbetaNL)
    println("")
    investigate( Ham, psiks )

    for isp = 1:Ham.atoms.Nspecies
        println( Ham.pspots[1] )
    end

end

main()
