using LinearAlgebra
using Random
using Printf

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures", "DATA_DeltaCodes", "PWINPUT")

function get_default_psp_uspp_pslib(atoms::Atoms; xcfunc="VWN")

    if xcfunc == "PBE"
        dir_psp = joinpath(DIR_PWDFT, "pseudopotentials", "PSLIB_US_PAW_PBE")
    elseif xcfunc == "VWN"
        dir_psp = joinpath(DIR_PWDFT, "pseudopotentials", "PSLIB_US_PAW_LDA")
    else
        error("Unsupported xcfunc = $xcfunc")
    end

    #XXX These can be cached ?
    listfiles = readdir(dir_psp)
    atsymb_list = String[]
    for l in listfiles
        atsymb = split(l, ".")[1]
        if !(atsymb in atsymb_list)
            append!(atsymb_list, [atsymb])
        end
    end
    dict_uspp = Dict{String,Vector{String}}()
    for atsymb in atsymb_list
        list_psp_atsymb = String[]
        for l in listfiles
            if occursin("rrkjus", l) && atsymb==split(l, ".")[1]
                append!(list_psp_atsymb, [l])
            end
        end
        merge!(dict_uspp, Dict(atsymb => list_psp_atsymb))
    end

    Nspecies = atoms.Nspecies
    pspfiles = Array{String}(undef,Nspecies)
    SpeciesSymbols = atoms.SpeciesSymbols
    for isp in 1:Nspecies
        atsymb = SpeciesSymbols[isp]
        # choose the one with the shortest filename
        idx_chosen = argmin(length.(dict_uspp[atsymb]))
        pspfiles[isp] = joinpath(dir_psp, dict_uspp[atsymb][idx_chosen])
    end
    return pspfiles
end



function do_calc_uspp_pslib( atom_symbol::String; xcfunc="VWN" )

    Random.seed!(1234)

    pwinput = PWSCFInput(joinpath(DIR_STRUCTURES, atom_symbol*".in"))
    atoms = pwinput.atoms

    # Initialize Hamiltonian
    pspfiles = get_default_psp_uspp_pslib(atoms, xcfunc=xcfunc)
    ecutwfc = 20.0
    dual = 5.0
    Ham = Hamiltonian(
        atoms, pspfiles, ecutwfc,
        dual = dual,
        xcfunc = xcfunc,
        use_smearing = true,
        smearing_kT = 0.01,
        meshk = collect(pwinput.meshk)
    )
    println(Ham)
    #
    psiks = rand_BlochWavefunc(Ham)
    electrons_scf_G!(Ham, psiks = psiks, betamix = 0.1, print_final_ebands = false)
    #
    # Additional info
    println("Final Focc = ")
    display(Ham.electrons.Focc); println()
    println("Final ebands (w.r.t) Fermi energy = ")
    display(Ham.electrons.ebands .- Ham.electrons.E_fermi); println()
    return
end


function do_calc_uspp_pslib_Emin( atom_symbol::String; xcfunc="VWN", use_initial_scf = true )

    Random.seed!(1234)

    pwinput = PWSCFInput(joinpath(DIR_STRUCTURES, atom_symbol*".in"))
    atoms = pwinput.atoms

    # Initialize Hamiltonian
    pspfiles = get_default_psp_uspp_pslib(atoms, xcfunc=xcfunc)
    ecutwfc = 20.0
    dual = 5.0
    Ham = Hamiltonian(
        atoms, pspfiles, ecutwfc,
        dual = dual,
        xcfunc = xcfunc,
        use_smearing = true,
        smearing_kT = 0.01,
        meshk = collect(pwinput.meshk),
    )
    println(Ham)
    # Prepare Haux
    Nstates = Ham.electrons.Nstates;
    Nspin = Ham.electrons.Nspin_wf;
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt;
    Nkspin = Nkpt*Nspin;
    Haux = Vector{Matrix{ComplexF64}}(undef, Nkspin);
    for ikspin in 1:Nkspin
        Haux[ikspin] = randn(ComplexF64, Nstates, Nstates);
        Haux[ikspin][:,:] = 0.5*( Haux[ikspin] + Haux[ikspin]' );
    end

    if use_initial_scf
        psiks = rand_BlochWavefunc(Ham)
        electrons_scf_G!(
            Ham,
            psiks = psiks,
            NiterMax = 5,
            betamix = 0.1,
            starting_magn = Ham.options.starting_magn,
            print_final_ebands = false
        )
        # Rewrite Haux
        for ikspin in 1:Nkspin
            Haux[ikspin][:,:] = diagm(0=>Ham.electrons.ebands[:,ikspin])
        end
    else
        psiks = PWDFT.rand_wfc(Ham);
        #psiks = zeros_BlochWavefunc(Ham);
        #initwfc!(Ham, psiks); # This not yet ready
        #
        #psiks = rand_BlochWavefunc(Ham);
    end

    Rhoe = calc_rhoe(Ham, psiks)

    electrons_Emin_Haux!(Ham, psiks=psiks, Haux=Haux, Rhoe=Rhoe, NiterMax=5)
    electrons_Emin_Haux!(Ham, psiks=psiks, Haux=Haux, Rhoe=Rhoe, NiterMax=100)

    return
end
