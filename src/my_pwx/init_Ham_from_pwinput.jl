function init_Ham_from_pwinput(; filename::Union{Nothing,String}=nothing)

    # XXX simply assign filename=ARGS[1] is length(ARGS) >= 1

    # Read filename from the command line argument
    if isnothing(filename)
        println("ARGS = ", ARGS)
        if length(ARGS) == 0
            println("No argument is given. Using PWINPUT as default filename.")
            @assert Base.Filesystem.isfile("PWINPUT")
            pwinput = PWSCFInput("PWINPUT")
        else
            pwinput = PWSCFInput(ARGS[1])
        end
    else
        pwinput = PWSCFInput(filename)
    end

    atoms = pwinput.atoms

    ecutwfc = pwinput.ecutwfc
    ecutrho = pwinput.ecutrho
    dual = ecutrho/ecutwfc

    pspfiles = pwinput.pspfiles
    # Need special treatement for GTH ?
    for isp in 1:atoms.Nspecies
        if is_using_extension_gth(pspfiles[isp])
            error("GTH pspot is not yet supported")
        end
    end

    xcfunc = decide_xcfunc(pwinput)
    # Note that other cases will be defaulted to VWN
    # XXX probably need to process input_dft string
    println("Using xcfunc = ", xcfunc)

    Ns = (pwinput.nr1, pwinput.nr2, pwinput.nr3)    
    if pwinput.nspin == 2
        @assert !pwinput.noncolin
    end
    if pwinput.lspinorb
        @assert pwinput.noncolin
    end

    options = HamiltonianOptions()
    options.meshk[1] = pwinput.meshk[1]
    options.meshk[2] = pwinput.meshk[2]
    options.meshk[3] = pwinput.meshk[3] 
    options.dual = dual
    if pwinput.nbnd != -1
        options.Nstates = pwinput.nbnd
    end
    options.xcfunc = xcfunc
    options.Ns = Ns
    options.lspinorb = pwinput.lspinorb
    options.noncollinear = pwinput.noncolin
    #
    # Determine Nspin_channel and Nspin_comp from pwinput
    if pwinput.lspinorb
        options.Nspin_channel = 1
        options.Nspin_comp = 4
        options.time_reversal = false
        @assert pwinput.noncolin
    end
    if pwinput.noncolin
        options.Nspin_channel = 1
        options.Nspin_comp = 4
        options.time_reversal = false
    else
        options.Nspin_channel = pwinput.nspin
        options.Nspin_comp = pwinput.nspin
    end

    if !isnothing(pwinput.starting_magnetization)
        options.starting_magn = pwinput.starting_magnetization
    end
    if !isnothing(pwinput.angle1)
        options.angle1 = pwinput.angle1
    end
    if !isnothing(pwinput.angle2)
        options.angle2 = pwinput.angle2
    end


    pspots = Vector{PsPot_UPF}(undef, atoms.Nspecies)
    for isp in 1:atoms.Nspecies
        pspots[isp] = PsPot_UPF(pspfiles[isp])
    end

    Ham = Hamiltonian(atoms, pspots, ecutwfc, options)

#    if pwinput.nbnd != -1
#        # nbnd is given from pwinput
#        @info "Pass here 45 in init_Ham_from_pwinput"
#        Ham = Hamiltonian(atoms, pspfiles, ecutwfc,
#            meshk=[meshk[1], meshk[2], meshk[3]], dual=dual, Nstates=pwinput.nbnd, xcfunc=xcfunc, Ns_=Ns,
#            Nspin=Nspin, use_soc=lspinorb, use_noncol_magn=noncolin)
#    else
#        Ham = Hamiltonian(atoms, pspfiles, ecutwfc,
#            meshk=[meshk[1], meshk[2], meshk[3]], dual=dual, xcfunc=xcfunc, Ns_=Ns,
#            Nspin=Nspin, use_soc=lspinorb, use_noncol_magn=noncolin)
#    end

    return Ham, pwinput
end