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
    Nspin = pwinput.nspin
    meshk = pwinput.meshk
    if pwinput.nbnd != -1
        Ham = Hamiltonian(atoms, pspfiles, ecutwfc,
            meshk=[meshk[1], meshk[2], meshk[3]], dual=dual, Nstates=pwinput.nbnd, xcfunc=xcfunc, Ns_=Ns,
            Nspin=Nspin)
    else
        Ham = Hamiltonian(atoms, pspfiles, ecutwfc,
            meshk=[meshk[1], meshk[2], meshk[3]], dual=dual, xcfunc=xcfunc, Ns_=Ns,
            Nspin=Nspin)
    end

    return Ham, pwinput
end