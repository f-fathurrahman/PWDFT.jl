function init_Ham_from_pwinput(; filename::Union{Nothing,String}=nothing)

    if filename == nothing
        println("ARGS = ", ARGS)
        @assert length(ARGS) == 1
        pwinput = PWSCFInput(ARGS[1])
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

    xcfunc = "VWN" # default in PWDFT.jl
    if uppercase(pwinput.input_dft) == "SCAN"
        xcfunc = "SCAN"
    elseif uppercase(pwinput.input_dft) in ["PBE", "GGA_X_PBE+GGA_C_PBE"]
        xcfunc = "PBE"
    end
    # Note that other cases will be defaulted to VWN
    # XXX probably need to process input_dft string
    println("Using xcfunc = ", xcfunc)

    Ns = (pwinput.nr1, pwinput.nr2, pwinput.nr3)
    meshk = pwinput.meshk
    if pwinput.nbnd != -1
        Ham = Hamiltonian(atoms, pspfiles, ecutwfc,
            meshk=[meshk[1], meshk[2], meshk[3]], dual=dual, Nstates=pwinput.nbnd, xcfunc=xcfunc, Ns_=Ns)
    else
        Ham = Hamiltonian(atoms, pspfiles, ecutwfc,
            meshk=[meshk[1], meshk[2], meshk[3]], dual=dual, xcfunc=xcfunc, Ns_=Ns)
    end

    return Ham, pwinput
end