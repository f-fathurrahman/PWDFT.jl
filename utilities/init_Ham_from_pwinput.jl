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
    end

    meshk = pwinput.meshk
    if pwinput.nbnd != -1
        Ham = Hamiltonian(atoms, pspfiles, ecutwfc,
            meshk=[meshk[1], meshk[2], meshk[3]], dual=dual, Nstates=pwinput.nbnd, xcfunc=xcfunc)
    else
        Ham = Hamiltonian(atoms, pspfiles, ecutwfc,
            meshk=[meshk[1], meshk[2], meshk[3]], dual=dual, xcfunc=xcfunc)
    end

    return Ham, pwinput
end