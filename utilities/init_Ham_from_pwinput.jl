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

    meshk = pwinput.meshk

    return Hamiltonian(atoms, pspfiles, ecutwfc, meshk=[meshk[1], meshk[2], meshk[3]], dual=dual)
end