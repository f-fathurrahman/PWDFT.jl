function test_02( Ham )

    Nspecies = Ham.atoms.Nspecies
    is_psp_using_nlcc = zeros(Bool, Nspecies) # by default we don't use any NLCC    
    for isp in 1:Nspecies
        if Ham.pspots[isp].is_nlcc
            is_psp_using_nlcc[isp] = true
        end
    end

    if any(is_psp_using_nlcc)
        res = @be calc_rhoe_core!(Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe_core)
        display(res)
    else
        println("No rhoe_core is needed")
    end

    return
end
