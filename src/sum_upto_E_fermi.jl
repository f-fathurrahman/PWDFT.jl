function sum_upto_E_fermi(Focc, evals, E_fermi)
    sFocc = 0.0
    Nstates = size(Focc)[1]
    for ist = 1:Nstates
        if evals[ist] <= E_fermi
            sFocc = sFocc + Focc[ist]
        end
    end
    return sFocc
end
